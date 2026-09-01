"""Reproducible point-based spatial analyses for chemokine/CAR-T studies.

The module operates on a long-format table with one row per segmented cell or
landmark point.  Calculations are performed independently within each patient
and tissue section so that spatial coordinates are never compared across
unrelated specimens.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse, FancyArrowPatch, PathPatch
from matplotlib.path import Path
from scipy.spatial import cKDTree

REQUIRED_SPATIAL_COLUMNS = (
    "patient_id",
    "section_id",
    "cell_id",
    "x",
    "y",
    "cell_type",
)

DEFAULT_QUERY_COLUMNS = {
    "ligand_source": "is_ligand_source",
    "car_t": "is_car_t",
}

DEFAULT_LANDMARKS = ("vessel", "tumor_nest", "stroma", "lymphatic")


def _as_column_tuple(columns: Iterable[str]) -> tuple[str, ...]:
    result = tuple(columns)
    if not result or any(not isinstance(column, str) or not column for column in result):
        raise ValueError("Column names must be non-empty strings.")
    if len(set(result)) != len(result):
        raise ValueError("Column names must be unique.")
    return result


def _coerce_boolean(series: pd.Series, column: str) -> pd.Series:
    """Return a strict nullable-free boolean series."""

    if series.isna().any():
        raise ValueError(f"Column '{column}' contains missing values.")
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.astype(bool)

    mapping = {
        True: True,
        False: False,
        "1": True,
        "0": False,
        "true": True,
        "false": False,
        "yes": True,
        "no": False,
    }

    def convert(value: Any) -> bool | None:
        normalized = value.strip().lower() if isinstance(value, str) else value
        return mapping.get(normalized)

    converted = series.map(convert)
    if converted.isna().any():
        examples = series.loc[converted.isna()].astype(str).drop_duplicates().head(5)
        raise ValueError(
            f"Column '{column}' must contain boolean-like values; invalid values: "
            + ", ".join(examples)
        )
    return converted.astype(bool)


def validate_spatial_table(
    table: pd.DataFrame,
    *,
    required_columns: Iterable[str] = REQUIRED_SPATIAL_COLUMNS,
    patient_col: str = "patient_id",
    section_col: str = "section_id",
    cell_id_col: str = "cell_id",
    coordinate_cols: Sequence[str] = ("x", "y"),
    boolean_columns: Iterable[str] = (),
    copy: bool = True,
) -> pd.DataFrame:
    """Validate and normalize a long-format cell/landmark table.

    Parameters
    ----------
    table:
        One row per segmented cell or landmark point.
    required_columns:
        Columns that must be present.  The defaults define the canonical input
        schema.
    patient_col, section_col, cell_id_col:
        Keys identifying a unique observation.  Cell identifiers only need to
        be unique within a patient and section.
    coordinate_cols:
        Two or more numeric coordinate columns.  The distance functions use
        the first two by default but can also operate in three dimensions.
    boolean_columns:
        Optional flag columns to coerce from strict boolean-like values.
    copy:
        If true, return a defensive copy.

    Returns
    -------
    pandas.DataFrame
        A validated table with numeric coordinates and normalized flag columns.
    """

    if not isinstance(table, pd.DataFrame):
        raise TypeError("table must be a pandas DataFrame.")
    if table.empty:
        raise ValueError("Spatial table is empty.")

    required = set(_as_column_tuple(required_columns))
    coordinates = _as_column_tuple(coordinate_cols)
    if len(coordinates) < 2:
        raise ValueError("At least two coordinate columns are required.")
    required.update((patient_col, section_col, cell_id_col, *coordinates))
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))

    result = table.copy(deep=True) if copy else table
    key_columns = (patient_col, section_col, cell_id_col)
    for column in key_columns:
        if result[column].isna().any():
            raise ValueError(f"Key column '{column}' contains missing values.")
        blank = result[column].astype(str).str.strip().eq("")
        if blank.any():
            raise ValueError(f"Key column '{column}' contains blank values.")

    descriptor_columns = required.difference((*key_columns, *coordinates))
    for column in sorted(descriptor_columns):
        if result[column].isna().any():
            raise ValueError(f"Required descriptor column '{column}' contains missing values.")
        if pd.api.types.is_object_dtype(result[column].dtype) or isinstance(
            result[column].dtype, pd.StringDtype
        ):
            if result[column].astype(str).str.strip().eq("").any():
                raise ValueError(f"Required descriptor column '{column}' contains blank values.")

    duplicate_mask = result.duplicated(list(key_columns), keep=False)
    if duplicate_mask.any():
        duplicate = result.loc[duplicate_mask, list(key_columns)].iloc[0].to_dict()
        raise ValueError(f"Duplicate patient/section/cell key: {duplicate}")

    for column in coordinates:
        numeric = pd.to_numeric(result[column], errors="coerce")
        invalid = numeric.isna() | ~np.isfinite(numeric.to_numpy(dtype=float))
        if invalid.any():
            raise ValueError(f"Coordinate column '{column}' contains non-finite values.")
        result[column] = numeric.astype(float)

    boolean_flags = tuple(boolean_columns)
    if boolean_flags:
        boolean_flags = _as_column_tuple(boolean_flags)
    for column in boolean_flags:
        if column not in result.columns:
            raise ValueError(f"Missing boolean flag column: {column}")
        result[column] = _coerce_boolean(result[column], column)

    return result


def _group_key(key: Any) -> tuple[Any, ...]:
    return key if isinstance(key, tuple) else (key,)


def distance_to_landmarks(
    table: pd.DataFrame,
    *,
    landmark_table: pd.DataFrame | None = None,
    query_columns: Mapping[str, str] | None = None,
    structure_col: str = "structure",
    landmark_values: Sequence[str] = DEFAULT_LANDMARKS,
    patient_col: str = "patient_id",
    section_col: str = "section_id",
    cell_id_col: str = "cell_id",
    coordinate_cols: Sequence[str] = ("x", "y"),
    coordinate_scale_um: float = 1.0,
) -> pd.DataFrame:
    """Measure nearest-point distances from ligand-source/CAR-T cells to landmarks.

    Query flags and landmark points may be stored in the same table.  A separate
    ``landmark_table`` is useful when vascular, stromal, tumor-nest or lymphatic
    boundaries have been sampled from segmentation masks.  Distances are
    calculated only within the same patient and section.

    ``coordinate_scale_um`` converts coordinate units to micrometres.  Point-to-
    point distances approximate distance to a structure only when the landmark
    table densely samples the segmented boundary or region.
    """

    if query_columns is None:
        query_columns = DEFAULT_QUERY_COLUMNS
    if not isinstance(query_columns, Mapping) or not query_columns:
        raise ValueError("query_columns must map at least one role to a flag column.")
    if any(not isinstance(role, str) or not role.strip() for role in query_columns):
        raise ValueError("Query role names must be non-empty strings.")
    if len(set(query_columns)) != len(query_columns):
        raise ValueError("Query role names must be unique.")
    flag_columns = tuple(query_columns.values())
    if len(set(flag_columns)) != len(flag_columns):
        raise ValueError("Each query role must use a distinct flag column.")
    if not np.isfinite(coordinate_scale_um) or coordinate_scale_um <= 0:
        raise ValueError("coordinate_scale_um must be a finite positive number.")

    coordinates = _as_column_tuple(coordinate_cols)
    cells = validate_spatial_table(
        table,
        patient_col=patient_col,
        section_col=section_col,
        cell_id_col=cell_id_col,
        coordinate_cols=coordinates,
        boolean_columns=flag_columns,
    )
    landmarks = (
        cells
        if landmark_table is None
        else validate_spatial_table(
            landmark_table,
            patient_col=patient_col,
            section_col=section_col,
            cell_id_col=cell_id_col,
            coordinate_cols=coordinates,
        )
    )
    if structure_col not in landmarks.columns:
        raise ValueError(f"Missing landmark classification column: {structure_col}")

    landmark_names = tuple(str(value) for value in landmark_values)
    if not landmark_names or len(set(landmark_names)) != len(landmark_names):
        raise ValueError("landmark_values must contain unique values.")

    group_columns = [patient_col, section_col]
    landmark_groups = {
        _group_key(key): group
        for key, group in landmarks.groupby(group_columns, sort=True, observed=True)
    }
    records: list[dict[str, Any]] = []

    for role, flag_col in query_columns.items():
        query_subset = cells.loc[cells[flag_col]].copy()
        for key, queries in query_subset.groupby(group_columns, sort=True, observed=True):
            normalized_key = _group_key(key)
            local_landmarks = landmark_groups.get(normalized_key)
            query_coordinates = queries.loc[:, coordinates].to_numpy(dtype=float)

            for landmark_name in landmark_names:
                if local_landmarks is None:
                    reference = landmarks.iloc[0:0]
                else:
                    reference = local_landmarks.loc[
                        local_landmarks[structure_col].astype(str).eq(landmark_name)
                    ]

                if reference.empty:
                    distances = np.full(len(queries), np.nan, dtype=float)
                    nearest_ids: list[Any] = [pd.NA] * len(queries)
                else:
                    reference_coordinates = reference.loc[:, coordinates].to_numpy(dtype=float)
                    tree = cKDTree(reference_coordinates)
                    raw_distances, nearest_indices = tree.query(query_coordinates, k=1)
                    distances = np.asarray(raw_distances, dtype=float) * coordinate_scale_um
                    reference_ids = reference[cell_id_col].to_numpy(dtype=object)
                    nearest_ids = reference_ids[np.asarray(nearest_indices, dtype=int)].tolist()

                for row_position, (_, query) in enumerate(queries.iterrows()):
                    record = {
                        patient_col: query[patient_col],
                        section_col: query[section_col],
                        cell_id_col: query[cell_id_col],
                        "query_role": str(role),
                        "landmark": landmark_name,
                        "distance_um": float(distances[row_position]),
                        "nearest_landmark_id": nearest_ids[row_position],
                        "n_landmarks": int(len(reference)),
                    }
                    for coordinate in coordinates:
                        record[coordinate] = float(query[coordinate])
                    records.append(record)

    output_columns = [
        patient_col,
        section_col,
        cell_id_col,
        "query_role",
        "landmark",
        "distance_um",
        "nearest_landmark_id",
        "n_landmarks",
        *coordinates,
    ]
    if not records:
        return pd.DataFrame(columns=output_columns)
    return (
        pd.DataFrame.from_records(records, columns=output_columns)
        .sort_values(
            [patient_col, section_col, "query_role", cell_id_col, "landmark"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )


def summarize_patient_distances(
    distances: pd.DataFrame,
    *,
    patient_col: str = "patient_id",
    section_col: str = "section_id",
    role_col: str = "query_role",
    landmark_col: str = "landmark",
    distance_col: str = "distance_um",
) -> pd.DataFrame:
    """Summarize distances per patient after section-level aggregation.

    The primary estimate, ``mean_section_median_distance_um``, gives each tissue
    section equal weight within a patient and therefore avoids allowing a large
    section to dominate solely because it contains more segmented cells.
    Pooled quantiles are retained as descriptive distribution summaries.
    """

    if not isinstance(distances, pd.DataFrame):
        raise TypeError("distances must be a pandas DataFrame.")
    required = {patient_col, section_col, role_col, landmark_col, distance_col}
    missing = sorted(required.difference(distances.columns))
    if missing:
        raise ValueError("Missing distance columns: " + ", ".join(missing))

    working = distances.copy()
    working[distance_col] = pd.to_numeric(working[distance_col], errors="coerce")
    finite = np.isfinite(working[distance_col].to_numpy(dtype=float))
    working = working.loc[finite]
    output_columns = [
        patient_col,
        role_col,
        landmark_col,
        "n_sections",
        "n_query_points",
        "pooled_median_distance_um",
        "pooled_q25_distance_um",
        "pooled_q75_distance_um",
        "mean_section_median_distance_um",
        "sd_section_median_distance_um",
    ]
    if working.empty:
        return pd.DataFrame(columns=output_columns)

    section_keys = [patient_col, section_col, role_col, landmark_col]
    section_summary = (
        working.groupby(section_keys, sort=True, observed=True)[distance_col]
        .agg(section_median="median", n_query_points="size")
        .reset_index()
    )
    patient_keys = [patient_col, role_col, landmark_col]
    hierarchical = (
        section_summary.groupby(patient_keys, sort=True, observed=True)
        .agg(
            n_sections=(section_col, "nunique"),
            n_query_points=("n_query_points", "sum"),
            mean_section_median_distance_um=("section_median", "mean"),
            sd_section_median_distance_um=("section_median", "std"),
        )
        .reset_index()
    )
    pooled = (
        working.groupby(patient_keys, sort=True, observed=True)[distance_col]
        .agg(
            pooled_median_distance_um="median",
            pooled_q25_distance_um=lambda values: values.quantile(0.25),
            pooled_q75_distance_um=lambda values: values.quantile(0.75),
        )
        .reset_index()
    )
    return (
        hierarchical.merge(pooled, on=patient_keys, validate="one_to_one")
        .loc[:, output_columns]
        .sort_values(patient_keys, kind="mergesort")
        .reset_index(drop=True)
    )


def _hierarchical_group_difference(
    values: np.ndarray,
    group_a_mask: np.ndarray,
    section_indices: Sequence[np.ndarray],
    section_patients: Sequence[Any],
) -> float:
    patient_differences: dict[Any, list[float]] = {}
    for indices, patient in zip(section_indices, section_patients, strict=True):
        local_groups = group_a_mask[indices]
        local_values = values[indices]
        if not local_groups.any() or local_groups.all():
            continue
        difference = float(local_values[local_groups].mean() - local_values[~local_groups].mean())
        patient_differences.setdefault(patient, []).append(difference)
    if not patient_differences:
        raise ValueError("No patient/section stratum contains both comparison groups.")
    patient_means = [float(np.mean(differences)) for differences in patient_differences.values()]
    return float(np.mean(patient_means))


def stratified_permutation_test(
    data: pd.DataFrame,
    *,
    value_col: str,
    group_col: str,
    group_a: Any,
    group_b: Any,
    patient_col: str = "patient_id",
    section_col: str = "section_id",
    n_permutations: int = 9_999,
    alternative: str = "two-sided",
    seed: int = 0,
) -> dict[str, Any]:
    """Test a group difference while permuting labels within patient/section.

    The statistic is a mean difference (``group_a - group_b``) calculated first
    per section, then averaged per patient, and finally averaged across patients.
    This hierarchy prevents patients with more sections or cells from receiving
    disproportionate weight.  Label counts are preserved exactly within every
    included patient/section stratum.
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    required = {value_col, group_col, patient_col, section_col}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError("Missing permutation-test columns: " + ", ".join(missing))
    if group_a == group_b:
        raise ValueError("group_a and group_b must differ.")
    if not isinstance(n_permutations, int) or n_permutations < 1:
        raise ValueError("n_permutations must be a positive integer.")
    if alternative not in {"two-sided", "greater", "less"}:
        raise ValueError("alternative must be 'two-sided', 'greater', or 'less'.")

    working = data.loc[data[group_col].isin([group_a, group_b])].copy()
    working[value_col] = pd.to_numeric(working[value_col], errors="coerce")
    finite = np.isfinite(working[value_col].to_numpy(dtype=float))
    working = working.loc[finite].reset_index(drop=True)
    if working.empty:
        raise ValueError("No finite observations remain for the selected groups.")

    section_indices: list[np.ndarray] = []
    section_patients: list[Any] = []
    for key, group in working.groupby([patient_col, section_col], sort=True, observed=True):
        if group[group_col].nunique(dropna=False) == 2:
            section_indices.append(group.index.to_numpy(dtype=int))
            section_patients.append(_group_key(key)[0])
    if not section_indices:
        raise ValueError("No patient/section stratum contains both comparison groups.")

    included_indices = np.concatenate(section_indices)
    working = working.loc[included_indices].reset_index(drop=True)
    section_indices = []
    section_patients = []
    for key, group in working.groupby([patient_col, section_col], sort=True, observed=True):
        section_indices.append(group.index.to_numpy(dtype=int))
        section_patients.append(_group_key(key)[0])

    values = working[value_col].to_numpy(dtype=float)
    observed_labels = working[group_col].eq(group_a).to_numpy(dtype=bool)
    observed = _hierarchical_group_difference(
        values, observed_labels, section_indices, section_patients
    )

    rng = np.random.default_rng(seed)
    null_distribution = np.empty(n_permutations, dtype=float)
    for permutation in range(n_permutations):
        permuted_labels = observed_labels.copy()
        for indices in section_indices:
            permuted_labels[indices] = rng.permutation(permuted_labels[indices])
        null_distribution[permutation] = _hierarchical_group_difference(
            values, permuted_labels, section_indices, section_patients
        )

    if alternative == "greater":
        extreme = int(np.count_nonzero(null_distribution >= observed))
    elif alternative == "less":
        extreme = int(np.count_nonzero(null_distribution <= observed))
    else:
        extreme = int(np.count_nonzero(np.abs(null_distribution) >= abs(observed)))
    p_value = (extreme + 1.0) / (n_permutations + 1.0)

    return {
        "group_a": group_a,
        "group_b": group_b,
        "observed_difference": observed,
        "alternative": alternative,
        "p_value": float(p_value),
        "null_mean": float(null_distribution.mean()),
        "null_sd": float(null_distribution.std(ddof=1)) if n_permutations > 1 else np.nan,
        "n_permutations": n_permutations,
        "seed": seed,
        "n_patients": int(working[patient_col].nunique()),
        "n_sections": int(len(section_indices)),
        "n_observations": int(len(working)),
        "null_distribution": null_distribution,
    }


def _bh_adjust(p_values: np.ndarray) -> np.ndarray:
    if p_values.size == 0:
        return p_values.copy()
    order = np.argsort(p_values, kind="mergesort")
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def _directed_radius_edges(coordinates: np.ndarray, radius: float) -> tuple[np.ndarray, np.ndarray]:
    tree = cKDTree(coordinates)
    source: list[int] = []
    target: list[int] = []
    for index, neighbours in enumerate(tree.query_ball_point(coordinates, r=radius)):
        for neighbour in sorted(neighbours):
            if neighbour != index:
                source.append(index)
                target.append(int(neighbour))
    return np.asarray(source, dtype=int), np.asarray(target, dtype=int)


def neighborhood_enrichment(
    table: pd.DataFrame,
    *,
    radius: float,
    cell_type_col: str = "cell_type",
    focal_types: Sequence[str] | None = None,
    neighbor_types: Sequence[str] | None = None,
    patient_col: str = "patient_id",
    section_col: str = "section_id",
    coordinate_cols: Sequence[str] = ("x", "y"),
    coordinate_scale_um: float = 1.0,
    n_permutations: int = 999,
    seed: int = 0,
    pseudocount: float = 0.5,
    fraction_pseudocount: float = 1e-12,
) -> pd.DataFrame:
    """Estimate directed cell-type neighborhood enrichment within a radius.

    Cell-type labels are permuted independently within each patient/section,
    preserving specimen composition and the observed spatial graph.  A directed
    edge is counted from every focal cell to every non-self neighbor within the
    specified micrometre radius.

    The inferential statistic is the focal-to-neighbor edge count divided by all
    directed edges in that section, averaged first across edge-bearing sections
    within a patient and then equally across patients with at least one such
    section.  Consequently, a patient with a larger tissue section or denser
    segmentation does not receive greater cohort weight.  Sections without a
    non-self edge do not define this fraction and are excluded with eligibility
    counts reported in the output.  Raw pooled edge counts are retained as
    descriptive columns only.  Enrichment and depletion tail probabilities,
    ``z_score``, ``log2_enrichment``, and the Benjamini-Hochberg adjusted value
    are all based on the patient-equal statistic.  These probabilities test a
    conditional random-label null in the observed tissues; they are not a
    population-level test of patient effects.
    """

    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be a finite positive number.")
    if not np.isfinite(coordinate_scale_um) or coordinate_scale_um <= 0:
        raise ValueError("coordinate_scale_um must be a finite positive number.")
    if not isinstance(n_permutations, int) or n_permutations < 1:
        raise ValueError("n_permutations must be a positive integer.")
    if not np.isfinite(pseudocount) or pseudocount <= 0:
        raise ValueError("pseudocount must be a finite positive number.")
    if not np.isfinite(fraction_pseudocount) or fraction_pseudocount <= 0:
        raise ValueError("fraction_pseudocount must be a finite positive number.")

    cells = validate_spatial_table(
        table,
        patient_col=patient_col,
        section_col=section_col,
        coordinate_cols=coordinate_cols,
    )
    if cell_type_col not in cells.columns:
        raise ValueError(f"Missing cell-type column: {cell_type_col}")
    missing_labels = cells[cell_type_col].isna()
    blank_labels = cells[cell_type_col].astype(str).str.strip().eq("")
    if missing_labels.any() or blank_labels.any():
        raise ValueError(f"Column '{cell_type_col}' contains missing or blank labels.")
    cells[cell_type_col] = cells[cell_type_col].astype(str)
    cells = cells.sort_values(
        [patient_col, section_col, "cell_id"],
        kind="mergesort",
        key=lambda values: values.astype(str),
    ).reset_index(drop=True)

    observed_types = sorted(cells[cell_type_col].unique())
    focal = tuple(observed_types if focal_types is None else map(str, focal_types))
    neighbors = tuple(observed_types if neighbor_types is None else map(str, neighbor_types))
    if not focal or not neighbors:
        raise ValueError("focal_types and neighbor_types must not be empty.")
    if len(set(focal)) != len(focal) or len(set(neighbors)) != len(neighbors):
        raise ValueError("Cell-type selections must contain unique labels.")

    focal_lookup = {name: index for index, name in enumerate(focal)}
    neighbor_lookup = {name: index for index, name in enumerate(neighbors)}
    shape = (len(focal), len(neighbors))
    strata: list[tuple[Any, np.ndarray, np.ndarray, np.ndarray]] = []
    n_sections_with_edges = 0
    for key, group in cells.groupby([patient_col, section_col], sort=True, observed=True):
        coordinates = group.loc[:, coordinate_cols].to_numpy(dtype=float) * coordinate_scale_um
        source, target = _directed_radius_edges(coordinates, radius)
        if source.size:
            n_sections_with_edges += 1
        patient = _group_key(key)[0]
        strata.append((patient, group[cell_type_col].to_numpy(dtype=object), source, target))

    if n_sections_with_edges == 0:
        raise ValueError("No non-self neighbors were found within the specified radius.")
    patients_with_edges = {patient for patient, _, source, _ in strata if source.size > 0}

    def count_edges(labels_by_stratum: Sequence[np.ndarray]) -> np.ndarray:
        counts = np.zeros(shape, dtype=float)
        for labels, (_, _, source, target) in zip(labels_by_stratum, strata, strict=True):
            for source_index, target_index in zip(source, target, strict=True):
                focal_index = focal_lookup.get(str(labels[source_index]))
                neighbor_index = neighbor_lookup.get(str(labels[target_index]))
                if focal_index is not None and neighbor_index is not None:
                    counts[focal_index, neighbor_index] += 1.0
        return counts

    def patient_mean_edge_fraction(labels_by_stratum: Sequence[np.ndarray]) -> np.ndarray:
        section_fractions_by_patient: dict[Any, list[np.ndarray]] = {}
        for labels, (patient, _, source, target) in zip(labels_by_stratum, strata, strict=True):
            if source.size == 0:
                continue
            counts = np.zeros(shape, dtype=float)
            for source_index, target_index in zip(source, target, strict=True):
                focal_index = focal_lookup.get(str(labels[source_index]))
                neighbor_index = neighbor_lookup.get(str(labels[target_index]))
                if focal_index is not None and neighbor_index is not None:
                    counts[focal_index, neighbor_index] += 1.0
            section_fractions_by_patient.setdefault(patient, []).append(counts / float(source.size))
        if not section_fractions_by_patient:  # guarded by n_sections_with_edges above
            raise ValueError("No patient has a section with non-self neighbors.")
        patient_fractions = np.stack(
            [
                np.mean(section_fractions, axis=0)
                for section_fractions in section_fractions_by_patient.values()
            ],
            axis=0,
        )
        return np.mean(patient_fractions, axis=0)

    original_labels = [labels.copy() for _, labels, _, _ in strata]
    observed_counts = count_edges(original_labels)
    observed_patient_fraction = patient_mean_edge_fraction(original_labels)
    rng = np.random.default_rng(seed)
    null_counts = np.empty((n_permutations, *shape), dtype=float)
    null_patient_fractions = np.empty((n_permutations, *shape), dtype=float)
    for permutation in range(n_permutations):
        permuted = [rng.permutation(labels) for labels in original_labels]
        null_counts[permutation] = count_edges(permuted)
        null_patient_fractions[permutation] = patient_mean_edge_fraction(permuted)

    expected_counts = null_counts.mean(axis=0)
    expected_patient_fraction = null_patient_fractions.mean(axis=0)
    null_sd = (
        null_patient_fractions.std(axis=0, ddof=1) if n_permutations > 1 else np.full(shape, np.nan)
    )
    p_enrichment = (np.sum(null_patient_fractions >= observed_patient_fraction, axis=0) + 1.0) / (
        n_permutations + 1.0
    )
    p_depletion = (np.sum(null_patient_fractions <= observed_patient_fraction, axis=0) + 1.0) / (
        n_permutations + 1.0
    )
    p_two_sided = np.minimum(1.0, 2.0 * np.minimum(p_enrichment, p_depletion))
    q_values = _bh_adjust(p_two_sided.ravel()).reshape(shape)
    log2_enrichment = np.log2(
        (observed_patient_fraction + fraction_pseudocount)
        / (expected_patient_fraction + fraction_pseudocount)
    )
    pooled_log2_edge_count_enrichment = np.log2(
        (observed_counts + pseudocount) / (expected_counts + pseudocount)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        z_scores = np.divide(
            observed_patient_fraction - expected_patient_fraction,
            null_sd,
            out=np.full(shape, np.nan),
            where=null_sd > 0,
        )

    records = []
    for focal_index, focal_name in enumerate(focal):
        for neighbor_index, neighbor_name in enumerate(neighbors):
            records.append(
                {
                    "focal_type": focal_name,
                    "neighbor_type": neighbor_name,
                    "observed_edges": int(observed_counts[focal_index, neighbor_index]),
                    "expected_edges": float(expected_counts[focal_index, neighbor_index]),
                    "observed_patient_mean_edge_fraction": float(
                        observed_patient_fraction[focal_index, neighbor_index]
                    ),
                    "expected_patient_mean_edge_fraction": float(
                        expected_patient_fraction[focal_index, neighbor_index]
                    ),
                    "log2_enrichment": float(log2_enrichment[focal_index, neighbor_index]),
                    "pooled_log2_edge_count_enrichment": float(
                        pooled_log2_edge_count_enrichment[focal_index, neighbor_index]
                    ),
                    "z_score": float(z_scores[focal_index, neighbor_index]),
                    "p_enrichment": float(p_enrichment[focal_index, neighbor_index]),
                    "p_depletion": float(p_depletion[focal_index, neighbor_index]),
                    "p_two_sided": float(p_two_sided[focal_index, neighbor_index]),
                    "q_value": float(q_values[focal_index, neighbor_index]),
                    "radius_um": float(radius),
                    "null_model": "labels_permuted_within_patient_section",
                    "inference_scope": "conditional_spatial_association_in_observed_tissues",
                    "n_patients": int(len(patients_with_edges)),
                    "n_patients_with_edges": int(len(patients_with_edges)),
                    "n_patients_total": int(cells[patient_col].nunique()),
                    "n_sections_with_edges": int(n_sections_with_edges),
                    "n_sections_total": int(len(strata)),
                    "n_permutations": int(n_permutations),
                    "seed": int(seed),
                }
            )
    return (
        pd.DataFrame.from_records(records)
        .sort_values(["focal_type", "neighbor_type"], kind="mergesort")
        .reset_index(drop=True)
    )


def plot_spatial_panel(
    table: pd.DataFrame,
    *,
    patient_id: Any | None = None,
    section_id: Any | None = None,
    patient_col: str = "patient_id",
    section_col: str = "section_id",
    cell_type_col: str = "cell_type",
    structure_col: str = "structure",
    coordinate_cols: Sequence[str] = ("x", "y"),
    car_t_col: str = "is_car_t",
    ligand_source_col: str = "is_ligand_source",
    palette: Mapping[str, str] | None = None,
    point_size: float = 18.0,
    invert_y: bool = True,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Plot one real patient/section as a publication-ready spatial panel."""

    cells = validate_spatial_table(
        table,
        patient_col=patient_col,
        section_col=section_col,
        coordinate_cols=coordinate_cols,
    )
    if cell_type_col not in cells.columns:
        raise ValueError(f"Missing cell-type column: {cell_type_col}")

    selected_patient = (
        sorted(cells[patient_col].unique(), key=str)[0] if patient_id is None else patient_id
    )
    patient_rows = cells.loc[cells[patient_col].eq(selected_patient)]
    if patient_rows.empty:
        raise ValueError(f"Patient not found: {selected_patient}")
    selected_section = (
        sorted(patient_rows[section_col].unique(), key=str)[0] if section_id is None else section_id
    )
    panel = patient_rows.loc[patient_rows[section_col].eq(selected_section)].copy()
    if panel.empty:
        raise ValueError(
            f"Section '{selected_section}' was not found for patient '{selected_patient}'."
        )

    if ax is None:
        figure, axis = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    else:
        axis = ax
        figure = axis.figure

    cell_types = sorted(panel[cell_type_col].astype(str).unique())
    if palette is None:
        color_map = plt.get_cmap("tab20")
        colors = {name: color_map(index % 20) for index, name in enumerate(cell_types)}
    else:
        missing_colors = sorted(set(cell_types).difference(palette))
        if missing_colors:
            raise ValueError("Palette lacks cell types: " + ", ".join(missing_colors))
        colors = dict(palette)

    x_col, y_col = coordinate_cols[:2]
    for cell_type in cell_types:
        subset = panel.loc[panel[cell_type_col].astype(str).eq(cell_type)]
        axis.scatter(
            subset[x_col],
            subset[y_col],
            s=point_size,
            color=colors[cell_type],
            edgecolors="none",
            alpha=0.78,
            label=cell_type,
            rasterized=True,
        )

    if structure_col in panel.columns:
        structure_styles = {
            "vessel": ("#C62828", "s"),
            "tumor_nest": ("#6A3D9A", "o"),
            "stroma": ("#8D6E63", "D"),
            "lymphatic": ("#00897B", "^"),
        }
        for structure, (color, marker) in structure_styles.items():
            subset = panel.loc[panel[structure_col].astype(str).eq(structure)]
            if not subset.empty:
                axis.scatter(
                    subset[x_col],
                    subset[y_col],
                    s=point_size * 2.1,
                    facecolors="none",
                    edgecolors=color,
                    marker=marker,
                    linewidths=1.0,
                    label=f"{structure} landmark",
                )

    overlay_specs = (
        (car_t_col, "CAR-T", "#111111", 2.0),
        (ligand_source_col, "Ligand source", "#FF8F00", 1.5),
    )
    for column, label, color, width in overlay_specs:
        if column in panel.columns:
            flags = _coerce_boolean(panel[column], column)
            subset = panel.loc[flags]
            if not subset.empty:
                axis.scatter(
                    subset[x_col],
                    subset[y_col],
                    s=point_size * 3.0,
                    facecolors="none",
                    edgecolors=color,
                    linewidths=width,
                    label=label,
                )

    axis.set_title(f"Patient {selected_patient} | section {selected_section}")
    axis.set_xlabel(f"{x_col} (coordinate units)")
    axis.set_ylabel(f"{y_col} (coordinate units)")
    axis.set_aspect("equal", adjustable="datalim")
    if invert_y:
        axis.invert_yaxis()
    axis.spines[["top", "right"]].set_visible(False)
    handles, labels = axis.get_legend_handles_labels()
    unique = dict(zip(labels, handles, strict=True))
    axis.legend(
        unique.values(),
        unique.keys(),
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        fontsize=8,
    )
    return figure, axis


def plot_conceptual_spatial_panel(
    *,
    seed: int = 0,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Draw a deterministic conceptual panel for a spatial-analysis workflow."""

    if ax is None:
        figure, axis = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    else:
        axis = ax
        figure = axis.figure
    rng = np.random.default_rng(seed)

    axis.set_facecolor("#FAF8F5")
    tumor = Ellipse(
        (57, 51),
        width=50,
        height=52,
        facecolor="#E7D9F2",
        edgecolor="#6A3D9A",
        linewidth=2.0,
        alpha=0.78,
        label="Tumor nest",
    )
    axis.add_patch(tumor)

    vessel_vertices = [(5, 10), (15, 30), (8, 57), (18, 90)]
    vessel_codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    vessel = PathPatch(
        Path(vessel_vertices, vessel_codes),
        fill=False,
        edgecolor="#C62828",
        linewidth=8,
        alpha=0.82,
        capstyle="round",
        label="Blood vessel",
    )
    axis.add_patch(vessel)

    lymph_vertices = [(93, 8), (81, 28), (96, 55), (84, 90)]
    lymphatic = PathPatch(
        Path(lymph_vertices, vessel_codes),
        fill=False,
        edgecolor="#00897B",
        linewidth=5,
        alpha=0.8,
        linestyle=(0, (5, 3)),
        capstyle="round",
        label="Lymphatic",
    )
    axis.add_patch(lymphatic)

    stroma_x = rng.uniform(20, 88, 65)
    stroma_y = rng.uniform(8, 92, 65)
    outside_tumor = ((stroma_x - 57) / 25) ** 2 + ((stroma_y - 51) / 26) ** 2 > 1
    axis.scatter(
        stroma_x[outside_tumor],
        stroma_y[outside_tumor],
        s=15,
        color="#A1887F",
        alpha=0.45,
        edgecolors="none",
        label="Stroma",
    )

    ligand_x = rng.normal(54, 8, 14)
    ligand_y = rng.normal(50, 10, 14)
    axis.scatter(
        ligand_x,
        ligand_y,
        s=45,
        color="#FFB300",
        edgecolors="#8D5D00",
        linewidths=0.7,
        label="Chemokine source",
        zorder=5,
    )

    car_x = np.array([17, 25, 33, 42, 49], dtype=float)
    car_y = np.array([29, 33, 39, 45, 50], dtype=float) + rng.normal(0, 1.0, 5)
    axis.scatter(
        car_x,
        car_y,
        s=70,
        color="#2E5AAC",
        edgecolors="white",
        linewidths=0.8,
        label="CAR-T cell",
        zorder=6,
    )
    for start, end in zip(
        zip(car_x[:-1], car_y[:-1], strict=False),
        zip(car_x[1:], car_y[1:], strict=False),
        strict=True,
    ):
        axis.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=11,
                linewidth=1.1,
                color="#2E5AAC",
                alpha=0.8,
            )
        )

    axis.annotate(
        "Entry / chemotaxis",
        xy=(27, 35),
        xytext=(24, 18),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
        fontsize=9,
    )
    axis.annotate(
        "Retention near ligand sources",
        xy=(50, 50),
        xytext=(48, 78),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
        fontsize=9,
    )
    axis.annotate(
        "Potential egress",
        xy=(85, 65),
        xytext=(69, 91),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
        fontsize=9,
    )

    axis.set_xlim(0, 100)
    axis.set_ylim(0, 100)
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.spines[:].set_visible(False)
    axis.set_title("Conceptual spatial readout (not quantitative)", loc="left", weight="bold")
    handles, labels = axis.get_legend_handles_labels()
    unique = dict(zip(labels, handles, strict=True))
    axis.legend(unique.values(), unique.keys(), loc="upper left", frameon=False, fontsize=8)
    return figure, axis


__all__ = [
    "DEFAULT_LANDMARKS",
    "DEFAULT_QUERY_COLUMNS",
    "REQUIRED_SPATIAL_COLUMNS",
    "distance_to_landmarks",
    "neighborhood_enrichment",
    "plot_conceptual_spatial_panel",
    "plot_spatial_panel",
    "stratified_permutation_test",
    "summarize_patient_distances",
    "validate_spatial_table",
]
