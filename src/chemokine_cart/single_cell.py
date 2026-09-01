"""Patient-level summaries for CAR-T single-cell and CITE-seq data.

The functions in this module operate on tabular cell metadata and count
matrices represented as :class:`pandas.DataFrame` objects.  They deliberately
do not depend on Scanpy or AnnData so that the statistical summaries can be
tested independently of a particular preprocessing workflow.

Cell states and CAR status must be assigned upstream using documented gates or
classification rules.  No function in this module infers a biological state
from expression data.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Literal

import numpy as np
import pandas as pd

DEFAULT_STATE_ORDER: tuple[str, ...] = (
    "memory",
    "effector",
    "dysfunctional",
)


def _normalise_columns(columns: Iterable[str], *, argument: str) -> tuple[str, ...]:
    """Return a validated, duplicate-free tuple of column names."""

    result = tuple(columns)
    if not result:
        raise ValueError(f"{argument} must contain at least one column name")
    if any(not isinstance(column, str) or not column for column in result):
        raise TypeError(f"{argument} must contain non-empty strings")
    if len(result) != len(set(result)):
        raise ValueError(f"{argument} contains duplicate column names")
    return result


def validate_required_columns(
    table: pd.DataFrame,
    required_columns: Iterable[str],
    *,
    table_name: str = "table",
) -> None:
    """Validate that a DataFrame contains an unambiguous required schema.

    Parameters
    ----------
    table:
        Input table to validate.
    required_columns:
        Column names required by the downstream operation.
    table_name:
        Human-readable name included in validation errors.

    Raises
    ------
    TypeError
        If ``table`` is not a DataFrame or column names are invalid.
    ValueError
        If the input has duplicate column names or required columns are absent.
    """

    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"{table_name} must be a pandas DataFrame")
    required = _normalise_columns(required_columns, argument="required_columns")
    duplicated = table.columns[table.columns.duplicated()].tolist()
    if duplicated:
        duplicate_text = ", ".join(map(str, duplicated))
        raise ValueError(f"{table_name} has duplicate column names: {duplicate_text}")
    missing = [column for column in required if column not in table.columns]
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(f"{table_name} is missing required columns: {missing_text}")


def _require_complete_metadata(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    table_name: str,
) -> None:
    missing_counts = table.loc[:, list(columns)].isna().sum()
    incomplete = missing_counts[missing_counts > 0]
    if not incomplete.empty:
        details = ", ".join(f"{column} ({int(count)})" for column, count in incomplete.items())
        raise ValueError(f"{table_name} contains missing values: {details}")


def split_car_endogenous(
    cells: pd.DataFrame,
    *,
    car_positive_col: str = "car_positive",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split cells into CAR-positive and endogenous T-cell tables.

    ``car_positive_col`` must be a precomputed Boolean gate.  Integer 0/1
    values are accepted, but free-text labels and continuous expression values
    are rejected to prevent an undocumented threshold from being introduced.
    The input table is never modified.

    Returns
    -------
    car_positive, endogenous:
        Independent copies retaining all original columns and row indices.
    """

    validate_required_columns(cells, (car_positive_col,), table_name="cells")
    status = cells[car_positive_col]
    if status.isna().any():
        raise ValueError(f"cells.{car_positive_col} contains missing values")

    if pd.api.types.is_bool_dtype(status.dtype):
        mask = status.astype(bool)
    elif pd.api.types.is_integer_dtype(status.dtype):
        observed = set(status.unique().tolist())
        if not observed.issubset({0, 1}):
            raise ValueError(f"cells.{car_positive_col} must contain only Boolean or 0/1 values")
        mask = status.astype(bool)
    else:
        raise TypeError(
            f"cells.{car_positive_col} must be Boolean or integer 0/1; "
            "apply and document the CAR-detection gate upstream"
        )

    return cells.loc[mask].copy(), cells.loc[~mask].copy()


def aggregate_pseudobulk(
    cells: pd.DataFrame,
    expression_columns: Iterable[str],
    *,
    group_columns: Iterable[str] = (
        "patient_id",
        "timepoint_id",
        "t_cell_origin",
    ),
    min_cells: int = 1,
) -> pd.DataFrame:
    """Sum raw RNA UMI or ADT counts within patient-level groups.

    This operation is suitable for pseudobulk inputs to downstream count-based
    models.  Expression columns must be finite, non-negative numeric values.
    Use separate calls, or a combined set of uniquely named columns, for RNA
    and CITE-seq antibody-derived tag counts.

    Parameters
    ----------
    cells:
        One row per cell, with grouping metadata and raw count columns.
    expression_columns:
        Raw RNA or ADT count columns to sum.
    group_columns:
        Biological aggregation unit.  The default preserves patient,
        timepoint, and CAR/endogenous origin.
    min_cells:
        Groups with fewer cells are omitted.  The retained cell count is
        reported in ``n_cells``.
    """

    expression = _normalise_columns(expression_columns, argument="expression_columns")
    groups = _normalise_columns(group_columns, argument="group_columns")
    overlap = sorted(set(expression).intersection(groups))
    if overlap:
        raise ValueError("expression_columns and group_columns overlap: " + ", ".join(overlap))
    if isinstance(min_cells, bool) or not isinstance(min_cells, int) or min_cells < 1:
        raise ValueError("min_cells must be a positive integer")
    validate_required_columns(
        cells,
        (*groups, *expression),
        table_name="cells",
    )
    _require_complete_metadata(cells, groups, table_name="cells grouping metadata")

    expression_frame = cells.loc[:, list(expression)]
    non_numeric = [
        column
        for column in expression
        if not pd.api.types.is_numeric_dtype(expression_frame[column].dtype)
    ]
    if non_numeric:
        raise TypeError("expression columns must be numeric: " + ", ".join(non_numeric))
    values = expression_frame.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError("expression columns contain missing or non-finite values")
    if (values < 0).any():
        raise ValueError("pseudobulk input counts must be non-negative")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("pseudobulk input must contain untransformed integer-valued counts")

    grouped = cells.groupby(list(groups), observed=True, sort=True, dropna=False)
    totals = grouped[list(expression)].sum().reset_index()
    cell_counts = grouped.size().rename("n_cells").reset_index()
    result = cell_counts.merge(totals, on=list(groups), validate="one_to_one")
    result = result.loc[result["n_cells"] >= min_cells].reset_index(drop=True)
    return result.loc[:, [*groups, "n_cells", *expression]]


def calculate_state_fractions(
    cells: pd.DataFrame,
    *,
    group_columns: Iterable[str] = (
        "patient_id",
        "timepoint_id",
        "t_cell_origin",
    ),
    state_col: str = "cell_state",
    states: Sequence[str] | None = None,
    denominator: Literal["all", "selected"] = "all",
    other_label: str | None = None,
    min_cells_per_group: int = 1,
) -> pd.DataFrame:
    """Calculate cell-state counts and fractions within patient-level groups.

    When ``states`` is supplied, every observed patient-level group is returned
    for every requested state; an absent state receives a count and fraction of
    zero.  If ``other_label`` is set, all observed labels outside ``states`` are
    explicitly pooled into that additional category.  ``passes_min_cells`` is
    reported rather than silently removing sparsely sampled groups.

    With ``denominator='all'`` (default), fractions use all cells in the group.
    With ``'selected'``, the denominator contains only the requested states and
    groups without any selected cell are rejected.  ``other_label`` cannot be
    combined with ``denominator='selected'`` because the explicit other category
    already makes the output composition exhaustive.
    """

    groups = _normalise_columns(group_columns, argument="group_columns")
    if state_col in groups:
        raise ValueError("state_col must not also appear in group_columns")
    if denominator not in {"all", "selected"}:
        raise ValueError("denominator must be 'all' or 'selected'")
    if (
        isinstance(min_cells_per_group, bool)
        or not isinstance(min_cells_per_group, int)
        or min_cells_per_group < 1
    ):
        raise ValueError("min_cells_per_group must be a positive integer")
    validate_required_columns(cells, (*groups, state_col), table_name="cells")
    _require_complete_metadata(
        cells,
        (*groups, state_col),
        table_name="cells state metadata",
    )
    if cells.empty:
        raise ValueError("cells must contain at least one row")

    working = cells.copy()
    if states is None:
        if other_label is not None:
            raise ValueError("other_label requires an explicit states sequence")
        state_levels = tuple(pd.unique(working[state_col]).tolist())
    else:
        state_levels = _normalise_columns(states, argument="states")
        if other_label is not None:
            if not isinstance(other_label, str) or not other_label:
                raise TypeError("other_label must be a non-empty string")
            if other_label in state_levels:
                raise ValueError("other_label must not duplicate a requested state")
            if denominator == "selected":
                raise ValueError("other_label cannot be used with denominator='selected'")
            selected_mask = working[state_col].isin(state_levels)
            working.loc[~selected_mask, state_col] = other_label
            state_levels = (*state_levels, other_label)
    if not state_levels:
        raise ValueError("no cell states are available for summarisation")

    observed_counts = (
        working.groupby([*groups, state_col], observed=True, sort=False, dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    observed_groups = working.loc[:, list(groups)].drop_duplicates(ignore_index=True)
    requested_states = pd.DataFrame({state_col: list(state_levels)})
    complete = observed_groups.merge(requested_states, how="cross")
    counts = complete.merge(
        observed_counts,
        on=[*groups, state_col],
        how="left",
        validate="one_to_one",
    )
    counts["n_cells"] = counts["n_cells"].fillna(0).astype("int64")

    if denominator == "all":
        totals = (
            working.groupby(list(groups), observed=True, sort=False, dropna=False)
            .size()
            .rename("n_total")
            .reset_index()
        )
    else:
        totals = (
            counts.groupby(list(groups), observed=True, sort=False, dropna=False)["n_cells"]
            .sum()
            .rename("n_total")
            .reset_index()
        )
        zero_groups = totals.loc[totals["n_total"] == 0, list(groups)]
        if not zero_groups.empty:
            raise ValueError("at least one patient-level group contains no selected states")

    result = counts.merge(totals, on=list(groups), validate="many_to_one")
    result["fraction"] = result["n_cells"] / result["n_total"]
    result["passes_min_cells"] = result["n_total"] >= min_cells_per_group
    return result.loc[:, [*groups, state_col, "n_cells", "n_total", "fraction", "passes_min_cells"]]


def bootstrap_patient_ci(
    patient_metrics: pd.DataFrame,
    value_col: str,
    *,
    patient_col: str = "patient_id",
    group_columns: Iterable[str] = (),
    statistic: Literal["mean", "median"] | Callable[[np.ndarray], float] = "mean",
    confidence_level: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | None = 0,
    min_patients: int = 2,
    on_insufficient: Literal["raise", "record"] = "raise",
) -> pd.DataFrame:
    """Estimate percentile bootstrap intervals by resampling patients.

    The input must contain one value per patient within each stratum defined by
    ``group_columns``.  This constraint prevents cells or technical replicates
    from being treated as independent biological replicates.  Aggregate such
    replicates before calling this function.  With ``on_insufficient='record'``,
    strata below ``min_patients`` are retained with an explicit status and blank
    confidence limits rather than aborting or silently omitting them.
    """

    groups = tuple(group_columns)
    if groups:
        groups = _normalise_columns(groups, argument="group_columns")
    if patient_col in groups:
        raise ValueError("patient_col must not also appear in group_columns")
    validate_required_columns(
        patient_metrics,
        (*groups, patient_col, value_col),
        table_name="patient_metrics",
    )
    _require_complete_metadata(
        patient_metrics,
        (*groups, patient_col, value_col),
        table_name="patient_metrics",
    )
    if patient_metrics.empty:
        raise ValueError("patient_metrics must contain at least one row")
    if not pd.api.types.is_numeric_dtype(patient_metrics[value_col].dtype):
        raise TypeError(f"patient_metrics.{value_col} must be numeric")
    values = patient_metrics[value_col].to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError(f"patient_metrics.{value_col} contains non-finite values")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between 0 and 1")
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 100:
        raise ValueError("n_boot must be an integer of at least 100")
    if isinstance(min_patients, bool) or not isinstance(min_patients, int) or min_patients < 2:
        raise ValueError("min_patients must be an integer of at least 2")
    if on_insufficient not in {"raise", "record"}:
        raise ValueError("on_insufficient must be 'raise' or 'record'")

    duplicate_key = [*groups, patient_col]
    duplicated = patient_metrics.duplicated(duplicate_key, keep=False)
    if duplicated.any():
        raise ValueError(
            "patient_metrics must contain one value per patient within each group; "
            "aggregate repeated samples before bootstrapping"
        )

    if statistic == "mean":
        statistic_function: Callable[[np.ndarray], float] = np.mean
        vectorised_statistic = "mean"
    elif statistic == "median":
        statistic_function = np.median
        vectorised_statistic = "median"
    elif callable(statistic):
        statistic_function = statistic
        vectorised_statistic = None
    else:
        raise ValueError("statistic must be 'mean', 'median', or a callable")

    rng = np.random.default_rng(random_state)
    alpha = (1.0 - confidence_level) / 2.0
    records: list[dict[str, object]] = []

    if groups:
        iterator = patient_metrics.groupby(list(groups), observed=True, sort=True, dropna=False)
    else:
        iterator = [((), patient_metrics)]

    for group_key, frame in iterator:
        if groups and not isinstance(group_key, tuple):
            group_key = (group_key,)
        ordered = frame.assign(
            __patient_sort=frame[patient_col].map(lambda value: f"{type(value).__name__}:{value}")
        ).sort_values("__patient_sort", kind="stable")
        patient_values = ordered[value_col].to_numpy(dtype=float)
        n_patients = patient_values.size
        if n_patients < min_patients:
            label = dict(zip(groups, group_key, strict=False)) if groups else {"all": "all"}
            if on_insufficient == "raise":
                raise ValueError(
                    f"group {label} has {n_patients} patients; at least {min_patients} are required"
                )
            record: dict[str, object] = {}
            if groups:
                record.update(dict(zip(groups, group_key, strict=False)))
            record.update(
                {
                    "n_patients": int(n_patients),
                    "estimate": float(statistic_function(patient_values)),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "confidence_level": float(confidence_level),
                    "n_boot": int(n_boot),
                    "ci_method": "patient_percentile_bootstrap",
                    "status": "insufficient_patients",
                }
            )
            records.append(record)
            continue

        sample_indices = rng.integers(
            0,
            n_patients,
            size=(n_boot, n_patients),
            endpoint=False,
        )
        sampled = patient_values[sample_indices]
        if vectorised_statistic == "mean":
            bootstrap_statistics = sampled.mean(axis=1)
        elif vectorised_statistic == "median":
            bootstrap_statistics = np.median(sampled, axis=1)
        else:
            bootstrap_statistics = np.asarray(
                [statistic_function(sample) for sample in sampled],
                dtype=float,
            )
        if not np.isfinite(bootstrap_statistics).all():
            raise ValueError("statistic returned non-finite bootstrap estimates")

        record: dict[str, object] = {}
        if groups:
            record.update(dict(zip(groups, group_key, strict=False)))
        record.update(
            {
                "n_patients": int(n_patients),
                "estimate": float(statistic_function(patient_values)),
                "ci_low": float(np.quantile(bootstrap_statistics, alpha)),
                "ci_high": float(np.quantile(bootstrap_statistics, 1.0 - alpha)),
                "confidence_level": float(confidence_level),
                "n_boot": int(n_boot),
                "ci_method": "patient_percentile_bootstrap",
                "status": "ok",
            }
        )
        records.append(record)

    return pd.DataFrame.from_records(records)


def bootstrap_paired_change_ci(
    patient_metrics: pd.DataFrame,
    value_col: str,
    *,
    reference_timepoint: str,
    timepoint_col: str = "timepoint_id",
    patient_col: str = "patient_id",
    group_columns: Iterable[str] = (),
    comparison_timepoints: Sequence[str] | None = None,
    statistic: Literal["mean", "median"] | Callable[[np.ndarray], float] = "mean",
    confidence_level: float = 0.95,
    n_boot: int = 2_000,
    random_state: int | None = 0,
    min_patients: int = 2,
    on_insufficient: Literal["raise", "record"] = "record",
) -> pd.DataFrame:
    """Estimate paired longitudinal changes by resampling complete patients.

    Each comparison is calculated as ``comparison - reference`` within a
    patient.  Only complete patient pairs contribute.  Technical or biological
    sample replicates must be resolved upstream; duplicate patient/timepoint
    rows within an analysis stratum are rejected.
    """

    groups = tuple(group_columns)
    if groups:
        groups = _normalise_columns(groups, argument="group_columns")
    forbidden = {patient_col, timepoint_col}.intersection(groups)
    if forbidden:
        raise ValueError(
            "group_columns must not contain patient or timepoint columns: "
            + ", ".join(sorted(forbidden))
        )
    if not isinstance(reference_timepoint, str) or not reference_timepoint:
        raise TypeError("reference_timepoint must be a non-empty string")
    validate_required_columns(
        patient_metrics,
        (*groups, patient_col, timepoint_col, value_col),
        table_name="patient_metrics",
    )
    _require_complete_metadata(
        patient_metrics,
        (*groups, patient_col, timepoint_col, value_col),
        table_name="patient_metrics",
    )
    if patient_metrics.empty:
        raise ValueError("patient_metrics must contain at least one row")
    if not pd.api.types.is_numeric_dtype(patient_metrics[value_col].dtype):
        raise TypeError(f"patient_metrics.{value_col} must be numeric")
    values = patient_metrics[value_col].to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError(f"patient_metrics.{value_col} contains non-finite values")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between 0 and 1")
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 100:
        raise ValueError("n_boot must be an integer of at least 100")
    if isinstance(min_patients, bool) or not isinstance(min_patients, int) or min_patients < 2:
        raise ValueError("min_patients must be an integer of at least 2")
    if on_insufficient not in {"raise", "record"}:
        raise ValueError("on_insufficient must be 'raise' or 'record'")

    duplicate_key = [*groups, patient_col, timepoint_col]
    if patient_metrics.duplicated(duplicate_key, keep=False).any():
        raise ValueError(
            "patient_metrics must contain one sample per patient/timepoint within "
            "each group before paired bootstrapping"
        )

    observed_timepoints = tuple(pd.unique(patient_metrics[timepoint_col]).tolist())
    if reference_timepoint not in observed_timepoints:
        raise ValueError(f"reference timepoint {reference_timepoint!r} is absent")
    if comparison_timepoints is None:
        comparisons = tuple(
            timepoint for timepoint in observed_timepoints if timepoint != reference_timepoint
        )
    else:
        comparisons = _normalise_columns(
            comparison_timepoints,
            argument="comparison_timepoints",
        )
        if reference_timepoint in comparisons:
            raise ValueError("comparison_timepoints must not include reference_timepoint")
        absent = set(comparisons).difference(observed_timepoints)
        if absent:
            raise ValueError("comparison_timepoints are absent: " + ", ".join(sorted(absent)))
    if not comparisons:
        raise ValueError("at least one comparison timepoint is required")

    if groups:
        iterator = patient_metrics.groupby(list(groups), observed=True, sort=True, dropna=False)
    else:
        iterator = [((), patient_metrics)]

    records: list[dict[str, object]] = []
    for group_key, frame in iterator:
        if groups and not isinstance(group_key, tuple):
            group_key = (group_key,)
        reference = frame.loc[
            frame[timepoint_col] == reference_timepoint,
            [patient_col, value_col],
        ].rename(columns={value_col: "reference_value"})
        for comparison in comparisons:
            follow_up = frame.loc[
                frame[timepoint_col] == comparison,
                [patient_col, value_col],
            ].rename(columns={value_col: "comparison_value"})
            paired = reference.merge(follow_up, on=patient_col, validate="one_to_one")
            paired["paired_change"] = paired["comparison_value"] - paired["reference_value"]

            base: dict[str, object] = {}
            if groups:
                base.update(dict(zip(groups, group_key, strict=False)))
            base.update(
                {
                    "reference_timepoint": reference_timepoint,
                    "comparison_timepoint": comparison,
                }
            )
            if paired.empty:
                if on_insufficient == "raise":
                    raise ValueError(f"group {base} contains no complete patient pairs")
                base.update(
                    {
                        "n_patients": 0,
                        "estimate": np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "confidence_level": float(confidence_level),
                        "n_boot": int(n_boot),
                        "ci_method": "paired_patient_percentile_bootstrap",
                        "status": "insufficient_patients",
                    }
                )
                records.append(base)
                continue

            interval = bootstrap_patient_ci(
                paired.loc[:, [patient_col, "paired_change"]],
                "paired_change",
                patient_col=patient_col,
                statistic=statistic,
                confidence_level=confidence_level,
                n_boot=n_boot,
                random_state=random_state,
                min_patients=min_patients,
                on_insufficient=on_insufficient,
            ).iloc[0]
            base.update(interval.to_dict())
            base["ci_method"] = "paired_patient_percentile_bootstrap"
            records.append(base)

    return pd.DataFrame.from_records(records)


def plot_state_fractions(
    state_fractions: pd.DataFrame,
    *,
    states: Sequence[str] = DEFAULT_STATE_ORDER,
    patient_col: str = "patient_id",
    timepoint_col: str = "timepoint_id",
    origin_col: str = "t_cell_origin",
    state_col: str = "cell_state",
    value_col: str = "fraction",
    timepoint_order: Sequence[object] | None = None,
    palette: Mapping[object, str] | None = None,
    connect_patient_samples: bool = True,
    axes: Sequence[object] | None = None,
) -> tuple[object, np.ndarray]:
    """Plot observed patient trajectories for CAR-T cell-state fractions.

    One panel is drawn for each requested state.  Every marker corresponds to
    an input observation; no interpolation, imputation, or synthetic example
    values are added.  Lines, when enabled, connect serial observations from
    the same patient and T-cell origin.

    Returns the Matplotlib figure and a one-dimensional array of axes.
    """

    import matplotlib.pyplot as plt

    requested_states = _normalise_columns(states, argument="states")
    required = (
        patient_col,
        timepoint_col,
        origin_col,
        state_col,
        value_col,
    )
    validate_required_columns(
        state_fractions,
        required,
        table_name="state_fractions",
    )
    _require_complete_metadata(
        state_fractions,
        required,
        table_name="state_fractions",
    )
    if not pd.api.types.is_numeric_dtype(state_fractions[value_col].dtype):
        raise TypeError(f"state_fractions.{value_col} must be numeric")
    values = state_fractions[value_col].to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError(f"state_fractions.{value_col} contains non-finite values")
    if ((values < 0) | (values > 1)).any():
        raise ValueError(f"state_fractions.{value_col} must be between 0 and 1")

    available_states = set(state_fractions[state_col].unique().tolist())
    absent = [state for state in requested_states if state not in available_states]
    if absent:
        raise ValueError("requested states are absent from state_fractions: " + ", ".join(absent))
    selected = state_fractions.loc[state_fractions[state_col].isin(requested_states)].copy()
    observation_key = [patient_col, timepoint_col, origin_col, state_col]
    if selected.duplicated(observation_key, keep=False).any():
        raise ValueError(
            "state_fractions contains repeated patient/timepoint/origin/state rows; "
            "aggregate technical replicates before plotting"
        )

    if timepoint_order is None:
        timepoints = tuple(pd.unique(selected[timepoint_col]).tolist())
    else:
        timepoints = tuple(timepoint_order)
        if len(timepoints) != len(set(timepoints)):
            raise ValueError("timepoint_order contains duplicates")
        missing_timepoints = set(selected[timepoint_col]).difference(timepoints)
        if missing_timepoints:
            missing_text = ", ".join(map(str, missing_timepoints))
            raise ValueError("timepoint_order omits observed timepoints: " + missing_text)
    x_positions = {timepoint: position for position, timepoint in enumerate(timepoints)}

    if axes is None:
        figure, axes_array = plt.subplots(
            1,
            len(requested_states),
            figsize=(4.2 * len(requested_states), 3.8),
            sharey=True,
            squeeze=False,
        )
        plot_axes = axes_array.ravel()
    else:
        plot_axes = np.asarray(tuple(axes), dtype=object).ravel()
        if plot_axes.size != len(requested_states):
            raise ValueError("axes must contain one axis per requested state")
        figure = plot_axes[0].figure

    origins = tuple(pd.unique(selected[origin_col]).tolist())
    if palette is None:
        default_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        colours = {
            origin: default_colours[index % len(default_colours)]
            for index, origin in enumerate(origins)
        }
    else:
        missing_origins = [origin for origin in origins if origin not in palette]
        if missing_origins:
            raise ValueError("palette is missing origins: " + ", ".join(map(str, missing_origins)))
        colours = dict(palette)

    offsets = {
        origin: (index - (len(origins) - 1) / 2.0) * 0.10 for index, origin in enumerate(origins)
    }
    legend_handles: dict[object, object] = {}
    for axis, state in zip(plot_axes, requested_states, strict=False):
        state_table = selected.loc[selected[state_col] == state]
        for origin in origins:
            origin_table = state_table.loc[state_table[origin_col] == origin]
            for _, patient_table in origin_table.groupby(patient_col, observed=True, sort=False):
                patient_table = patient_table.assign(
                    __x=patient_table[timepoint_col].map(x_positions).astype(float)
                    + offsets[origin]
                ).sort_values("__x", kind="stable")
                x = patient_table["__x"].to_numpy(dtype=float)
                y = patient_table[value_col].to_numpy(dtype=float)
                if connect_patient_samples and len(patient_table) > 1:
                    axis.plot(x, y, color=colours[origin], alpha=0.28, linewidth=0.8)
                handle = axis.scatter(
                    x,
                    y,
                    color=colours[origin],
                    edgecolor="white",
                    linewidth=0.45,
                    s=28,
                    alpha=0.9,
                    label=str(origin),
                    zorder=3,
                )
                legend_handles.setdefault(origin, handle)
        axis.set_title(str(state))
        axis.set_xticks(range(len(timepoints)), [str(value) for value in timepoints])
        axis.set_xlabel(str(timepoint_col).replace("_", " ").title())
        axis.set_ylim(-0.02, 1.02)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)

    plot_axes[0].set_ylabel("Fraction of cells")
    if legend_handles:
        figure.legend(
            list(legend_handles.values()),
            [str(origin) for origin in legend_handles],
            title=str(origin_col).replace("_", " ").title(),
            loc="upper center",
            bbox_to_anchor=(0.5, 1.04),
            ncol=max(1, len(legend_handles)),
            frameon=False,
        )
    figure.tight_layout()
    return figure, plot_axes


__all__ = [
    "DEFAULT_STATE_ORDER",
    "aggregate_pseudobulk",
    "bootstrap_patient_ci",
    "bootstrap_paired_change_ci",
    "calculate_state_fractions",
    "plot_state_fractions",
    "split_car_endogenous",
    "validate_required_columns",
]
