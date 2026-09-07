"""Patient-level aggregation and locked longitudinal contrasts."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import benjamini_hochberg, exact_two_sided_sign_test, holm_adjust, quantile
from .constants import (
    CELL_COLUMNS,
    EARLY_FAMILY_ID,
    MARKER_COUNT_COLUMNS,
    MARKER_SETS,
    SAMPLE_COLUMNS,
)


def validate_cells(cells: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize the standard cell-level interchange table."""
    missing = set(CELL_COLUMNS).difference(cells.columns)
    if missing:
        raise ValueError(f"Cell table is missing columns: {sorted(missing)}")
    marker_columns = [column for column in MARKER_COUNT_COLUMNS if column in cells]
    normalized = cells.loc[:, [*CELL_COLUMNS, *marker_columns]].copy()
    identity = ["dataset", "biological_sample_id", "cell_id"]
    if normalized.duplicated(identity).any():
        raise ValueError("Cell identifiers must be unique within each biological sample")
    for column in ("total_umi", "cxcr6_umi", *marker_columns):
        values = pd.to_numeric(normalized[column], errors="raise")
        if (~np.isfinite(values)).any() or (values < 0).any():
            raise ValueError(f"{column} must contain finite non-negative counts")
        if not np.equal(values, np.floor(values)).all():
            raise ValueError(f"{column} must contain raw integer UMI counts")
        normalized[column] = values.astype(np.int64)
    if (normalized["cxcr6_umi"] > normalized["total_umi"]).any():
        raise ValueError("CXCR6 UMI cannot exceed the cell total UMI")
    return normalized


def summarize_markers(cells: pd.DataFrame) -> pd.DataFrame:
    """Aggregate locked marker-set detection without assigning cell states."""
    cells = validate_cells(cells)
    biological_key = ["dataset", "patient_id", "product", "stage"]
    rows: list[dict[str, object]] = []
    for key, group in cells.groupby(biological_key, sort=True, observed=True, dropna=False):
        record = dict(zip(biological_key, key, strict=True))
        for module, genes in MARKER_SETS.items():
            columns = [f"gene_{gene}_umi" for gene in genes]
            missing = [column for column in columns if column not in group]
            if missing:
                continue
            detected = group.loc[:, columns].to_numpy() > 0
            rows.append(
                {
                    **record,
                    "module": module,
                    "genes": ";".join(genes),
                    "n_cells": len(group),
                    "n_genes": len(genes),
                    "detected_gene_cell_pairs": int(detected.sum()),
                    "possible_gene_cell_pairs": int(detected.size),
                    "mean_marker_detection_fraction": float(detected.mean()),
                }
            )
    return pd.DataFrame(rows)


def summarize_samples(cells: pd.DataFrame) -> pd.DataFrame:
    """Collapse cells and technical records to one row per patient and stage."""
    cells = validate_cells(cells)
    biological_key = ["dataset", "patient_id", "product", "stage"]
    sample_identity = ["dataset", "biological_sample_id"]
    sample_mapping_counts = cells.groupby(sample_identity, dropna=False)[
        ["patient_id", "product", "stage"]
    ].nunique(dropna=False)
    if sample_mapping_counts.ne(1).any(axis=1).any():
        raise ValueError("A biological sample must map to exactly one patient, product, and stage")
    grouped = (
        cells.assign(cxcr6_detected=cells["cxcr6_umi"].gt(0))
        .groupby(biological_key, sort=True, observed=True, dropna=False)
        .agg(
            biological_sample_id=(
                "biological_sample_id",
                lambda values: ";".join(sorted(set(map(str, values)))),
            ),
            n_technical_records=("biological_sample_id", "nunique"),
            eligible_cells=("cell_id", "size"),
            cxcr6_positive_cells=("cxcr6_detected", "sum"),
            cxcr6_umi=("cxcr6_umi", "sum"),
            total_umi=("total_umi", "sum"),
        )
        .reset_index()
    )
    grouped["cxcr6_fraction"] = grouped["cxcr6_positive_cells"] / grouped["eligible_cells"]
    grouped["cxcr6_cpm"] = np.divide(
        grouped["cxcr6_umi"] * 1_000_000.0,
        grouped["total_umi"],
        out=np.full(len(grouped), np.nan),
        where=grouped["total_umi"].to_numpy() > 0,
    )
    result = grouped.loc[:, SAMPLE_COLUMNS]
    if result.duplicated(biological_key).any():
        raise AssertionError("Aggregation produced duplicate patient-stage rows")
    return result.sort_values(biological_key, kind="mergesort").reset_index(drop=True)


def validate_samples(samples: pd.DataFrame) -> pd.DataFrame:
    """Validate an already aggregated patient-stage table."""
    missing = set(SAMPLE_COLUMNS).difference(samples.columns)
    if missing:
        raise ValueError(f"Sample table is missing columns: {sorted(missing)}")
    samples = samples.loc[:, SAMPLE_COLUMNS].copy()
    key = ["dataset", "patient_id", "product", "stage"]
    if samples.duplicated(key).any():
        raise ValueError("Sample table has duplicate patient-stage rows")
    integer_columns = (
        "n_technical_records",
        "eligible_cells",
        "cxcr6_positive_cells",
        "cxcr6_umi",
        "total_umi",
    )
    for column in integer_columns:
        values = pd.to_numeric(samples[column], errors="raise")
        if (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"{column} must contain non-negative integers")
        samples[column] = values.astype(np.int64)
    expected_fraction = samples["cxcr6_positive_cells"] / samples["eligible_cells"]
    if not np.allclose(samples["cxcr6_fraction"], expected_fraction, rtol=0, atol=1e-14):
        raise ValueError("cxcr6_fraction does not match the count-derived fraction")
    return samples


def _validate_gse162975_stage_aliases(aliases: pd.DataFrame) -> pd.DataFrame:
    required = {
        "deposited_sampling_stage",
        "canonical_stage",
        "selection_priority",
    }
    missing = required.difference(aliases.columns)
    if missing:
        raise ValueError(f"GSE162975 stage-alias table is missing columns: {sorted(missing)}")
    rules = aliases.loc[aliases["selection_priority"].astype(str).ne("")].copy()
    if rules.empty:
        raise ValueError("GSE162975 stage-alias table contains no prioritized alias rules")
    priority = pd.to_numeric(rules["selection_priority"], errors="raise")
    if not np.equal(priority, np.floor(priority)).all() or priority.le(0).any():
        raise ValueError("GSE162975 alias priorities must be positive integers")
    rules["selection_priority"] = priority.astype(int)
    if rules.duplicated(["canonical_stage", "selection_priority"]).any():
        raise ValueError("GSE162975 alias priorities must be unique within each canonical stage")
    if rules.duplicated(["canonical_stage", "deposited_sampling_stage"]).any():
        raise ValueError("A deposited stage may appear only once within each GSE162975 alias")
    if rules["canonical_stage"].eq(rules["deposited_sampling_stage"]).any():
        raise ValueError("Prioritized aliases must differ from their deposited stages")
    return rules.sort_values(
        ["canonical_stage", "selection_priority"], kind="mergesort"
    ).reset_index(drop=True)


def add_gse162975_stage_aliases(frame: pd.DataFrame, aliases: pd.DataFrame) -> pd.DataFrame:
    """Add one configured post-peak and extended stage per GSE162975 patient.

    The same operation is applied to patient-stage CXCR6 aggregates and marker
    aggregates so that both analyses use identical preselected source stages.
    """
    if frame.empty:
        return frame.copy()
    required = {"dataset", "patient_id", "product", "stage"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Stage table is missing columns: {sorted(missing)}")
    rules = _validate_gse162975_stage_aliases(aliases)
    result = frame.copy()
    cohort = result.loc[result["dataset"].eq("GSE162975")]
    additions: list[pd.DataFrame] = []
    for _, patient_rows in cohort.groupby(
        ["patient_id", "product"], sort=True, observed=True, dropna=False
    ):
        available_stages = set(patient_rows["stage"].astype(str))
        for alias, alias_rules in rules.groupby("canonical_stage", sort=True, observed=True):
            source_stage = next(
                (
                    stage
                    for stage in alias_rules["deposited_sampling_stage"].astype(str)
                    if stage in available_stages
                ),
                None,
            )
            if source_stage is None:
                continue
            selected = patient_rows.loc[patient_rows["stage"].eq(source_stage)].copy()
            selected["stage"] = alias
            additions.append(selected)
    if additions:
        result = pd.concat([result, *additions], ignore_index=True)
    identity = ["dataset", "patient_id", "product", "stage"]
    if "module" in result:
        identity.append("module")
    if result.duplicated(identity).any():
        raise ValueError("Stage aliasing produced duplicate patient-stage rows")
    if set(SAMPLE_COLUMNS).issubset(result.columns):
        return validate_samples(result)
    return result.sort_values(identity, kind="mergesort").reset_index(drop=True)


def build_contrast_inclusion(samples: pd.DataFrame, contrasts: pd.DataFrame) -> pd.DataFrame:
    """Create a complete, auditable patient inclusion table for every contrast."""
    samples = validate_samples(samples)
    required = {
        "dataset",
        "contrast",
        "baseline_stage",
        "followup_stage",
        "min_cells",
        "role",
        "family_id",
    }
    missing = required.difference(contrasts.columns)
    if missing:
        raise ValueError(f"Contrast table is missing columns: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for contrast in contrasts.to_dict("records"):
        cohort = samples.loc[samples["dataset"].eq(contrast["dataset"])]
        patients = cohort[["patient_id", "product"]].drop_duplicates()
        for patient in patients.to_dict("records"):
            selection = cohort[
                cohort["patient_id"].eq(patient["patient_id"])
                & cohort["product"].eq(patient["product"])
            ]
            baseline = selection.loc[selection["stage"].eq(contrast["baseline_stage"])]
            followup = selection.loc[selection["stage"].eq(contrast["followup_stage"])]
            reason = "included"
            if len(baseline) != 1 or len(followup) != 1:
                reason = "missing_or_nonunique_stage"
            elif int(baseline.iloc[0]["eligible_cells"]) < int(contrast["min_cells"]) or int(
                followup.iloc[0]["eligible_cells"]
            ) < int(contrast["min_cells"]):
                reason = "below_minimum_eligible_cells"
            included = reason == "included"
            row: dict[str, object] = {
                **contrast,
                **patient,
                "baseline_eligible_cells": (
                    int(baseline.iloc[0]["eligible_cells"]) if len(baseline) == 1 else pd.NA
                ),
                "followup_eligible_cells": (
                    int(followup.iloc[0]["eligible_cells"]) if len(followup) == 1 else pd.NA
                ),
                "baseline_fraction": (
                    float(baseline.iloc[0]["cxcr6_fraction"]) if len(baseline) == 1 else np.nan
                ),
                "followup_fraction": (
                    float(followup.iloc[0]["cxcr6_fraction"]) if len(followup) == 1 else np.nan
                ),
                "included": included,
                "reason": reason,
            }
            row["paired_fraction_change"] = (
                row["followup_fraction"] - row["baseline_fraction"] if included else np.nan
            )
            rows.append(row)
    result = pd.DataFrame(rows)
    order = [
        "dataset",
        "contrast",
        "role",
        "family_id",
        "baseline_stage",
        "followup_stage",
        "min_cells",
        "patient_id",
        "product",
        "baseline_eligible_cells",
        "followup_eligible_cells",
        "baseline_fraction",
        "followup_fraction",
        "paired_fraction_change",
        "included",
        "reason",
    ]
    if result.empty:
        return pd.DataFrame(columns=order)
    return (
        result.loc[:, order]
        .sort_values(["dataset", "contrast", "patient_id"], kind="mergesort")
        .reset_index(drop=True)
    )


def summarize_contrast(inclusion: pd.DataFrame) -> dict[str, object]:
    """Summarize one locked contrast at the independent-patient level."""
    identity = inclusion[["dataset", "contrast", "role", "family_id", "min_cells"]]
    if len(identity.drop_duplicates()) != 1:
        raise ValueError("A contrast summary must contain exactly one contrast definition")
    used = inclusion.loc[inclusion["included"]].copy()
    if used["patient_id"].duplicated().any():
        raise ValueError("A patient appears more than once in a contrast")
    changes = used["paired_fraction_change"].astype(float)
    increases = int(changes.gt(0).sum())
    decreases = int(changes.lt(0).sum())
    unchanged = int(changes.eq(0).sum())
    definition = identity.iloc[0].to_dict()
    return {
        **definition,
        "n_pairs": len(changes),
        "n_nonzero": increases + decreases,
        "n_increase": increases,
        "n_decrease": decreases,
        "n_unchanged": unchanged,
        "median_change_fraction": float(changes.median()) if len(changes) else np.nan,
        "median_change_percentage_points": (
            float(changes.median() * 100.0) if len(changes) else np.nan
        ),
        "q1_change_fraction": quantile(changes, 0.25) if len(changes) else np.nan,
        "q3_change_fraction": quantile(changes, 0.75) if len(changes) else np.nan,
        "sign_test_p": exact_two_sided_sign_test(increases, decreases),
    }


def summarize_gse197268_product_strata(inclusion: pd.DataFrame) -> pd.DataFrame:
    """Summarize the locked GSE197268 contrast by product without testing an interaction."""
    cohort = inclusion.loc[
        inclusion["dataset"].eq("GSE197268") & inclusion["contrast"].eq("D7-CART_minus_IP")
    ]
    rows: list[dict[str, object]] = []
    for product, group in cohort.groupby("product", sort=True, observed=True):
        summary = summarize_contrast(group)
        rows.append(
            {
                "dataset": summary["dataset"],
                "contrast": summary["contrast"],
                "subgroup": product,
                "cell_definition": "author-QC CAR-positive CD8 T cells",
                "role": "descriptive product stratum",
                "family_id": "none",
                "min_cells": summary["min_cells"],
                "n_pairs": summary["n_pairs"],
                "n_nonzero": summary["n_nonzero"],
                "n_increase": summary["n_increase"],
                "n_decrease": summary["n_decrease"],
                "n_unchanged": summary["n_unchanged"],
                "median_change_fraction": summary["median_change_fraction"],
                "median_change_percentage_points": summary["median_change_percentage_points"],
                "q1_change_fraction": summary["q1_change_fraction"],
                "q3_change_fraction": summary["q3_change_fraction"],
                "sign_test_p": summary["sign_test_p"],
                "holm_p": np.nan,
            }
        )
    return pd.DataFrame(rows)


def analyze_contrasts(
    samples: pd.DataFrame, contrasts: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return patient inclusion rows and separately summarized contrasts."""
    inclusion = build_contrast_inclusion(samples, contrasts)
    if inclusion.empty:
        return inclusion, pd.DataFrame()
    results = pd.DataFrame(
        [
            summarize_contrast(group)
            for _, group in inclusion.groupby(["dataset", "contrast"], sort=True, observed=True)
        ]
    )
    results["holm_p"] = np.nan
    early = results["family_id"].eq(EARLY_FAMILY_ID)
    if early.any():
        results.loc[early, "holm_p"] = holm_adjust(results.loc[early, "sign_test_p"])
    return inclusion, results.sort_values(["dataset", "contrast"], kind="mergesort").reset_index(
        drop=True
    )


def leave_one_out(inclusion: pd.DataFrame) -> pd.DataFrame:
    """Report leave-one-patient-out medians for eligible contrasts with n >= 5."""
    rows: list[dict[str, object]] = []
    for (dataset, contrast), group in inclusion.groupby(
        ["dataset", "contrast"], sort=True, observed=True
    ):
        used = group.loc[group["included"]]
        if len(used) < 5:
            continue
        for patient in used["patient_id"]:
            retained = used.loc[used["patient_id"].ne(patient), "paired_fraction_change"].astype(
                float
            )
            rows.append(
                {
                    "dataset": dataset,
                    "contrast": contrast,
                    "omitted_patient_id": patient,
                    "n_retained": len(retained),
                    "median_change_fraction": float(retained.median()),
                }
            )
    return pd.DataFrame(rows)


def analyze_marker_contrasts(markers: pd.DataFrame, inclusion: pd.DataFrame) -> pd.DataFrame:
    """Summarize marker changes in the locked CXCR6 patient-inclusion frame."""
    empty = pd.DataFrame(
        columns=[
            "dataset",
            "contrast",
            "module",
            "n_pairs",
            "n_increase",
            "n_decrease",
            "n_unchanged",
            "median_change_fraction",
            "sign_test_p",
            "bh_q_within_contrast",
        ]
    )
    if markers.empty or inclusion.empty:
        return empty
    marker_required = {
        "dataset",
        "patient_id",
        "product",
        "stage",
        "module",
        "n_cells",
        "mean_marker_detection_fraction",
    }
    missing = marker_required.difference(markers.columns)
    if missing:
        raise ValueError(f"Marker table is missing columns: {sorted(missing)}")
    inclusion_required = {
        "dataset",
        "contrast",
        "baseline_stage",
        "followup_stage",
        "patient_id",
        "product",
        "baseline_eligible_cells",
        "followup_eligible_cells",
        "included",
    }
    missing = inclusion_required.difference(inclusion.columns)
    if missing:
        raise ValueError(f"Contrast-inclusion table is missing columns: {sorted(missing)}")
    if markers.duplicated(["dataset", "patient_id", "product", "stage", "module"]).any():
        raise ValueError("Marker table has duplicate patient-stage-module rows")
    rows: list[dict[str, object]] = []
    for (dataset, contrast_name), contrast_rows in inclusion.groupby(
        ["dataset", "contrast"], sort=True, observed=True
    ):
        definition = contrast_rows[["baseline_stage", "followup_stage"]].drop_duplicates()
        if len(definition) != 1:
            raise ValueError("Marker contrast must contain exactly one stage definition")
        used = contrast_rows.loc[contrast_rows["included"]].copy()
        if used.empty:
            continue
        if used.duplicated(["patient_id", "product"]).any():
            raise ValueError("A patient appears more than once in a marker contrast")
        baseline_stage = str(definition.iloc[0]["baseline_stage"])
        followup_stage = str(definition.iloc[0]["followup_stage"])
        cohort = markers.loc[markers["dataset"].eq(dataset)]
        for module, module_rows in cohort.groupby("module", sort=True, observed=True):
            baseline = module_rows.loc[
                module_rows["stage"].eq(baseline_stage),
                [
                    "patient_id",
                    "product",
                    "n_cells",
                    "mean_marker_detection_fraction",
                ],
            ].rename(
                columns={
                    "n_cells": "marker_baseline_cells",
                    "mean_marker_detection_fraction": "baseline_marker_fraction",
                }
            )
            followup = module_rows.loc[
                module_rows["stage"].eq(followup_stage),
                [
                    "patient_id",
                    "product",
                    "n_cells",
                    "mean_marker_detection_fraction",
                ],
            ].rename(
                columns={
                    "n_cells": "marker_followup_cells",
                    "mean_marker_detection_fraction": "followup_marker_fraction",
                }
            )
            paired = used.merge(
                baseline,
                on=["patient_id", "product"],
                how="left",
                validate="one_to_one",
            ).merge(
                followup,
                on=["patient_id", "product"],
                how="left",
                validate="one_to_one",
            )
            marker_columns = [
                "marker_baseline_cells",
                "baseline_marker_fraction",
                "marker_followup_cells",
                "followup_marker_fraction",
            ]
            if paired[marker_columns].isna().any().any():
                raise ValueError(
                    f"Missing marker aggregate for an included pair: {dataset} {contrast_name} {module}"
                )
            if not np.array_equal(
                paired["marker_baseline_cells"].astype(int),
                paired["baseline_eligible_cells"].astype(int),
            ) or not np.array_equal(
                paired["marker_followup_cells"].astype(int),
                paired["followup_eligible_cells"].astype(int),
            ):
                raise ValueError(
                    f"Marker and CXCR6 eligible-cell frames differ: {dataset} {contrast_name} {module}"
                )
            change = paired["followup_marker_fraction"] - paired["baseline_marker_fraction"]
            increases = int(change.gt(0).sum())
            decreases = int(change.lt(0).sum())
            rows.append(
                {
                    "dataset": dataset,
                    "contrast": contrast_name,
                    "module": module,
                    "n_pairs": len(change),
                    "n_increase": increases,
                    "n_decrease": decreases,
                    "n_unchanged": int(change.eq(0).sum()),
                    "median_change_fraction": float(change.median()),
                    "sign_test_p": exact_two_sided_sign_test(increases, decreases),
                }
            )
    results = pd.DataFrame(rows)
    if results.empty:
        return empty
    results["bh_q_within_contrast"] = np.nan
    for _, indices in results.groupby(["dataset", "contrast"], sort=True).groups.items():
        results.loc[indices, "bh_q_within_contrast"] = benjamini_hochberg(
            results.loc[indices, "sign_test_p"]
        )
    return results.sort_values(["dataset", "contrast", "module"], kind="mergesort").reset_index(
        drop=True
    )
