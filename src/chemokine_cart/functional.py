"""Analysis utilities for functional chemokine assays.

The functions in this module operate on long-format observations and make the
biological replicate, rather than the well, image, field, or cell, the unit of
analysis.  They are suitable for chemotaxis, retention, egress, CXCL16-form,
and receptor-matched product comparisons without imposing a biological
direction on the result.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from scipy import stats

SUPPORTED_ASSAYS = frozenset({"chemotaxis", "retention", "egress", "cxcl16"})
SUPPORTED_CXCL16_FORMS = frozenset({"soluble", "membrane"})

_TECHNICAL_METADATA_COLUMNS = frozenset(
    {
        "technical_replicate",
        "well",
        "well_id",
        "field",
        "field_id",
        "cell_id",
        "event_id",
        "plate",
        "plate_id",
        "image",
        "image_id",
    }
)


@dataclass(frozen=True)
class FunctionalQC:
    """Structured quality-control result for a functional-assay table."""

    passed: bool
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    n_rows: int
    n_biological_replicates: int
    n_missing_response: int
    n_nonfinite_response: int
    n_exact_duplicate_rows: int
    n_repeated_unit_condition_rows: int

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable representation."""

        return asdict(self)


@dataclass(frozen=True)
class ComparisonResult:
    """Result of a two-group comparison at the biological-replicate level.

    ``mean_difference`` and its interval are always group B minus group A.
    For unpaired data, ``standardized_effect`` is Hedges' g.  For paired data,
    it is the small-sample-corrected standardized mean paired difference.
    """

    group_a: str
    group_b: str
    paired: bool
    n_a: int
    n_b: int
    n_pairs: int | None
    mean_a: float
    mean_b: float
    median_a: float
    median_b: float
    mean_difference: float
    mean_difference_ci_low: float
    mean_difference_ci_high: float
    standardized_effect: float
    standardized_effect_ci_low: float
    standardized_effect_ci_high: float
    test: str
    statistic: float
    p_value: float
    bootstrap_iterations: int
    confidence_level: float
    dropped_incomplete_pairs: int
    notes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable representation."""

        return asdict(self)


def _require_columns(data: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def _as_numeric(data: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(data[column], errors="coerce")
    invalid = data[column].notna() & values.isna()
    if invalid.any():
        raise ValueError(f"Column '{column}' contains {int(invalid.sum())} non-numeric value(s).")
    return values.astype(float)


def _infer_design_columns(
    data: pd.DataFrame,
    *,
    response_col: str,
    biological_replicate_col: str,
    technical_replicate_col: str | None,
) -> list[str]:
    excluded = {
        response_col,
        biological_replicate_col,
        *(_TECHNICAL_METADATA_COLUMNS & set(data.columns)),
    }
    if technical_replicate_col is not None:
        excluded.add(technical_replicate_col)
    return [column for column in data.columns if column not in excluded]


def validate_functional_data(
    data: pd.DataFrame,
    *,
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    assay_col: str = "assay",
    technical_replicate_col: str | None = "technical_replicate",
    design_cols: Sequence[str] | None = None,
    allowed_assays: Sequence[str] | None = tuple(sorted(SUPPORTED_ASSAYS)),
    cxcl16_form_col: str = "cxcl16_form",
    allowed_cxcl16_forms: Sequence[str] = tuple(sorted(SUPPORTED_CXCL16_FORMS)),
    response_bounds: tuple[float | None, float | None] | None = None,
) -> FunctionalQC:
    """Validate a long-format functional-assay table without altering it.

    Replicate identifiers and responses are required.  Assay labels are
    checked only when ``assay_col`` is present.  Bounds are optional because
    valid response scales differ between counts, fractions, percentages, and
    normalized indices.
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")

    errors: list[str] = []
    warnings: list[str] = []
    required = [biological_replicate_col, response_col]
    missing = [column for column in required if column not in data.columns]
    if missing:
        return FunctionalQC(
            passed=False,
            errors=(f"Missing required columns: {', '.join(missing)}",),
            warnings=(),
            n_rows=len(data),
            n_biological_replicates=0,
            n_missing_response=0,
            n_nonfinite_response=0,
            n_exact_duplicate_rows=int(data.duplicated().sum()),
            n_repeated_unit_condition_rows=0,
        )

    if data.empty:
        errors.append("The input table contains no observations.")

    missing_units = int(data[biological_replicate_col].isna().sum())
    if missing_units:
        errors.append(f"{missing_units} row(s) have no biological-replicate identifier.")

    numeric = pd.to_numeric(data[response_col], errors="coerce")
    missing_response = int(data[response_col].isna().sum())
    nonnumeric = int((data[response_col].notna() & numeric.isna()).sum())
    nonfinite = int((numeric.notna() & ~np.isfinite(numeric)).sum())
    if missing_response:
        errors.append(f"{missing_response} row(s) have a missing response.")
    if nonnumeric:
        errors.append(f"{nonnumeric} row(s) have a non-numeric response.")
    if nonfinite:
        errors.append(f"{nonfinite} row(s) have a non-finite response.")

    if response_bounds is not None:
        lower, upper = response_bounds
        finite = numeric[np.isfinite(numeric)]
        if lower is not None:
            below = int((finite < lower).sum())
            if below:
                errors.append(f"{below} response value(s) are below {lower}.")
        if upper is not None:
            above = int((finite > upper).sum())
            if above:
                errors.append(f"{above} response value(s) are above {upper}.")

    if assay_col in data.columns:
        missing_assays = int(data[assay_col].isna().sum())
        if missing_assays:
            errors.append(f"{missing_assays} row(s) have no assay label.")
        if allowed_assays is not None:
            observed = set(data[assay_col].dropna().astype(str).str.lower())
            unexpected = sorted(observed - {str(value).lower() for value in allowed_assays})
            if unexpected:
                errors.append(f"Unsupported assay label(s): {', '.join(unexpected)}")
    else:
        warnings.append(f"Column '{assay_col}' is absent; assay-specific labels were not checked.")

    if cxcl16_form_col in data.columns:
        observed_forms = set(data[cxcl16_form_col].dropna().astype(str).str.lower().str.strip())
        allowed_forms = {str(value).lower() for value in allowed_cxcl16_forms}
        unexpected_forms = sorted(observed_forms - allowed_forms)
        if unexpected_forms:
            errors.append("Unsupported CXCL16 form label(s): " + ", ".join(unexpected_forms))

    exact_duplicates = int(data.duplicated().sum())
    if exact_duplicates:
        warnings.append(f"{exact_duplicates} exact duplicate row(s) were detected.")

    if design_cols is None:
        design = _infer_design_columns(
            data,
            response_col=response_col,
            biological_replicate_col=biological_replicate_col,
            technical_replicate_col=technical_replicate_col,
        )
    else:
        design = list(dict.fromkeys(design_cols))
        missing_design = [column for column in design if column not in data.columns]
        if missing_design:
            errors.append(f"Missing design columns: {', '.join(missing_design)}")
            design = [column for column in design if column in data.columns]

    unit_condition_cols = [biological_replicate_col, *design]
    repeated = int(data.duplicated(unit_condition_cols, keep=False).sum())
    if repeated:
        if technical_replicate_col and technical_replicate_col in data.columns:
            warnings.append(
                f"{repeated} row(s) share a biological unit and condition and will be "
                "collapsed before inference."
            )
            technical_key = [*unit_condition_cols, technical_replicate_col]
            duplicate_technical_ids = int(data.duplicated(technical_key, keep=False).sum())
            if duplicate_technical_ids:
                warnings.append(
                    f"{duplicate_technical_ids} row(s) repeat a technical-replicate "
                    "identifier within the same biological unit and condition."
                )
        else:
            warnings.append(
                f"{repeated} row(s) repeat a biological unit and condition without an "
                "explicit technical-replicate column; they will be collapsed before "
                "inference."
            )

    n_units = int(data[biological_replicate_col].dropna().nunique())
    if 0 < n_units < 2:
        warnings.append(
            "Fewer than two biological replicates are present; dispersion and "
            "inferential statistics are not estimable."
        )

    return FunctionalQC(
        passed=not errors,
        errors=tuple(errors),
        warnings=tuple(warnings),
        n_rows=len(data),
        n_biological_replicates=n_units,
        n_missing_response=missing_response,
        n_nonfinite_response=nonfinite + nonnumeric,
        n_exact_duplicate_rows=exact_duplicates,
        n_repeated_unit_condition_rows=repeated,
    )


def validate_protein_data(
    data: pd.DataFrame,
    *,
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    analyte_col: str = "analyte",
    form_col: str = "molecular_form",
    unit_col: str = "unit",
    assay_col: str = "assay",
    technical_replicate_col: str | None = "technical_replicate",
    allowed_forms: Sequence[str] = tuple(sorted(SUPPORTED_CXCL16_FORMS)),
    response_bounds: tuple[float | None, float | None] | None = None,
) -> FunctionalQC:
    """Validate long-format soluble and membrane protein measurements.

    Units are required and must be internally consistent within each
    analyte-form-assay stratum.  Different forms and assays may retain their
    native units because this function does not compare them numerically.
    """

    base = validate_functional_data(
        data,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        assay_col=assay_col,
        technical_replicate_col=technical_replicate_col,
        allowed_assays=None,
        cxcl16_form_col=form_col,
        allowed_cxcl16_forms=allowed_forms,
        response_bounds=response_bounds,
    )
    errors = list(base.errors)
    warnings = list(base.warnings)
    required = [analyte_col, form_col, unit_col, assay_col]
    missing = [column for column in required if column not in data.columns]
    if missing:
        errors.append(f"Missing protein metadata columns: {', '.join(missing)}")
    else:
        for column in required:
            n_missing = int(data[column].isna().sum())
            if n_missing:
                errors.append(f"{n_missing} row(s) have no value in '{column}'.")
        complete = data.dropna(subset=required)
        if not complete.empty:
            unit_counts = complete.groupby(
                [analyte_col, form_col, assay_col], dropna=False, observed=True
            )[unit_col].nunique()
            inconsistent = unit_counts[unit_counts > 1]
            if not inconsistent.empty:
                errors.append(
                    f"{len(inconsistent)} analyte-form-assay stratum/strata use "
                    "multiple response units."
                )

    return FunctionalQC(
        passed=not errors,
        errors=tuple(dict.fromkeys(errors)),
        warnings=tuple(dict.fromkeys(warnings)),
        n_rows=base.n_rows,
        n_biological_replicates=base.n_biological_replicates,
        n_missing_response=base.n_missing_response,
        n_nonfinite_response=base.n_nonfinite_response,
        n_exact_duplicate_rows=base.n_exact_duplicate_rows,
        n_repeated_unit_condition_rows=base.n_repeated_unit_condition_rows,
    )


def aggregate_biological_replicates(
    data: pd.DataFrame,
    *,
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    design_cols: Sequence[str] | None = None,
    method: Literal["mean", "median"] = "mean",
) -> pd.DataFrame:
    """Collapse subsamples to one response per biological unit and condition.

    When ``design_cols`` is omitted, non-identifier columns are retained as
    experimental factors while conventional technical metadata (well, plate,
    field, image, cell, and event identifiers) are excluded from grouping.
    """

    _require_columns(data, [biological_replicate_col, response_col])
    if method not in {"mean", "median"}:
        raise ValueError("method must be 'mean' or 'median'")

    working = data.copy()
    working[response_col] = _as_numeric(working, response_col)
    if working[biological_replicate_col].isna().any():
        raise ValueError("Biological-replicate identifiers must not be missing.")
    if working[response_col].isna().any() or (~np.isfinite(working[response_col])).any():
        raise ValueError("Responses must be complete and finite before aggregation.")

    if design_cols is None:
        design = _infer_design_columns(
            working,
            response_col=response_col,
            biological_replicate_col=biological_replicate_col,
            technical_replicate_col=technical_replicate_col,
        )
    else:
        design = list(dict.fromkeys(design_cols))
        _require_columns(working, design)

    group_cols = list(dict.fromkeys([biological_replicate_col, *design]))
    grouped = working.groupby(group_cols, dropna=False, observed=True)[response_col]
    if method == "mean":
        aggregate = grouped.mean()
    else:
        aggregate = grouped.median()

    result = aggregate.rename(response_col).reset_index()
    counts = grouped.size().rename("n_subsamples").reset_index()
    result = result.merge(counts, on=group_cols, how="left", validate="one_to_one")
    result.attrs["analysis_unit"] = biological_replicate_col
    result.attrs["aggregation"] = method
    return result


def summarize_groups(
    data: pd.DataFrame,
    *,
    group_cols: str | Sequence[str],
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    technical_aggregation: Literal["mean", "median"] = "mean",
) -> pd.DataFrame:
    """Summarize responses after collapsing technical subsamples."""

    groups = [group_cols] if isinstance(group_cols, str) else list(group_cols)
    _require_columns(data, groups)
    aggregated = aggregate_biological_replicates(
        data,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        technical_replicate_col=technical_replicate_col,
        design_cols=groups,
        method=technical_aggregation,
    )
    summary = (
        aggregated.groupby(groups, dropna=False, observed=True)[response_col]
        .agg(
            n="count",
            mean="mean",
            sd="std",
            median="median",
            minimum="min",
            maximum="max",
        )
        .reset_index()
    )
    summary["sem"] = summary["sd"] / np.sqrt(summary["n"])
    quantiles = (
        aggregated.groupby(groups, dropna=False, observed=True)[response_col]
        .quantile([0.25, 0.75])
        .unstack()
        .rename(columns={0.25: "q25", 0.75: "q75"})
        .reset_index()
    )
    return summary.merge(quantiles, on=groups, how="left", validate="one_to_one")


def _standardized_unpaired(group_a: np.ndarray, group_b: np.ndarray) -> float:
    n_a, n_b = len(group_a), len(group_b)
    if n_a < 2 or n_b < 2:
        return float("nan")
    pooled_variance = (
        (n_a - 1) * np.var(group_a, ddof=1) + (n_b - 1) * np.var(group_b, ddof=1)
    ) / (n_a + n_b - 2)
    if not np.isfinite(pooled_variance) or pooled_variance <= 0:
        return float("nan")
    cohen_d = (np.mean(group_b) - np.mean(group_a)) / np.sqrt(pooled_variance)
    correction = 1.0 - 3.0 / (4.0 * (n_a + n_b) - 9.0)
    return float(correction * cohen_d)


def _standardized_paired(differences: np.ndarray) -> float:
    n_pairs = len(differences)
    if n_pairs < 2:
        return float("nan")
    sd = np.std(differences, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return float("nan")
    correction = 1.0 - 3.0 / (4.0 * n_pairs - 5.0)
    return float(correction * np.mean(differences) / sd)


def _bootstrap_effects(
    group_a: np.ndarray,
    group_b: np.ndarray,
    *,
    paired: bool,
    iterations: int,
    confidence_level: float,
    seed: int | None,
) -> tuple[float, float, float, float]:
    if iterations < 100:
        raise ValueError("bootstrap_iterations must be at least 100")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between 0 and 1")

    rng = np.random.default_rng(seed)
    differences = np.empty(iterations, dtype=float)
    standardized = np.empty(iterations, dtype=float)
    if paired:
        paired_differences = group_b - group_a
        size = len(paired_differences)
        for index in range(iterations):
            sampled = paired_differences[rng.integers(0, size, size=size)]
            differences[index] = np.mean(sampled)
            standardized[index] = _standardized_paired(sampled)
    else:
        n_a, n_b = len(group_a), len(group_b)
        for index in range(iterations):
            sample_a = group_a[rng.integers(0, n_a, size=n_a)]
            sample_b = group_b[rng.integers(0, n_b, size=n_b)]
            differences[index] = np.mean(sample_b) - np.mean(sample_a)
            standardized[index] = _standardized_unpaired(sample_a, sample_b)

    alpha = (1.0 - confidence_level) / 2.0
    diff_low, diff_high = np.quantile(differences, [alpha, 1.0 - alpha])
    finite_standardized = standardized[np.isfinite(standardized)]
    if finite_standardized.size:
        effect_low, effect_high = np.quantile(finite_standardized, [alpha, 1.0 - alpha])
    else:
        effect_low = effect_high = float("nan")
    return float(diff_low), float(diff_high), float(effect_low), float(effect_high)


def _paired_values(
    aggregated: pd.DataFrame,
    *,
    group_col: str,
    group_a: object,
    group_b: object,
    pair_col: str,
    biological_replicate_col: str,
    response_col: str,
) -> tuple[pd.DataFrame, int]:
    selected = aggregated[aggregated[group_col].isin([group_a, group_b])]
    _validate_pair_mapping(
        selected,
        pair_col=pair_col,
        identity_columns=(biological_replicate_col,),
    )
    counts = selected.groupby([pair_col, group_col], dropna=False, observed=True).size()
    if (counts > 1).any():
        raise ValueError(
            "More than one biological-unit response exists per pair and group. "
            "Subset to one assay context or revise design_cols."
        )
    wide = selected.pivot(index=pair_col, columns=group_col, values=response_col)
    total_pairs = len(wide)
    for group in (group_a, group_b):
        if group not in wide.columns:
            wide[group] = np.nan
    complete = wide[[group_a, group_b]].dropna()
    return complete, total_pairs - len(complete)


def _validate_pair_mapping(
    data: pd.DataFrame,
    *,
    pair_col: str,
    identity_columns: Sequence[str],
) -> None:
    """Require every pair identifier to resolve to one biological identity."""

    _require_columns(data, [pair_col, *identity_columns])
    missing_pair = data[pair_col].isna() | data[pair_col].astype("string").str.strip().eq(
        ""
    ).fillna(True)
    if missing_pair.any():
        raise ValueError(f"Paired analysis requires a non-missing {pair_col} for every row.")

    identities = list(dict.fromkeys(identity_columns))
    for identity_col in identities:
        missing_identity = data[identity_col].isna() | data[identity_col].astype(
            "string"
        ).str.strip().eq("").fillna(True)
        if missing_identity.any():
            raise ValueError(
                f"Paired analysis requires a non-missing {identity_col} for every row."
            )
        identities_per_pair = data.groupby(pair_col, dropna=False, observed=True)[
            identity_col
        ].nunique(dropna=False)
        invalid_pairs = identities_per_pair[identities_per_pair != 1]
        if not invalid_pairs.empty:
            examples = ", ".join(map(str, invalid_pairs.index[:5]))
            raise ValueError(
                f"Each {pair_col} must link exactly one {identity_col}; invalid pair(s): {examples}"
            )

    pairs_per_identity = data.groupby(identities, dropna=False, observed=True)[pair_col].nunique(
        dropna=False
    )
    duplicated_identities = pairs_per_identity[pairs_per_identity != 1]
    if not duplicated_identities.empty:
        examples = ", ".join(
            "/".join(map(str, key if isinstance(key, tuple) else (key,)))
            for key in duplicated_identities.index[:5]
        )
        raise ValueError(
            "Each biological identity must link exactly one "
            f"{pair_col}; invalid biological identity/identities: {examples}"
        )

    for cohort_identity in (
        column for column in identities if column in {"patient_id", "donor_id"}
    ):
        pairs_per_cohort_unit = data.groupby(cohort_identity, dropna=False, observed=True)[
            pair_col
        ].nunique(dropna=False)
        duplicated_cohort_units = pairs_per_cohort_unit[pairs_per_cohort_unit != 1]
        if not duplicated_cohort_units.empty:
            examples = ", ".join(map(str, duplicated_cohort_units.index[:5]))
            raise ValueError(
                f"Each {cohort_identity} must link exactly one {pair_col} within a paired "
                f"comparison; invalid {cohort_identity}(s): {examples}"
            )


def compare_groups(
    data: pd.DataFrame,
    *,
    group_col: str,
    group_a: object,
    group_b: object,
    paired: bool,
    pair_col: str | None = None,
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    design_cols: Sequence[str] | None = None,
    technical_aggregation: Literal["mean", "median"] = "mean",
    test: str | None = None,
    alternative: Literal["two-sided", "less", "greater"] = "two-sided",
    bootstrap_iterations: int = 5000,
    confidence_level: float = 0.95,
    seed: int | None = 0,
) -> ComparisonResult:
    """Compare two groups after collapsing subsamples within biological units.

    For a paired design, ``pair_col`` identifies matched biological units and
    defaults to ``biological_replicate_col``.  Missing members of a pair are
    excluded and counted.  The contrast is always group B minus group A.
    """

    _require_columns(data, [group_col, biological_replicate_col, response_col])
    if group_a == group_b:
        raise ValueError("group_a and group_b must be different")

    selected = data[data[group_col].isin([group_a, group_b])].copy()
    if selected.empty:
        raise ValueError("Neither requested group is present in the input table.")
    missing_groups = [
        str(group) for group in (group_a, group_b) if group not in set(selected[group_col].dropna())
    ]
    if missing_groups:
        raise ValueError(f"Requested group(s) not found: {', '.join(missing_groups)}")

    pairing = pair_col or biological_replicate_col
    if paired:
        _require_columns(selected, [pairing])
        pair_identity_columns = [biological_replicate_col]
        pair_identity_columns.extend(
            column for column in ("patient_id", "donor_id") if column in selected.columns
        )
        _validate_pair_mapping(
            selected,
            pair_col=pairing,
            identity_columns=pair_identity_columns,
        )

    if design_cols is None:
        design = _infer_design_columns(
            selected,
            response_col=response_col,
            biological_replicate_col=biological_replicate_col,
            technical_replicate_col=technical_replicate_col,
        )
    else:
        design = list(dict.fromkeys(design_cols))
    for essential in (group_col, pairing if paired else None):
        if essential is not None and essential not in design:
            design.append(essential)

    aggregated = aggregate_biological_replicates(
        selected,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        technical_replicate_col=technical_replicate_col,
        design_cols=design,
        method=technical_aggregation,
    )

    notes: list[str] = []
    dropped_pairs = 0
    if paired:
        wide, dropped_pairs = _paired_values(
            aggregated,
            group_col=group_col,
            group_a=group_a,
            group_b=group_b,
            pair_col=pairing,
            biological_replicate_col=biological_replicate_col,
            response_col=response_col,
        )
        if wide.empty:
            raise ValueError("No complete biological pairs are available.")
        if len(wide) < 2:
            raise ValueError("At least two complete biological pairs are required.")
        values_a = wide[group_a].to_numpy(dtype=float)
        values_b = wide[group_b].to_numpy(dtype=float)
        if dropped_pairs:
            notes.append(f"Excluded {dropped_pairs} incomplete pair(s).")
        selected_test = test or "paired_t"
        if selected_test == "paired_t":
            test_result = stats.ttest_rel(values_b, values_a, alternative=alternative)
        elif selected_test == "wilcoxon":
            test_result = stats.wilcoxon(
                values_b, values_a, alternative=alternative, zero_method="wilcox"
            )
        else:
            raise ValueError("Paired test must be 'paired_t' or 'wilcoxon'.")
        standardized_effect = _standardized_paired(values_b - values_a)
        n_pairs: int | None = len(wide)
    else:
        group_counts = aggregated.groupby(
            [biological_replicate_col, group_col], dropna=False, observed=True
        ).size()
        if (group_counts > 1).any():
            raise ValueError(
                "More than one response exists per biological unit and group. "
                "Subset to one assay context or revise design_cols."
            )
        values_a = aggregated.loc[aggregated[group_col] == group_a, response_col].to_numpy(
            dtype=float
        )
        values_b = aggregated.loc[aggregated[group_col] == group_b, response_col].to_numpy(
            dtype=float
        )
        if len(values_a) < 2 or len(values_b) < 2:
            raise ValueError("Each group requires at least two biological replicates.")
        selected_test = test or "welch_t"
        if selected_test == "welch_t":
            test_result = stats.ttest_ind(
                values_b,
                values_a,
                equal_var=False,
                alternative=alternative,
            )
        elif selected_test == "mannwhitney":
            test_result = stats.mannwhitneyu(values_b, values_a, alternative=alternative)
        else:
            raise ValueError("Unpaired test must be 'welch_t' or 'mannwhitney'.")
        standardized_effect = _standardized_unpaired(values_a, values_b)
        n_pairs = None

    diff_ci_low, diff_ci_high, effect_ci_low, effect_ci_high = _bootstrap_effects(
        values_a,
        values_b,
        paired=paired,
        iterations=bootstrap_iterations,
        confidence_level=confidence_level,
        seed=seed,
    )
    return ComparisonResult(
        group_a=str(group_a),
        group_b=str(group_b),
        paired=paired,
        n_a=len(values_a),
        n_b=len(values_b),
        n_pairs=n_pairs,
        mean_a=float(np.mean(values_a)),
        mean_b=float(np.mean(values_b)),
        median_a=float(np.median(values_a)),
        median_b=float(np.median(values_b)),
        mean_difference=float(np.mean(values_b) - np.mean(values_a)),
        mean_difference_ci_low=diff_ci_low,
        mean_difference_ci_high=diff_ci_high,
        standardized_effect=standardized_effect,
        standardized_effect_ci_low=effect_ci_low,
        standardized_effect_ci_high=effect_ci_high,
        test=selected_test,
        statistic=float(test_result.statistic),
        p_value=float(test_result.pvalue),
        bootstrap_iterations=bootstrap_iterations,
        confidence_level=confidence_level,
        dropped_incomplete_pairs=dropped_pairs,
        notes=tuple(notes),
    )


def compare_cxcl16_forms(
    data: pd.DataFrame,
    *,
    soluble_label: object = "soluble",
    membrane_label: object = "membrane",
    form_col: str = "cxcl16_form",
    response_unit_col: str = "unit",
    assay_col: str = "assay",
    comparable_scale_confirmed: bool = False,
    paired: bool = True,
    **kwargs: object,
) -> ComparisonResult:
    """Compare membrane-associated with soluble CXCL16 on a common scale.

    The reported contrast is membrane minus soluble.  This wrapper does not
    permit a direct effect size between measurements such as supernatant
    concentration and cell-surface fluorescence.  The caller must explicitly
    confirm a common normalized scale, and both forms must have exactly the
    same response unit and assay label.  Use :func:`summarize_cxcl16_forms`
    when the forms were measured on different scales.
    """

    _require_columns(data, [form_col, response_unit_col, assay_col])
    selected = data[data[form_col].isin([soluble_label, membrane_label])]
    if not comparable_scale_confirmed:
        raise ValueError(
            "A cross-form effect size requires comparable_scale_confirmed=True "
            "after soluble and membrane CXCL16 have been placed on a common, "
            "biologically interpretable scale."
        )
    if selected.empty:
        raise ValueError("Neither requested CXCL16 form is present.")
    units = selected[response_unit_col].dropna().astype(str).unique()
    if selected[response_unit_col].isna().any() or len(units) != 1:
        raise ValueError(
            "Soluble and membrane CXCL16 cannot be compared across different or "
            "missing response units; use summarize_cxcl16_forms instead."
        )
    assays = selected[assay_col].dropna().astype(str).unique()
    if selected[assay_col].isna().any() or len(assays) != 1:
        raise ValueError(
            "Soluble and membrane CXCL16 effect sizes require one common assay; "
            "use summarize_cxcl16_forms for assay-specific descriptions."
        )

    return compare_groups(
        selected,
        group_col=form_col,
        group_a=soluble_label,
        group_b=membrane_label,
        paired=paired,
        **kwargs,
    )


def summarize_cxcl16_forms(
    data: pd.DataFrame,
    *,
    form_col: str = "cxcl16_form",
    response_unit_col: str = "unit",
    assay_col: str = "assay",
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    technical_aggregation: Literal["mean", "median"] = "mean",
) -> pd.DataFrame:
    """Describe each CXCL16 form separately without a cross-scale effect size.

    Results are stratified by form, response unit, and assay.  This is the safe
    default for soluble protein concentrations and membrane abundance because
    their numerical scales are generally not commensurable.
    """

    _require_columns(data, [form_col, response_unit_col, assay_col])
    return summarize_groups(
        data,
        group_cols=[form_col, response_unit_col, assay_col],
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        technical_replicate_col=technical_replicate_col,
        technical_aggregation=technical_aggregation,
    )


def compare_receptor_matched_product(
    data: pd.DataFrame,
    *,
    control_label: object = "identical_control",
    matched_label: object = "receptor_matched",
    product_col: str = "product",
    paired: bool = True,
    **kwargs: object,
) -> ComparisonResult:
    """Compare a receptor-matched product with its otherwise identical control.

    The reported contrast is receptor-matched minus identical control.  Product
    identity must be defined in the input metadata; the function cannot verify
    construct equivalence from response data alone.
    """

    return compare_groups(
        data,
        group_col=product_col,
        group_a=control_label,
        group_b=matched_label,
        paired=paired,
        **kwargs,
    )


def _new_axis(ax: Axes | None, *, figsize: tuple[float, float]) -> Axes:
    if ax is not None:
        return ax
    _, created = plt.subplots(figsize=figsize, constrained_layout=True)
    return created


def plot_dose_response(
    data: pd.DataFrame,
    *,
    dose_col: str = "dose",
    group_col: str | None = "product",
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    technical_aggregation: Literal["mean", "median"] = "mean",
    error: Literal["sd", "sem", "none"] = "sd",
    log_x: bool = False,
    response_label: str | None = None,
    ax: Axes | None = None,
) -> tuple[Axes, pd.DataFrame]:
    """Plot biological-replicate dose responses and descriptive summaries."""

    required = [dose_col, biological_replicate_col, response_col]
    if group_col is not None:
        required.append(group_col)
    _require_columns(data, required)
    if error not in {"sd", "sem", "none"}:
        raise ValueError("error must be 'sd', 'sem', or 'none'")

    working = data.copy()
    working[dose_col] = _as_numeric(working, dose_col)
    if working[dose_col].isna().any() or (~np.isfinite(working[dose_col])).any():
        raise ValueError("Dose values must be complete and finite.")
    if log_x and (working[dose_col] <= 0).any():
        raise ValueError("All doses must be positive when log_x=True.")

    design = [dose_col]
    if group_col is not None:
        design.append(group_col)
    aggregated = aggregate_biological_replicates(
        working,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        technical_replicate_col=technical_replicate_col,
        design_cols=design,
        method=technical_aggregation,
    )
    summary_groups = design
    summary = (
        aggregated.groupby(summary_groups, dropna=False, observed=True)[response_col]
        .agg(n="count", mean="mean", sd="std")
        .reset_index()
    )
    summary["sem"] = summary["sd"] / np.sqrt(summary["n"])

    axis = _new_axis(ax, figsize=(6.4, 4.2))
    if group_col is None:
        series = [("All", aggregated, summary)]
    else:
        series = []
        for label, group_data in aggregated.groupby(group_col, sort=False, observed=True):
            group_summary = summary[summary[group_col] == label]
            series.append((str(label), group_data, group_summary))

    for label, observations, group_summary in series:
        group_summary = group_summary.sort_values(dose_col)
        line = axis.plot(
            group_summary[dose_col],
            group_summary["mean"],
            marker="o",
            linewidth=1.8,
            label=label,
        )[0]
        color = line.get_color()
        axis.scatter(
            observations[dose_col],
            observations[response_col],
            color=color,
            s=24,
            alpha=0.45,
            linewidths=0,
        )
        if error != "none":
            axis.errorbar(
                group_summary[dose_col],
                group_summary["mean"],
                yerr=group_summary[error],
                fmt="none",
                ecolor=color,
                elinewidth=1.1,
                capsize=3,
            )

    if log_x:
        axis.set_xscale("log")
    axis.set_xlabel(dose_col.replace("_", " ").title())
    axis.set_ylabel(response_label or response_col.replace("_", " ").title())
    if group_col is not None:
        axis.legend(frameon=False, title=group_col.replace("_", " ").title())
    axis.spines[["top", "right"]].set_visible(False)
    return axis, summary


def plot_migration(
    data: pd.DataFrame,
    *,
    condition_col: str = "condition",
    group_col: str | None = "product",
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    technical_aggregation: Literal["mean", "median"] = "mean",
    response_label: str = "Migration response",
    ax: Axes | None = None,
) -> tuple[Axes, pd.DataFrame]:
    """Plot migration observations with mean and standard deviation."""

    design = [condition_col]
    if group_col is not None:
        design.append(group_col)
    _require_columns(data, [biological_replicate_col, response_col, *design])
    aggregated = aggregate_biological_replicates(
        data,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        technical_replicate_col=technical_replicate_col,
        design_cols=design,
        method=technical_aggregation,
    )
    summary = (
        aggregated.groupby(design, dropna=False, observed=True)[response_col]
        .agg(n="count", mean="mean", sd="std")
        .reset_index()
    )
    summary["sem"] = summary["sd"] / np.sqrt(summary["n"])

    axis = _new_axis(ax, figsize=(6.4, 4.2))
    conditions = list(pd.unique(aggregated[condition_col]))
    if group_col is None:
        groups: list[object | None] = [None]
    else:
        groups = list(pd.unique(aggregated[group_col]))
    width = 0.7 / max(len(groups), 1)

    for group_index, group in enumerate(groups):
        offset = (group_index - (len(groups) - 1) / 2) * width
        label = "All" if group is None else str(group)
        color = f"C{group_index % 10}"
        for condition_index, condition in enumerate(conditions):
            mask = aggregated[condition_col] == condition
            if group is not None:
                mask &= aggregated[group_col] == group
            values = aggregated.loc[mask, response_col].to_numpy(dtype=float)
            if values.size == 0:
                continue
            center = condition_index + offset
            point_offsets = np.linspace(-0.08, 0.08, values.size) if values.size > 1 else [0]
            axis.scatter(
                center + np.asarray(point_offsets),
                values,
                color=color,
                s=30,
                alpha=0.65,
                linewidths=0,
                label=label if condition_index == 0 else None,
            )
            axis.errorbar(
                center,
                np.mean(values),
                yerr=np.std(values, ddof=1) if values.size > 1 else np.nan,
                fmt="D",
                color=color,
                markeredgecolor="white",
                markeredgewidth=0.7,
                capsize=3,
                zorder=4,
            )

    axis.set_xticks(range(len(conditions)), [str(value) for value in conditions])
    axis.set_xlabel(condition_col.replace("_", " ").title())
    axis.set_ylabel(response_label)
    if group_col is not None:
        axis.legend(frameon=False, title=group_col.replace("_", " ").title())
    axis.spines[["top", "right"]].set_visible(False)
    return axis, summary


def plot_paired_comparison(
    data: pd.DataFrame,
    *,
    group_col: str,
    group_a: object,
    group_b: object,
    pair_col: str | None = None,
    biological_replicate_col: str = "biological_replicate",
    response_col: str = "response",
    technical_replicate_col: str | None = "technical_replicate",
    technical_aggregation: Literal["mean", "median"] = "mean",
    response_label: str | None = None,
    ax: Axes | None = None,
) -> tuple[Axes, pd.DataFrame]:
    """Plot complete biological pairs after technical-replicate aggregation."""

    pairing = pair_col or biological_replicate_col
    _require_columns(data, [group_col, pairing, biological_replicate_col, response_col])
    selected = data[data[group_col].isin([group_a, group_b])].copy()
    pair_identity_columns = [biological_replicate_col]
    pair_identity_columns.extend(
        column for column in ("patient_id", "donor_id") if column in selected.columns
    )
    _validate_pair_mapping(
        selected,
        pair_col=pairing,
        identity_columns=pair_identity_columns,
    )
    aggregated = aggregate_biological_replicates(
        selected,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
        technical_replicate_col=technical_replicate_col,
        design_cols=list(dict.fromkeys([group_col, pairing])),
        method=technical_aggregation,
    )
    wide, _ = _paired_values(
        aggregated,
        group_col=group_col,
        group_a=group_a,
        group_b=group_b,
        pair_col=pairing,
        biological_replicate_col=biological_replicate_col,
        response_col=response_col,
    )
    if wide.empty:
        raise ValueError("No complete biological pairs are available.")

    axis = _new_axis(ax, figsize=(4.5, 4.2))
    for _, row in wide.iterrows():
        axis.plot(
            [0, 1],
            [row[group_a], row[group_b]],
            color="0.72",
            linewidth=1,
            zorder=1,
        )
    axis.scatter(np.zeros(len(wide)), wide[group_a], color="C0", s=34, zorder=2, label=str(group_a))
    axis.scatter(np.ones(len(wide)), wide[group_b], color="C1", s=34, zorder=2, label=str(group_b))
    axis.scatter(
        [0, 1],
        [wide[group_a].mean(), wide[group_b].mean()],
        marker="D",
        color=["C0", "C1"],
        edgecolor="white",
        linewidth=0.8,
        s=62,
        zorder=3,
    )
    axis.set_xticks([0, 1], [str(group_a), str(group_b)])
    axis.set_xlim(-0.35, 1.35)
    axis.set_ylabel(response_label or response_col.replace("_", " ").title())
    axis.spines[["top", "right"]].set_visible(False)
    wide = wide.reset_index()
    wide["difference_b_minus_a"] = wide[group_b] - wide[group_a]
    return axis, wide


__all__ = [
    "ComparisonResult",
    "FunctionalQC",
    "SUPPORTED_ASSAYS",
    "SUPPORTED_CXCL16_FORMS",
    "aggregate_biological_replicates",
    "compare_cxcl16_forms",
    "compare_groups",
    "compare_receptor_matched_product",
    "plot_dose_response",
    "plot_migration",
    "plot_paired_comparison",
    "summarize_groups",
    "summarize_cxcl16_forms",
    "validate_functional_data",
    "validate_protein_data",
]
