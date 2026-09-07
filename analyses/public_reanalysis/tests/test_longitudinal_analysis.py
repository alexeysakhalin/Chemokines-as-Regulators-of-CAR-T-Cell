from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from longitudinal_extension.analysis import (
    analyze_contrasts,
    build_contrast_inclusion,
    leave_one_out,
    summarize_gse197268_product_strata,
    summarize_samples,
    validate_cells,
)
from longitudinal_extension.common import exact_two_sided_sign_test, holm_adjust


@pytest.mark.parametrize(
    ("up", "down", "expected"),
    [(7, 2, 0.1796875), (4, 1, 0.375), (4, 6, 0.75390625), (1, 6, 0.125), (1, 3, 0.625)],
)
def test_sign_test_matches_locked_results(up: int, down: int, expected: float) -> None:
    assert exact_two_sided_sign_test(up, down) == expected


def test_sign_test_all_zeros_is_missing() -> None:
    assert math.isnan(exact_two_sided_sign_test(0, 0))


def test_holm_three_test_family() -> None:
    observed = holm_adjust(pd.Series([0.1796875, 0.375, 0.75390625]))
    np.testing.assert_allclose(observed, [0.5390625, 0.75, 0.75390625])


def test_cells_must_be_raw_integer_counts(synthetic_cells: pd.DataFrame) -> None:
    synthetic_cells["cxcr6_umi"] = synthetic_cells["cxcr6_umi"].astype(float)
    synthetic_cells.loc[0, "cxcr6_umi"] = 0.5
    with pytest.raises(ValueError, match="raw integer"):
        validate_cells(synthetic_cells)


def test_duplicate_cells_are_rejected(synthetic_cells: pd.DataFrame) -> None:
    duplicate = pd.concat([synthetic_cells, synthetic_cells.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="unique"):
        validate_cells(duplicate)


def test_sample_aggregation_preserves_patient_unit(synthetic_cells: pd.DataFrame) -> None:
    samples = summarize_samples(synthetic_cells)
    assert len(samples) == 6
    assert not samples.duplicated(["dataset", "patient_id", "product", "stage"]).any()
    p1_d7 = samples.query("patient_id == 'P1' and stage == 'D7'").iloc[0]
    assert p1_d7["eligible_cells"] == 4
    assert p1_d7["cxcr6_positive_cells"] == 3
    assert p1_d7["cxcr6_fraction"] == 0.75


def test_biological_sample_has_one_patient_product_and_stage(
    synthetic_cells: pd.DataFrame,
) -> None:
    conflicting = synthetic_cells.iloc[[0]].copy()
    conflicting["cell_id"] = "conflicting-cell"
    conflicting["patient_id"] = "another-patient"
    with pytest.raises(ValueError, match="exactly one patient, product, and stage"):
        summarize_samples(pd.concat([synthetic_cells, conflicting], ignore_index=True))


def test_contrast_is_one_row_per_patient(
    synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    samples = summarize_samples(synthetic_cells)
    inclusion = build_contrast_inclusion(samples, synthetic_contrasts)
    assert len(inclusion) == 3
    assert inclusion["patient_id"].is_unique


def test_contrast_summary_uses_directions_not_cells(
    synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    _, results = analyze_contrasts(summarize_samples(synthetic_cells), synthetic_contrasts)
    row = results.iloc[0]
    assert row["n_pairs"] == 3
    assert row["n_increase"] == 1
    assert row["n_decrease"] == 1
    assert row["n_unchanged"] == 1
    assert row["sign_test_p"] == 1.0


def test_minimum_cell_rule_is_applied_before_summary(
    synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    synthetic_contrasts.loc[0, "min_cells"] = 5
    inclusion, results = analyze_contrasts(summarize_samples(synthetic_cells), synthetic_contrasts)
    assert not inclusion["included"].any()
    assert results.iloc[0]["n_pairs"] == 0


def test_leave_one_out_requires_five_pairs(
    synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    inclusion = build_contrast_inclusion(summarize_samples(synthetic_cells), synthetic_contrasts)
    assert leave_one_out(inclusion).empty


def test_gse197268_product_strata_remain_descriptive(
    synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    synthetic_cells["dataset"] = "GSE197268"
    synthetic_cells.loc[synthetic_cells["patient_id"].eq("P3"), "product"] = "Tisa-cel"
    synthetic_cells.loc[~synthetic_cells["patient_id"].eq("P3"), "product"] = "Axi-cel"
    synthetic_contrasts.loc[0, ["dataset", "contrast"]] = [
        "GSE197268",
        "D7-CART_minus_IP",
    ]
    inclusion = build_contrast_inclusion(summarize_samples(synthetic_cells), synthetic_contrasts)
    strata = summarize_gse197268_product_strata(inclusion)
    assert strata["subgroup"].tolist() == ["Axi-cel", "Tisa-cel"]
    assert strata["role"].eq("descriptive product stratum").all()
    assert strata["holm_p"].isna().all()
