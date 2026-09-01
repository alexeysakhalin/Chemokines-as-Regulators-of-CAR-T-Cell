from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chemokine_cart.functional import (
    aggregate_biological_replicates,
    compare_cxcl16_forms,
    compare_groups,
    compare_receptor_matched_product,
    plot_dose_response,
    plot_migration,
    plot_paired_comparison,
    summarize_cxcl16_forms,
    summarize_groups,
    validate_functional_data,
    validate_protein_data,
)


@pytest.fixture
def technical_replicate_data() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    values = {
        ("D1", "identical_control"): [10.0, 12.0],
        ("D1", "receptor_matched"): [16.0, 18.0],
        ("D2", "identical_control"): [20.0, 22.0],
        ("D2", "receptor_matched"): [24.0, 26.0],
        ("D3", "identical_control"): [30.0, 32.0],
        ("D3", "receptor_matched"): [33.0, 35.0],
    }
    for (donor, product), responses in values.items():
        for technical_index, response in enumerate(responses, start=1):
            rows.append(
                {
                    "assay": "chemotaxis",
                    "biological_replicate": donor,
                    "technical_replicate": f"T{technical_index}",
                    "condition": "CXCL10",
                    "product": product,
                    "response": response,
                }
            )
    return pd.DataFrame(rows)


def test_qc_recognizes_valid_long_format_and_repeated_subsamples(
    technical_replicate_data: pd.DataFrame,
) -> None:
    qc = validate_functional_data(technical_replicate_data)

    assert qc.passed
    assert qc.n_rows == 12
    assert qc.n_biological_replicates == 3
    assert qc.n_repeated_unit_condition_rows == 12
    assert any("collapsed before inference" in warning for warning in qc.warnings)


def test_qc_reports_invalid_assay_cxcl16_form_response_and_bounds() -> None:
    data = pd.DataFrame(
        {
            "assay": ["unknown", "cxcl16", "cxcl16"],
            "biological_replicate": ["D1", "D2", None],
            "cxcl16_form": ["soluble", "secreted", "membrane"],
            "response": [20.0, "not_numeric", 150.0],
        }
    )

    qc = validate_functional_data(data, response_bounds=(0.0, 100.0))

    assert not qc.passed
    assert any("Unsupported assay" in error for error in qc.errors)
    assert any("Unsupported CXCL16" in error for error in qc.errors)
    assert any("non-numeric" in error for error in qc.errors)
    assert any("above 100.0" in error for error in qc.errors)
    assert any("no biological-replicate" in error for error in qc.errors)


def test_protein_qc_keeps_forms_on_native_units_and_detects_mixed_units() -> None:
    data = pd.DataFrame(
        {
            "biological_replicate": ["D1", "D2", "D1", "D2"],
            "analyte": ["CXCL16"] * 4,
            "molecular_form": ["soluble", "soluble", "membrane", "membrane"],
            "unit": ["pg_per_ml", "pg_per_ml", "MFI", "MFI"],
            "assay": ["immunoassay", "immunoassay", "flow_cytometry", "flow_cytometry"],
            "response": [40.0, 35.0, 800.0, 920.0],
        }
    )

    assert validate_protein_data(data).passed

    inconsistent = data.copy()
    inconsistent.loc[1, "unit"] = "ng_per_ml"
    qc = validate_protein_data(inconsistent)
    assert not qc.passed
    assert any("multiple response units" in error for error in qc.errors)


def test_aggregation_makes_biological_replicate_the_analysis_unit(
    technical_replicate_data: pd.DataFrame,
) -> None:
    aggregated = aggregate_biological_replicates(
        technical_replicate_data,
        design_cols=["assay", "condition", "product"],
    )

    assert len(aggregated) == 6
    assert set(aggregated["n_subsamples"]) == {2}
    d1_control = aggregated.loc[
        (aggregated["biological_replicate"] == "D1")
        & (aggregated["product"] == "identical_control"),
        "response",
    ].iat[0]
    assert d1_control == pytest.approx(11.0)


def test_group_summary_counts_biological_not_technical_replicates(
    technical_replicate_data: pd.DataFrame,
) -> None:
    summary = summarize_groups(technical_replicate_data, group_cols="product")

    assert set(summary["n"]) == {3}
    control_mean = summary.loc[summary["product"] == "identical_control", "mean"].iat[0]
    assert control_mean == pytest.approx(21.0)


def test_paired_comparison_uses_matched_biological_units(
    technical_replicate_data: pd.DataFrame,
) -> None:
    incomplete = pd.concat(
        [
            technical_replicate_data,
            pd.DataFrame(
                [
                    {
                        "assay": "chemotaxis",
                        "biological_replicate": "D4",
                        "technical_replicate": "T1",
                        "condition": "CXCL10",
                        "product": "identical_control",
                        "response": 15.0,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    result = compare_groups(
        incomplete,
        group_col="product",
        group_a="identical_control",
        group_b="receptor_matched",
        paired=True,
        bootstrap_iterations=400,
        seed=19,
    )

    assert result.n_pairs == 3
    assert result.n_a == result.n_b == 3
    assert result.dropped_incomplete_pairs == 1
    assert result.mean_a == pytest.approx(21.0)
    assert result.mean_b == pytest.approx(25.3333333333)
    assert result.mean_difference == pytest.approx(4.3333333333)
    assert result.mean_difference_ci_low <= result.mean_difference
    assert result.mean_difference_ci_high >= result.mean_difference
    assert result.standardized_effect > 0


def test_pair_id_cannot_link_different_biological_units() -> None:
    mismatched = pd.DataFrame(
        {
            "biological_replicate": ["D1", "D2", "D3", "D4"],
            "pair_id": ["PAIR_1", "PAIR_1", "PAIR_2", "PAIR_2"],
            "product": [
                "identical_control",
                "receptor_matched",
                "identical_control",
                "receptor_matched",
            ],
            "response": [1.0, 2.0, 3.0, 4.0],
        }
    )

    with pytest.raises(ValueError, match="must link exactly one biological_replicate"):
        compare_groups(
            mismatched,
            group_col="product",
            group_a="identical_control",
            group_b="receptor_matched",
            paired=True,
            pair_col="pair_id",
            bootstrap_iterations=200,
            seed=8,
        )


def test_pair_id_cannot_link_different_patients_with_reused_replicate_label() -> None:
    mismatched = pd.DataFrame(
        {
            "patient_id": ["P1", "P2", "P3", "P4"],
            "biological_replicate": ["R1", "R1", "R2", "R2"],
            "pair_id": ["PAIR_1", "PAIR_1", "PAIR_2", "PAIR_2"],
            "product": [
                "identical_control",
                "receptor_matched",
                "identical_control",
                "receptor_matched",
            ],
            "response": [1.0, 2.0, 3.0, 4.0],
        }
    )

    with pytest.raises(ValueError, match="must link exactly one patient_id"):
        compare_groups(
            mismatched,
            group_col="product",
            group_a="identical_control",
            group_b="receptor_matched",
            paired=True,
            pair_col="pair_id",
            design_cols=["product", "pair_id"],
            bootstrap_iterations=200,
            seed=9,
        )


def test_one_biological_identity_cannot_be_split_across_pair_ids() -> None:
    duplicated_unit = pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P1", "P1"],
            "biological_replicate": ["R1", "R1", "R1", "R1"],
            "pair_id": ["PAIR_1", "PAIR_1", "PAIR_2", "PAIR_2"],
            "product": [
                "identical_control",
                "receptor_matched",
                "identical_control",
                "receptor_matched",
            ],
            "response": [1.0, 2.0, 1.5, 2.5],
        }
    )

    with pytest.raises(ValueError, match="biological identity must link exactly one pair_id"):
        compare_groups(
            duplicated_unit,
            group_col="product",
            group_a="identical_control",
            group_b="receptor_matched",
            paired=True,
            pair_col="pair_id",
            design_cols=["product", "pair_id"],
            bootstrap_iterations=200,
            seed=10,
        )


def test_one_patient_cannot_contribute_multiple_independent_pair_ids() -> None:
    duplicated_patient = pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P1", "P1"],
            "biological_replicate": ["R1", "R1", "R2", "R2"],
            "pair_id": ["PAIR_1", "PAIR_1", "PAIR_2", "PAIR_2"],
            "product": [
                "identical_control",
                "receptor_matched",
                "identical_control",
                "receptor_matched",
            ],
            "response": [1.0, 2.0, 1.5, 2.5],
        }
    )

    with pytest.raises(ValueError, match="patient_id must link exactly one pair_id"):
        compare_groups(
            duplicated_patient,
            group_col="product",
            group_a="identical_control",
            group_b="receptor_matched",
            paired=True,
            pair_col="pair_id",
            design_cols=["product", "pair_id"],
            bootstrap_iterations=200,
            seed=13,
        )


@pytest.mark.parametrize("identity_column", ["biological_replicate", "patient_id"])
def test_pairing_rejects_blank_biological_identity(identity_column: str) -> None:
    blank_identity = pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P2", "P2"],
            "biological_replicate": ["R1", "R1", "R2", "R2"],
            "pair_id": ["PAIR_1", "PAIR_1", "PAIR_2", "PAIR_2"],
            "product": [
                "identical_control",
                "receptor_matched",
                "identical_control",
                "receptor_matched",
            ],
            "response": [1.0, 2.0, 3.0, 4.0],
        }
    )
    blank_identity.loc[0, identity_column] = ""

    with pytest.raises(ValueError, match=f"non-missing {identity_column}"):
        compare_groups(
            blank_identity,
            group_col="product",
            group_a="identical_control",
            group_b="receptor_matched",
            paired=True,
            pair_col="pair_id",
            bootstrap_iterations=200,
            seed=11,
        )


def test_pairing_rejects_blank_categorical_pair_id() -> None:
    blank_pair = pd.DataFrame(
        {
            "biological_replicate": ["R1", "R1", "R2", "R2"],
            "pair_id": pd.Categorical(["", "", "PAIR_2", "PAIR_2"]),
            "product": [
                "identical_control",
                "receptor_matched",
                "identical_control",
                "receptor_matched",
            ],
            "response": [1.0, 2.0, 3.0, 4.0],
        }
    )

    with pytest.raises(ValueError, match="non-missing pair_id"):
        compare_groups(
            blank_pair,
            group_col="product",
            group_a="identical_control",
            group_b="receptor_matched",
            paired=True,
            pair_col="pair_id",
            bootstrap_iterations=200,
            seed=12,
        )


def test_unpaired_comparison_does_not_count_wells_as_independent(
    technical_replicate_data: pd.DataFrame,
) -> None:
    result = compare_groups(
        technical_replicate_data,
        group_col="product",
        group_a="identical_control",
        group_b="receptor_matched",
        paired=False,
        bootstrap_iterations=300,
        seed=3,
    )

    assert result.n_a == 3
    assert result.n_b == 3
    assert result.n_pairs is None
    assert result.test == "welch_t"
    assert result.mean_difference > 0
    assert np.isfinite(result.p_value)


def test_comparison_requires_single_assay_context() -> None:
    data = pd.DataFrame(
        {
            "assay": [
                "chemotaxis",
                "retention",
                "chemotaxis",
                "retention",
            ]
            * 2,
            "biological_replicate": np.repeat(["D1", "D2"], 4),
            "product": ["A", "A", "B", "B"] * 2,
            "response": np.arange(8, dtype=float),
        }
    )

    with pytest.raises(ValueError, match="Subset to one assay context"):
        compare_groups(
            data,
            group_col="product",
            group_a="A",
            group_b="B",
            paired=False,
            bootstrap_iterations=100,
        )


def test_cxcl16_wrapper_reports_membrane_minus_soluble() -> None:
    data = pd.DataFrame(
        {
            "assay": ["cxcl16"] * 6,
            "biological_replicate": ["D1", "D1", "D2", "D2", "D3", "D3"],
            "cxcl16_form": ["soluble", "membrane"] * 3,
            "unit": ["normalized_common_scale"] * 6,
            "response": [1.0, 3.0, 2.0, 5.0, 4.0, 6.0],
        }
    )

    result = compare_cxcl16_forms(
        data,
        comparable_scale_confirmed=True,
        bootstrap_iterations=250,
        seed=4,
    )

    assert result.group_a == "soluble"
    assert result.group_b == "membrane"
    assert result.mean_difference == pytest.approx(7.0 / 3.0)


def test_cxcl16_cross_scale_effect_is_rejected_and_forms_are_summarized() -> None:
    data = pd.DataFrame(
        {
            "assay": ["cxcl16"] * 8,
            "biological_replicate": ["D1", "D2", "D3", "D4"] * 2,
            "cxcl16_form": ["soluble"] * 4 + ["membrane"] * 4,
            "unit": ["pg_per_ml"] * 4 + ["median_fluorescence"] * 4,
            "response": [30.0, 42.0, 35.0, 39.0, 800.0, 920.0, 850.0, 880.0],
        }
    )

    with pytest.raises(ValueError, match="different or missing response units"):
        compare_cxcl16_forms(data, comparable_scale_confirmed=True)

    summary = summarize_cxcl16_forms(data)
    assert len(summary) == 2
    assert set(summary["unit"]) == {"pg_per_ml", "median_fluorescence"}
    assert set(summary["n"]) == {4}


def test_product_wrapper_reports_matched_minus_identical_control(
    technical_replicate_data: pd.DataFrame,
) -> None:
    result = compare_receptor_matched_product(
        technical_replicate_data,
        bootstrap_iterations=250,
        seed=9,
    )

    assert result.n_pairs == 3
    assert result.mean_difference == pytest.approx(13.0 / 3.0)


def test_dose_response_plot_summarizes_biological_replicates() -> None:
    rows = []
    for donor, base in [("D1", 1.0), ("D2", 2.0), ("D3", 3.0)]:
        for dose in (1.0, 10.0):
            for technical in (1, 2):
                rows.append(
                    {
                        "biological_replicate": donor,
                        "technical_replicate": technical,
                        "product": "receptor_matched",
                        "dose": dose,
                        "response": base + np.log10(dose) + 0.1 * technical,
                    }
                )
    data = pd.DataFrame(rows)

    axis, summary = plot_dose_response(data, log_x=True)

    assert axis.get_xscale() == "log"
    assert set(summary["n"]) == {3}
    assert len(summary) == 2
    plt.close(axis.figure)


def test_migration_and_paired_plots_use_aggregated_units(
    technical_replicate_data: pd.DataFrame,
) -> None:
    migration_axis, migration_summary = plot_migration(
        technical_replicate_data,
        condition_col="condition",
        response_label="Migrated cells",
    )
    paired_axis, paired_table = plot_paired_comparison(
        technical_replicate_data,
        group_col="product",
        group_a="identical_control",
        group_b="receptor_matched",
    )

    assert set(migration_summary["n"]) == {3}
    assert len(paired_table) == 3
    assert "difference_b_minus_a" in paired_table.columns
    assert len(paired_axis.lines) >= 3
    plt.close(migration_axis.figure)
    plt.close(paired_axis.figure)


def test_log_dose_plot_rejects_nonpositive_doses() -> None:
    data = pd.DataFrame(
        {
            "biological_replicate": ["D1", "D2"],
            "dose": [0.0, 1.0],
            "response": [1.0, 2.0],
        }
    )

    with pytest.raises(ValueError, match="positive"):
        plot_dose_response(data, group_col=None, log_x=True)


def test_inference_requires_replicated_biological_units() -> None:
    data = pd.DataFrame(
        {
            "biological_replicate": ["D1", "D1"],
            "product": ["A", "B"],
            "response": [1.0, 2.0],
        }
    )

    with pytest.raises(ValueError, match="two complete biological pairs"):
        compare_groups(
            data,
            group_col="product",
            group_a="A",
            group_b="B",
            paired=True,
            bootstrap_iterations=100,
        )
