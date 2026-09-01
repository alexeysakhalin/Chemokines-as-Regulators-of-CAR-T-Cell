from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from chemokine_cart.single_cell import (
    aggregate_pseudobulk,
    bootstrap_paired_change_ci,
    bootstrap_patient_ci,
    calculate_state_fractions,
    plot_state_fractions,
    split_car_endogenous,
    validate_required_columns,
)


def test_validate_required_columns_reports_missing_schema() -> None:
    table = pd.DataFrame({"patient_id": ["P01"]})

    with pytest.raises(ValueError, match="timepoint, cell_state"):
        validate_required_columns(
            table,
            ("patient_id", "timepoint", "cell_state"),
            table_name="cells",
        )


def test_split_car_endogenous_uses_only_documented_binary_gate() -> None:
    cells = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "car_positive": pd.Series([1, 0, 1, 0], dtype="int8"),
        }
    )

    car_cells, endogenous_cells = split_car_endogenous(cells)

    assert car_cells["cell_id"].tolist() == ["c1", "c3"]
    assert endogenous_cells["cell_id"].tolist() == ["c2", "c4"]
    car_cells.loc[:, "cell_id"] = "changed"
    assert cells.loc[0, "cell_id"] == "c1"

    invalid = cells.assign(car_positive=["CAR+", "CAR-", "CAR+", "CAR-"])
    with pytest.raises(TypeError, match="document the CAR-detection gate upstream"):
        split_car_endogenous(invalid)


def test_aggregate_pseudobulk_sums_counts_and_filters_low_cell_groups() -> None:
    cells = pd.DataFrame(
        {
            "patient_id": ["P01", "P01", "P01", "P02"],
            "timepoint_id": ["D0", "D0", "D7", "D0"],
            "t_cell_origin": ["CAR+", "CAR+", "CAR+", "endogenous"],
            "CXCR3": [1, 2, 4, 8],
            "CXCL16_ADT": [3, 5, 7, 9],
        }
    )

    result = aggregate_pseudobulk(
        cells,
        ("CXCR3", "CXCL16_ADT"),
        min_cells=2,
    )

    assert result.to_dict("records") == [
        {
            "patient_id": "P01",
            "timepoint_id": "D0",
            "t_cell_origin": "CAR+",
            "n_cells": 2,
            "CXCR3": 3,
            "CXCL16_ADT": 8,
        }
    ]


def test_aggregate_pseudobulk_rejects_normalised_negative_values() -> None:
    cells = pd.DataFrame(
        {
            "patient_id": ["P01"],
            "timepoint_id": ["D0"],
            "t_cell_origin": ["CAR+"],
            "gene": [-0.25],
        }
    )

    with pytest.raises(ValueError, match="must be non-negative"):
        aggregate_pseudobulk(cells, ("gene",))

    cells.loc[:, "gene"] = 0.25
    with pytest.raises(ValueError, match="untransformed integer-valued counts"):
        aggregate_pseudobulk(cells, ("gene",))


def test_calculate_state_fractions_is_patient_level_and_completes_zero_states() -> None:
    cells = pd.DataFrame(
        {
            "patient_id": ["P01", "P01", "P01", "P02", "P02"],
            "timepoint_id": ["D7"] * 5,
            "t_cell_origin": ["CAR+"] * 5,
            "cell_state": [
                "memory",
                "memory",
                "effector",
                "effector",
                "cycling",
            ],
        }
    )

    result = calculate_state_fractions(
        cells,
        states=("memory", "effector", "dysfunctional"),
        denominator="all",
    )

    p01 = result.loc[result["patient_id"] == "P01"].set_index("cell_state")
    p02 = result.loc[result["patient_id"] == "P02"].set_index("cell_state")
    assert p01.loc["memory", "n_cells"] == 2
    assert p01.loc["memory", "fraction"] == pytest.approx(2 / 3)
    assert p01.loc["dysfunctional", "fraction"] == 0
    assert p02.loc["effector", "fraction"] == pytest.approx(1 / 2)
    assert p02.loc["memory", "fraction"] == 0
    assert (p02["n_total"] == 2).all()


def test_state_fractions_report_other_and_low_cell_qc() -> None:
    cells = pd.DataFrame(
        {
            "patient_id": ["P01", "P01", "P01"],
            "timepoint_id": ["D7"] * 3,
            "t_cell_origin": ["CAR_T"] * 3,
            "cell_state": ["memory", "cycling", "interferon_responsive"],
        }
    )

    result = calculate_state_fractions(
        cells,
        states=("memory", "effector", "dysfunctional"),
        other_label="other",
        min_cells_per_group=4,
    ).set_index("cell_state")

    assert result.loc["other", "n_cells"] == 2
    assert result["fraction"].sum() == pytest.approx(1.0)
    assert not result["passes_min_cells"].any()


def test_bootstrap_patient_ci_is_reproducible_and_rejects_pseudoreplication() -> None:
    patient_metrics = pd.DataFrame(
        {
            "patient_id": ["P01", "P02", "P03", "P01", "P02", "P03"],
            "cell_state": ["memory"] * 3 + ["effector"] * 3,
            "fraction": [0.2, 0.4, 0.6, 0.7, 0.5, 0.3],
        }
    )

    first = bootstrap_patient_ci(
        patient_metrics,
        "fraction",
        group_columns=("cell_state",),
        n_boot=500,
        random_state=19,
    )
    second = bootstrap_patient_ci(
        patient_metrics.sample(frac=1, random_state=4),
        "fraction",
        group_columns=("cell_state",),
        n_boot=500,
        random_state=19,
    )

    pd.testing.assert_frame_equal(first, second)
    assert first.set_index("cell_state").loc["memory", "estimate"] == pytest.approx(0.4)
    assert (first["n_patients"] == 3).all()

    duplicated = pd.concat([patient_metrics, patient_metrics.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="one value per patient"):
        bootstrap_patient_ci(
            duplicated,
            "fraction",
            group_columns=("cell_state",),
        )


def test_bootstrap_records_insufficient_strata_without_a_false_interval() -> None:
    result = bootstrap_patient_ci(
        pd.DataFrame({"patient_id": ["P01", "P02"], "fraction": [0.2, 0.8]}),
        "fraction",
        min_patients=3,
        n_boot=100,
        on_insufficient="record",
    ).iloc[0]

    assert result["status"] == "insufficient_patients"
    assert result["estimate"] == pytest.approx(0.5)
    assert np.isnan(result["ci_low"])
    assert np.isnan(result["ci_high"])


def test_paired_bootstrap_uses_complete_within_patient_changes() -> None:
    metrics = pd.DataFrame(
        {
            "patient_id": ["P01", "P01", "P02", "P02", "P03"],
            "timepoint_id": ["PRODUCT", "D7", "PRODUCT", "D7", "PRODUCT"],
            "cell_state": ["memory"] * 5,
            "fraction": [0.2, 0.5, 0.4, 0.5, 0.9],
        }
    )

    result = bootstrap_paired_change_ci(
        metrics,
        "fraction",
        reference_timepoint="PRODUCT",
        group_columns=("cell_state",),
        min_patients=2,
        n_boot=500,
        random_state=17,
    ).iloc[0]

    assert result["n_patients"] == 2
    assert result["estimate"] == pytest.approx(0.2)
    assert result["status"] == "ok"


def test_plot_state_fractions_draws_only_input_observations() -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fractions = pd.DataFrame(
        {
            "patient_id": ["P01", "P01", "P01", "P01", "P01", "P01"],
            "timepoint_id": ["D0", "D7"] * 3,
            "t_cell_origin": ["CAR+"] * 6,
            "cell_state": ["memory"] * 2 + ["effector"] * 2 + ["dysfunctional"] * 2,
            "fraction": [0.7, 0.5, 0.2, 0.3, 0.1, 0.2],
        }
    )

    figure, axes = plot_state_fractions(fractions)

    assert len(axes) == 3
    for axis in axes:
        plotted_offsets = np.concatenate(
            [collection.get_offsets() for collection in axis.collections]
        )
        assert plotted_offsets.shape[0] == 2
    observed_y = sorted(
        float(value)
        for axis in axes
        for collection in axis.collections
        for value in collection.get_offsets()[:, 1]
    )
    assert observed_y == pytest.approx(sorted(fractions["fraction"].tolist()))
    plt.close(figure)
