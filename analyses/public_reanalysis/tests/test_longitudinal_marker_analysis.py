from __future__ import annotations

from pathlib import Path

import pandas as pd
from longitudinal_extension.analysis import (
    add_gse162975_stage_aliases,
    analyze_contrasts,
    analyze_marker_contrasts,
    summarize_markers,
    summarize_samples,
)
from longitudinal_extension.constants import MARKER_SETS

CONFIG = Path(__file__).resolve().parents[1] / "config"
MEMORY_COLUMNS = tuple(f"gene_{gene}_umi" for gene in MARKER_SETS["memory_associated"])


def _cells(
    dataset: str,
    patient_id: str,
    stage: str,
    n_cells: int,
    *,
    marker_detected: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index in range(n_cells):
        row: dict[str, object] = {
            "dataset": dataset,
            "patient_id": patient_id,
            "product": "test-product",
            "stage": stage,
            "biological_sample_id": f"{patient_id}_{stage}",
            "cell_id": f"{stage}_cell_{index}",
            "total_umi": 100,
            "cxcr6_umi": index % 2,
        }
        row.update({column: int(marker_detected) for column in MEMORY_COLUMNS})
        rows.append(row)
    return rows


def _contrast(dataset: str, baseline: str, followup: str, *, min_cells: int) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "dataset": dataset,
                "contrast": f"{followup}_minus_{baseline}",
                "baseline_stage": baseline,
                "followup_stage": followup,
                "min_cells": min_cells,
                "role": "test",
                "family_id": "none",
            }
        ]
    )


def test_marker_contrast_excludes_ot02_below_early_cell_threshold() -> None:
    rows = [
        *_cells("GSE162975", "OT02", "T0", 7, marker_detected=False),
        *_cells("GSE162975", "OT02", "T1", 30, marker_detected=True),
        *_cells("GSE162975", "CR07", "T0", 25, marker_detected=False),
        *_cells("GSE162975", "CR07", "T1", 25, marker_detected=True),
    ]
    cells = pd.DataFrame(rows)
    inclusion, _ = analyze_contrasts(
        summarize_samples(cells), _contrast("GSE162975", "T0", "T1", min_cells=25)
    )

    ot02 = inclusion.loc[inclusion["patient_id"].eq("OT02")].iloc[0]
    assert not ot02["included"]
    assert ot02["reason"] == "below_minimum_eligible_cells"
    marker_results = analyze_marker_contrasts(summarize_markers(cells), inclusion)
    assert marker_results.iloc[0]["n_pairs"] == 1


def test_gse273170_marker_contrast_uses_only_four_locked_pairs() -> None:
    eligible_counts = {
        "CAR116": (1, 6),
        "CAR118": (2, 17),
        "CAR122": (7, 0),
        "CAR128": (83, 26),
        "CAR137": (39, 102),
        "CAR138": (11, 1),
        "CAR139": (76, 14),
        "CAR141": (3, 19),
        "CAR142": (122, 20),
        "CAR147": (3, 19),
    }
    rows: list[dict[str, object]] = []
    for patient_id, (d7_cells, d14_cells) in eligible_counts.items():
        rows.extend(_cells("GSE273170", patient_id, "D7", d7_cells, marker_detected=False))
        rows.extend(_cells("GSE273170", patient_id, "D14", d14_cells, marker_detected=True))
    cells = pd.DataFrame(rows)
    inclusion, _ = analyze_contrasts(
        summarize_samples(cells), _contrast("GSE273170", "D7", "D14", min_cells=10)
    )

    assert inclusion["included"].sum() == 4
    marker_results = analyze_marker_contrasts(summarize_markers(cells), inclusion)
    assert marker_results.iloc[0]["n_pairs"] == 4


def test_configured_late_aliases_are_added_to_samples_and_markers() -> None:
    aliases = pd.read_csv(
        CONFIG / "gse162975_crosswalk.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    extended = aliases.loc[aliases["canonical_stage"].eq("T6-plus")]
    assert extended["deposited_sampling_stage"].tolist() == ["T6", "T9", "T12", "T15", "T16"]
    rows = [
        *_cells("GSE162975", "P1", "T1", 25, marker_detected=False),
        *_cells("GSE162975", "P1", "T15", 25, marker_detected=True),
        *_cells("GSE162975", "P1", "T16", 25, marker_detected=False),
    ]
    cells = pd.DataFrame(rows)
    samples = add_gse162975_stage_aliases(summarize_samples(cells), aliases)
    markers = add_gse162975_stage_aliases(summarize_markers(cells), aliases)

    late_sample = samples.loc[samples["stage"].eq("T6-plus")]
    late_marker = markers.loc[markers["stage"].eq("T6-plus")]
    assert len(late_sample) == 1
    assert len(late_marker) == 1
    assert late_marker.iloc[0]["mean_marker_detection_fraction"] == 1.0
    inclusion, _ = analyze_contrasts(samples, _contrast("GSE162975", "T1", "T6-plus", min_cells=25))
    marker_results = analyze_marker_contrasts(markers, inclusion)
    assert marker_results.iloc[0]["n_pairs"] == 1
    assert marker_results.iloc[0]["median_change_fraction"] == 1.0
