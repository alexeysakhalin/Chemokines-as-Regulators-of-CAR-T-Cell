from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_public_reanalysis import (  # noqa: E402
    CHEMOKINE_RECEPTORS,
    STATE_MARKER_SETS,
    _matched_null,
    analyze_gse125881,
    benjamini_hochberg,
)


def test_benjamini_hochberg_preserves_index_and_monotonicity() -> None:
    values = pd.Series([0.01, 0.04, 0.03], index=["a", "b", "c"])
    observed = benjamini_hochberg(values)
    assert list(observed.index) == ["a", "b", "c"]
    np.testing.assert_allclose(observed.to_numpy(), [0.03, 0.04, 0.04])


def test_gse125881_summary_uses_patient_phase_units(tmp_path: Path) -> None:
    metadata = pd.DataFrame(
        {
            "Cell": ["c1", "c2", "c3", "c4"],
            "Patient": ["P1", "P1", "P1", "P1"],
            "Group": ["IP", "Early", "Late", "Very late"],
            "nUMI": [10, 20, 30, 40],
        }
    )
    metadata_path = tmp_path / "metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    matrix_path = tmp_path / "matrix.csv.gz"
    with gzip.open(matrix_path, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "c1", "c2", "c3", "c4"])
        values = {gene: [0, 0, 0, 0] for gene in CHEMOKINE_RECEPTORS}
        values["CXCR6"] = [0, 2, 1, 0]
        for genes in STATE_MARKER_SETS.values():
            for gene in genes:
                values.setdefault(gene, [0, 0, 0, 0])
        values["CCR7"] = [1, 1, 1, 1]
        for gene in sorted(values):
            writer.writerow([gene, *values[gene]])
    result = analyze_gse125881(matrix_path, metadata_path)
    summary = result["summary"].set_index("Group")
    assert summary.loc["IP", "cxcr6_detected_fraction"] == 0
    assert summary.loc["Early", "cxcr6_detected_fraction"] == 1
    assert summary.loc["Early", "cxcr6_counts_per_million_umis"] == 100000


def test_depth_matched_randomization_is_seeded() -> None:
    distances = np.linspace(0, 99, 100)
    observed_mask = np.zeros(100, dtype=bool)
    observed_mask[[1, 11, 21, 31, 41, 51, 61, 71, 81, 91]] = True
    strata = np.repeat(np.arange(10), 10)
    first = _matched_null(distances, observed_mask, strata, 50, 99, 17)
    second = _matched_null(distances, observed_mask, strata, 50, 99, 17)
    assert first[0] == second[0]
    np.testing.assert_array_equal(first[1], second[1])
    assert first[2] == second[2]


def test_frozen_results_match_declared_invariants() -> None:
    frozen = ROOT / "results" / "frozen" / "tables"
    donor = pd.read_csv(frozen / "scRNA_donor_summary.tsv", sep="\t")
    onset = donor[donor["analysis_group"].eq("ICANS onset")]
    comparator = donor[donor["analysis_group"].eq("IIH comparator")]
    assert onset["n_t_cells"].sum() == 4201
    assert comparator["n_t_cells"].sum() == 3105
    primary = pd.read_csv(frozen / "scRNA_primary_tests.tsv", sep="\t")
    assert primary.iloc[0]["p_value_exact"] == 0.015873016
    longitudinal = pd.read_csv(frozen / "gse125881_summary.tsv", sep="\t")
    assert longitudinal["n_cells"].sum() == 62167
    assert longitudinal.groupby("Patient").size().eq(4).all()
    spatial = pd.read_csv(frozen / "spatial_marker_summary.tsv", sep="\t").set_index("marker")
    assert spatial.loc["CXCL16", "positive_bins"] == 960
    assert spatial.loc["CXCR6", "positive_bins"] == 14


def test_gse125881_crosswalk_matches_frozen_phase_counts() -> None:
    crosswalk = pd.read_csv(ROOT / "config" / "gse125881_sample_crosswalk.tsv", sep="\t")
    frozen = pd.read_csv(ROOT / "results" / "frozen" / "tables" / "gse125881_summary.tsv", sep="\t")
    expected = (
        crosswalk.groupby(["patient", "authors_metadata_group"], as_index=False)["n_cells"]
        .sum()
        .rename(columns={"patient": "Patient", "authors_metadata_group": "Group"})
        .sort_values(["Patient", "Group"])
        .reset_index(drop=True)
    )
    observed = (
        frozen[["Patient", "Group", "n_cells"]]
        .sort_values(["Patient", "Group"])
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(observed, expected)
