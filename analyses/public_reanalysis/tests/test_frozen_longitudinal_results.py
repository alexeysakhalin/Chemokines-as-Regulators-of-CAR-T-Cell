from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from longitudinal_extension.verification import verify_working_tree_manifest

MODULE_ROOT = Path(__file__).resolve().parents[1]
FROZEN = MODULE_ROOT / "results/longitudinal_extension/frozen"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def test_frozen_manifest_digests() -> None:
    for line in (FROZEN / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        assert _sha256(FROZEN / relative) == expected


def test_frozen_provenance_matches_current_working_tree() -> None:
    manifest = json.loads((FROZEN / "analysis_manifest.json").read_text(encoding="utf-8"))
    verify_working_tree_manifest(manifest, MODULE_ROOT)
    assert len(manifest["canonical_outputs"]) == 9
    assert all(not Path(path).is_absolute() for path in manifest["code"])
    assert all(not Path(path).is_absolute() for path in manifest["inputs"])


def test_frozen_primary_results() -> None:
    results = pd.read_csv(FROZEN / "contrast_results.tsv", sep="\t")
    expected = {
        ("GSE197268", "D7-CART_minus_IP"): (
            9,
            7,
            2,
            0.23626373626373628,
            0.1796875,
            0.5390625,
        ),
        ("GSE235760", "Peak_minus_IP"): (
            5,
            4,
            1,
            0.0136159449741511,
            0.375,
            0.75,
        ),
        ("GSE162975", "T1_minus_T0"): (
            10,
            4,
            6,
            -0.00926880659763755,
            0.75390625,
            0.75390625,
        ),
    }
    for key, values in expected.items():
        row = results.set_index(["dataset", "contrast"]).loc[key]
        observed = (
            int(row["n_pairs"]),
            int(row["n_increase"]),
            int(row["n_decrease"]),
            float(row["median_change_fraction"]),
            float(row["sign_test_p"]),
            float(row["holm_p"]),
        )
        assert observed[:3] == values[:3]
        np.testing.assert_allclose(observed[3:], values[3:], rtol=0, atol=1e-15)


def test_frozen_supportive_results() -> None:
    results = pd.read_csv(FROZEN / "contrast_results.tsv", sep="\t")
    indexed = results.set_index(["dataset", "contrast"])
    post_peak = indexed.loc[("GSE162975", "T3-or-T4.5_minus_T1")]
    assert (post_peak["n_pairs"], post_peak["n_decrease"], post_peak["sign_test_p"]) == (
        7,
        6,
        0.125,
    )
    extended = indexed.loc[("GSE162975", "T6-plus_minus_T1")]
    assert (extended["n_pairs"], extended["n_decrease"], extended["sign_test_p"]) == (
        5,
        3,
        1.0,
    )
    gse273 = indexed.loc[("GSE273170", "D14_minus_D7")]
    assert (gse273["n_pairs"], gse273["n_decrease"], gse273["sign_test_p"]) == (
        4,
        3,
        0.625,
    )


def test_frozen_report_rejects_universal_increase() -> None:
    report = (FROZEN / "REPORT.md").read_text(encoding="utf-8")
    assert "none of the three early tests had Holm-adjusted P<0.05" in report
    assert "do not support a universal early CXCR6 increase" in report
    assert "nine of 21 paired patients" in report


def test_gse197268_frozen_screen_and_product_strata() -> None:
    screening = pd.read_csv(FROZEN / "audit/gse197268_screening_exclusions.tsv", sep="\t")
    assert len(screening) == 21
    assert screening["included_at_min25"].sum() == 9
    assert set(screening.loc[screening["included_at_min25"].eq(0), "reason"]) == {
        "below_minimum_eligible_cells"
    }
    strata = pd.read_csv(
        FROZEN / "audit/gse197268_product_stratified_results.tsv", sep="\t"
    ).set_index("subgroup")
    assert (strata.loc["Axi-cel", "n_pairs"], strata.loc["Axi-cel", "n_increase"]) == (
        5,
        5,
    )
    assert (strata.loc["Tisa-cel", "n_pairs"], strata.loc["Tisa-cel", "n_increase"]) == (
        4,
        2,
    )
    np.testing.assert_allclose(
        strata.loc[["Axi-cel", "Tisa-cel"], "median_change_fraction"],
        [0.26120981387478853, -0.007292857971591002],
        rtol=0,
        atol=1e-15,
    )
    patients = pd.read_csv(FROZEN / "audit/gse197268_verified_patient_aggregates.tsv", sep="\t")
    changes = patients["paired_fraction_change"]
    primary = (
        pd.read_csv(FROZEN / "contrast_results.tsv", sep="\t")
        .query("dataset == 'GSE197268'")
        .iloc[0]
    )
    np.testing.assert_allclose(
        [changes.quantile(0.25), changes.median(), changes.quantile(0.75)],
        [
            primary["q1_change_fraction"],
            primary["median_change_fraction"],
            primary["q3_change_fraction"],
        ],
        rtol=0,
        atol=1e-15,
    )
