from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from longitudinal_extension.analysis import (
    analyze_contrasts,
    leave_one_out,
    summarize_samples,
)
from longitudinal_extension.common import sha256
from longitudinal_extension.output import build_report, write_outputs


def test_report_states_patient_unit() -> None:
    results = pd.DataFrame(
        [
            {
                "dataset": "GSE197268",
                "contrast": "D7_minus_IP",
                "family_id": "early_cxcr6_fraction_three_cohort",
                "n_pairs": 9,
                "n_increase": 7,
                "n_decrease": 2,
                "n_unchanged": 0,
                "median_change_percentage_points": 23.626,
                "sign_test_p": 0.1796875,
                "holm_p": 0.5390625,
            },
            {
                "dataset": "GSE235760",
                "contrast": "Peak_minus_IP",
                "family_id": "early_cxcr6_fraction_three_cohort",
                "n_pairs": 5,
                "n_increase": 4,
                "n_decrease": 1,
                "n_unchanged": 0,
                "median_change_percentage_points": 1.3616,
                "sign_test_p": 0.375,
                "holm_p": 0.75,
            },
            {
                "dataset": "GSE162975",
                "contrast": "T1_minus_T0",
                "family_id": "early_cxcr6_fraction_three_cohort",
                "n_pairs": 10,
                "n_increase": 4,
                "n_decrease": 6,
                "n_unchanged": 0,
                "median_change_percentage_points": -0.9269,
                "sign_test_p": 0.75390625,
                "holm_p": 0.75390625,
            },
        ]
    )
    report = build_report(results)
    assert "patient is the independent unit" in report
    assert "none of the three early tests had Holm-adjusted P<0.05" in report
    assert "not a universal early increase" in report


def test_outputs_are_byte_identical(
    tmp_path: Path, synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    samples = summarize_samples(synthetic_cells)
    inclusion, results = analyze_contrasts(samples, synthetic_contrasts)
    empty = pd.DataFrame()
    first = tmp_path / "first"
    second = tmp_path / "second"
    script = Path(__file__)
    kwargs = {
        "samples": samples,
        "markers": empty,
        "inclusion": inclusion,
        "results": results,
        "marker_results": empty,
        "leave_one_out": leave_one_out(inclusion),
        "input_paths": [],
        "code_paths": [script],
        "seed": 123,
    }
    write_outputs(first, **kwargs)
    write_outputs(second, **kwargs)
    assert {path.name: sha256(path) for path in first.iterdir()} == {
        path.name: sha256(path) for path in second.iterdir()
    }


def test_checksum_manifest_excludes_itself(
    tmp_path: Path, synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    samples = summarize_samples(synthetic_cells)
    inclusion, results = analyze_contrasts(samples, synthetic_contrasts)
    output = tmp_path / "run"
    write_outputs(
        output,
        samples=samples,
        markers=pd.DataFrame(),
        inclusion=inclusion,
        results=results,
        marker_results=pd.DataFrame(),
        leave_one_out=pd.DataFrame(),
        input_paths=[],
        code_paths=[],
        seed=1,
    )
    text = (output / "SHA256SUMS").read_text()
    assert "SHA256SUMS" not in text
    assert "REPORT.md" in text


def test_manifest_is_independent_of_checkout_and_raw_cache_paths(
    tmp_path: Path, synthetic_cells: pd.DataFrame, synthetic_contrasts: pd.DataFrame
) -> None:
    samples = summarize_samples(synthetic_cells)
    inclusion, results = analyze_contrasts(samples, synthetic_contrasts)
    empty = pd.DataFrame()
    roots = [tmp_path / "checkout-a", tmp_path / "checkout-b"]
    outputs = []
    for index, root in enumerate(roots):
        config = root / "config" / "longitudinal_sources.tsv"
        code = root / "longitudinal_extension" / "workflow.py"
        raw = (
            root / "data" / "raw" / "TEST" / "matrix.bin"
            if index == 0
            else tmp_path / "renamed-cache" / "download.bin"
        )
        config.parent.mkdir(parents=True)
        code.parent.mkdir(parents=True)
        raw.parent.mkdir(parents=True, exist_ok=True)
        config.write_text("locked source manifest\n", encoding="utf-8")
        code.write_text("SEED = 123\n", encoding="utf-8")
        raw.write_bytes(b"same public matrix")
        output = root / "results" / "local"
        write_outputs(
            output,
            samples=samples,
            markers=empty,
            inclusion=inclusion,
            results=results,
            marker_results=empty,
            leave_one_out=leave_one_out(inclusion),
            input_paths=[config, raw],
            code_paths=[code],
            seed=123,
            manifest_root=root,
            source_manifest_path=config,
        )
        outputs.append(output)

    first = json.loads((outputs[0] / "analysis_manifest.json").read_text(encoding="utf-8"))
    second = json.loads((outputs[1] / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert first == second
    assert first["python_requirement"] == ">=3.11,<3.12"
    assert set(first["inputs"]) == {
        "config/longitudinal_sources.tsv",
        f"source/sha256/{sha256(roots[0] / 'data/raw/TEST/matrix.bin')}",
    }
    assert all(not Path(path).is_absolute() for path in first["code"])
