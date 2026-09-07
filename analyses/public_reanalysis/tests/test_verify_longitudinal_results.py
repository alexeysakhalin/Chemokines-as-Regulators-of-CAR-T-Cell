from __future__ import annotations

import json
from pathlib import Path

import pytest
from longitudinal_extension.common import sha256
from longitudinal_extension.verification import (
    verify_checksums,
    verify_observed_manifest,
    verify_working_tree_manifest,
)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")


def _expected_manifest(module_root: Path, observed: Path) -> dict[str, object]:
    code = module_root / "longitudinal_extension" / "workflow.py"
    source = module_root / "config" / "sources.tsv"
    audit = module_root / "results" / "frozen" / "audit.tsv"
    for path, content in (
        (code, "locked code\n"),
        (source, "locked sources\n"),
        (audit, "patient\tchange\nP1\t0.1\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    canonical_outputs = {
        name: sha256(observed / name) for name in ("contrast_results.tsv", "REPORT.md")
    }
    return {
        "analysis": "test longitudinal analysis",
        "seed": 1,
        "python_requirement": ">=3.11,<3.12",
        "inputs": {"source/sha256/" + "a" * 64: "a" * 64},
        "code": {"longitudinal_extension/workflow.py": sha256(code)},
        "source_manifest": {"path": "config/sources.tsv", "sha256": sha256(source)},
        "audited_patient_rows": {"results/frozen/audit.tsv": sha256(audit)},
        "canonical_outputs": canonical_outputs,
        "independent_unit": "patient",
        "cell_level_p_values": False,
        "cross_study_pooling": False,
    }


def test_verifier_checks_working_tree_and_complete_outputs(tmp_path: Path) -> None:
    module_root = tmp_path / "module"
    observed = tmp_path / "observed"
    observed.mkdir()
    (observed / "contrast_results.tsv").write_text(
        "dataset\tcontrast\trole\nTEST\tD7_minus_IP\tprincipal\n", encoding="utf-8"
    )
    (observed / "REPORT.md").write_text("# Locked report\n", encoding="utf-8")
    expected = _expected_manifest(module_root, observed)
    _write_json(observed / "analysis_manifest.json", expected)

    verify_working_tree_manifest(expected, module_root)
    verify_observed_manifest(observed, expected)

    (observed / "contrast_results.tsv").write_text(
        "dataset\tcontrast\trole\nTEST\tD7_minus_IP\tsupportive\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="canonical_outputs checksum mismatch"):
        verify_observed_manifest(observed, expected)


def test_verifier_rejects_stale_code_digest(tmp_path: Path) -> None:
    module_root = tmp_path / "module"
    observed = tmp_path / "observed"
    observed.mkdir()
    (observed / "contrast_results.tsv").write_text("locked\n", encoding="utf-8")
    (observed / "REPORT.md").write_text("locked\n", encoding="utf-8")
    expected = _expected_manifest(module_root, observed)
    (module_root / "longitudinal_extension" / "workflow.py").write_text(
        "changed code\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="code checksum mismatch"):
        verify_working_tree_manifest(expected, module_root)


def test_verifier_rejects_absolute_manifest_paths(tmp_path: Path) -> None:
    source = tmp_path / "sources.tsv"
    source.write_text("locked\n", encoding="utf-8")
    manifest = {
        "code": {str(source): sha256(source)},
        "source_manifest": {"path": "sources.tsv", "sha256": sha256(source)},
    }
    with pytest.raises(ValueError, match="Absolute path is forbidden"):
        verify_working_tree_manifest(manifest, tmp_path)


def test_checksum_inventory_must_cover_every_frozen_file(tmp_path: Path) -> None:
    frozen = tmp_path / "frozen"
    frozen.mkdir()
    result = frozen / "contrast_results.tsv"
    result.write_text("locked\n", encoding="utf-8")
    (frozen / "unlisted.tsv").write_text("not in inventory\n", encoding="utf-8")
    (frozen / "SHA256SUMS").write_text(
        f"{sha256(result)}  contrast_results.tsv\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Incomplete SHA256SUMS inventory"):
        verify_checksums(frozen)
