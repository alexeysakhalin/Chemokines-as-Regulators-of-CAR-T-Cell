"""Verification of frozen provenance and canonical longitudinal outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

KEY = ["dataset", "contrast"]
NUMERIC = [
    "n_pairs",
    "n_increase",
    "n_decrease",
    "n_unchanged",
    "median_change_fraction",
    "sign_test_p",
    "holm_p",
]
MANIFEST_PARITY_FIELDS = (
    "analysis",
    "seed",
    "python_requirement",
    "inputs",
    "code",
    "source_manifest",
    "independent_unit",
    "cell_level_p_values",
    "cross_study_pooling",
)


def sha256(path: Path) -> str:
    """Return the streaming SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict[str, Any]:
    """Read a JSON-object manifest."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Manifest must contain a JSON object: {path}")
    return payload


def _validated_digest(value: object, *, label: str) -> str:
    digest = str(value)
    if len(digest) != 64 or set(digest).difference("0123456789abcdef"):
        raise ValueError(f"Invalid SHA-256 digest for {label}: {digest}")
    return digest


def _resolve_declared_path(root: Path, relative: object, *, label: str) -> Path:
    candidate = Path(str(relative))
    if candidate.is_absolute():
        raise ValueError(f"Absolute path is forbidden in {label}: {candidate}")
    resolved_root = root.resolve()
    resolved = (resolved_root / candidate).resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as error:
        raise ValueError(f"Path escapes {label} root: {candidate}") from error
    return resolved


def _verify_digest_mapping(mapping: object, root: Path, *, label: str) -> None:
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{label} must be a non-empty path-to-SHA-256 mapping")
    for relative, expected_value in sorted(mapping.items()):
        expected = _validated_digest(expected_value, label=f"{label}:{relative}")
        path = _resolve_declared_path(root, relative, label=label)
        if not path.is_file():
            raise ValueError(f"Declared {label} file is missing: {relative}")
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"{label} checksum mismatch: {relative}: {observed} != {expected}")


def verify_checksums(frozen: Path) -> None:
    """Verify every frozen artifact and require complete checksum coverage."""
    checksum_path = frozen / "SHA256SUMS"
    declared: dict[str, str] = {}
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        try:
            expected_value, relative = line.split("  ", 1)
        except ValueError as error:
            raise ValueError(f"Malformed SHA256SUMS line {line_number}") from error
        if relative in declared:
            raise ValueError(f"Duplicate SHA256SUMS entry: {relative}")
        expected = _validated_digest(expected_value, label=f"SHA256SUMS:{relative}")
        path = _resolve_declared_path(frozen, relative, label="SHA256SUMS")
        if not path.is_file():
            raise ValueError(f"Frozen file listed but missing: {relative}")
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"Frozen-file checksum mismatch: {relative}: {observed} != {expected}")
        declared[relative] = expected
    actual = {
        path.relative_to(frozen).as_posix()
        for path in frozen.rglob("*")
        if path.is_file() and path != checksum_path
    }
    if set(declared) != actual:
        missing = sorted(actual.difference(declared))
        stale = sorted(set(declared).difference(actual))
        raise ValueError(f"Incomplete SHA256SUMS inventory; missing={missing}, stale={stale}")


def verify_working_tree_manifest(manifest: dict[str, Any], module_root: Path) -> None:
    """Reject a freeze whose declared code, source manifest, or audit digests are stale."""
    _verify_digest_mapping(manifest.get("code"), module_root, label="code")
    source_manifest = manifest.get("source_manifest")
    if not isinstance(source_manifest, dict) or set(source_manifest) != {"path", "sha256"}:
        raise ValueError("source_manifest must contain exactly path and sha256")
    _verify_digest_mapping(
        {source_manifest["path"]: source_manifest["sha256"]},
        module_root,
        label="source_manifest",
    )
    audited = manifest.get("audited_patient_rows")
    if audited is not None:
        _verify_digest_mapping(audited, module_root, label="audited_patient_rows")


def verify_results(observed: Path, frozen: Path) -> None:
    """Provide a focused diagnostic for the locked cohort-level numeric endpoints."""
    expected = pd.read_csv(frozen / "contrast_results.tsv", sep="\t").set_index(KEY)
    actual = pd.read_csv(observed / "contrast_results.tsv", sep="\t").set_index(KEY)
    missing = expected.index.difference(actual.index)
    if len(missing):
        raise ValueError(f"Observed run is missing contrasts: {list(missing)}")
    for index, expected_row in expected.iterrows():
        actual_row = actual.loc[index]
        for column in NUMERIC:
            expected_value = float(expected_row[column])
            actual_value = float(actual_row[column])
            if np.isnan(expected_value) and np.isnan(actual_value):
                continue
            if not np.isclose(actual_value, expected_value, rtol=0, atol=1e-12):
                raise ValueError(
                    f"Result mismatch for {index}, {column}: {actual_value} != {expected_value}"
                )


def verify_observed_manifest(observed: Path, expected_manifest: dict[str, Any]) -> None:
    """Check provenance parity and every byte of each declared canonical output."""
    actual_manifest = load_manifest(observed / "analysis_manifest.json")
    for field in MANIFEST_PARITY_FIELDS:
        if field not in expected_manifest:
            raise ValueError(f"Frozen manifest does not declare required field: {field}")
        if actual_manifest.get(field) != expected_manifest[field]:
            raise ValueError(f"Observed manifest differs from the freeze for field: {field}")

    expected_outputs = expected_manifest.get("canonical_outputs")
    actual_outputs = actual_manifest.get("canonical_outputs")
    if not isinstance(expected_outputs, dict) or not expected_outputs:
        raise ValueError("Frozen manifest must declare canonical_outputs")
    if actual_outputs != expected_outputs:
        raise ValueError("Observed canonical-output digest map differs from the freeze")
    _verify_digest_mapping(expected_outputs, observed, label="canonical_outputs")
