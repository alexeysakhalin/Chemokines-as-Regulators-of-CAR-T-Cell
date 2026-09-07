from __future__ import annotations

from pathlib import Path

import pytest
from longitudinal_extension.sources import load_source_manifest, verify_source


def test_source_manifest_requires_sha256(tmp_path: Path) -> None:
    path = tmp_path / "sources.tsv"
    path.write_text(
        "dataset\turl\tlocal_path\tbytes\tsha256\trole\n"
        "TEST\thttps://example.org/a\ta\t1\tpending\tdata\n"
    )
    with pytest.raises(ValueError, match="Invalid SHA-256"):
        load_source_manifest(path)


def test_source_verification_detects_byte_mismatch(tmp_path: Path) -> None:
    path = tmp_path / "input"
    path.write_bytes(b"abc")
    with pytest.raises(ValueError, match="Size mismatch"):
        verify_source(path, 4, "0" * 64)


def test_source_verification_detects_digest_mismatch(tmp_path: Path) -> None:
    path = tmp_path / "input"
    path.write_bytes(b"abc")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_source(path, 3, "0" * 64)
