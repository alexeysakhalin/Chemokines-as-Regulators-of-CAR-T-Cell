from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest

from scripts.fetch_public_data import PublicFile, read_manifest, resolve_destination, verify


def test_verify_checks_size_and_digest(tmp_path: Path) -> None:
    payload = b"counts\n"
    path = tmp_path / "input.tsv"
    path.write_bytes(payload)
    record = PublicFile(
        dataset_id="GSETEST",
        file_id="input.tsv",
        url="https://example.invalid/input.tsv",
        size_bytes=len(payload),
        sha256=hashlib.sha256(payload).hexdigest(),
        access_class="public",
        local_relative_path=Path("data/raw/GSETEST/input.tsv"),
    )
    assert verify(path, record) == (True, "verified")
    path.write_bytes(payload + b"x")
    valid, message = verify(path, record)
    assert not valid
    assert "size mismatch" in message


def test_manifest_rejects_unsafe_path(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    fields = [
        "dataset_id",
        "file_id",
        "exact_file_url",
        "size_bytes",
        "sha256",
        "access_class",
        "local_relative_path",
    ]
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "dataset_id": "GSETEST",
                "file_id": "input.tsv",
                "exact_file_url": "https://example.invalid/input.tsv",
                "size_bytes": "1",
                "sha256": "0" * 64,
                "access_class": "public",
                "local_relative_path": "../input.tsv",
            }
        )
    with pytest.raises(ValueError, match="Unsafe local path"):
        read_manifest(manifest)


def test_destination_must_remain_under_raw(tmp_path: Path) -> None:
    record = PublicFile(
        dataset_id="GSETEST",
        file_id="input.tsv",
        url="https://example.invalid/input.tsv",
        size_bytes=1,
        sha256="0" * 64,
        access_class="public",
        local_relative_path=Path("results/input.tsv"),
    )
    with pytest.raises(ValueError, match="data/raw"):
        resolve_destination(tmp_path, record)
