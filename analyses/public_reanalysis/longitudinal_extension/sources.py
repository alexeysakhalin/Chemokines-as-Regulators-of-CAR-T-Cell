"""Checksum-verified retrieval helpers for public longitudinal inputs."""

from __future__ import annotations

import csv
import hashlib
import shutil
import time
import urllib.request
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_source_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"dataset", "url", "local_path", "bytes", "sha256", "role"}
    if not rows or required.difference(rows[0]):
        raise ValueError(f"Source manifest must contain {sorted(required)}")
    for row in rows:
        if not row["url"].startswith("https://"):
            raise ValueError("Only HTTPS sources are accepted")
        if len(row["sha256"]) != 64 or set(row["sha256"]).difference("0123456789abcdef"):
            raise ValueError(f"Invalid SHA-256 for {row['local_path']}")
        if int(row["bytes"]) <= 0:
            raise ValueError(f"Invalid byte size for {row['local_path']}")
    return rows


def verify_source(path: Path, expected_bytes: int, expected_sha256: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size != expected_bytes:
        raise ValueError(f"Size mismatch for {path}: {path.stat().st_size} != {expected_bytes}")
    observed = digest(path)
    if observed != expected_sha256:
        raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected_sha256}")


def download_source(
    url: str, target: Path, expected_bytes: int, expected_sha256: str, retries: int = 4
) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        verify_source(target, expected_bytes, expected_sha256)
        return
    temporary = target.with_suffix(target.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "cart-public-reanalysis/1.0"})
    error: Exception | None = None
    for attempt in range(retries):
        try:
            with (
                urllib.request.urlopen(request, timeout=300) as response,
                temporary.open("wb") as output,
            ):
                shutil.copyfileobj(response, output, length=1024 * 1024)
            verify_source(temporary, expected_bytes, expected_sha256)
            temporary.replace(target)
            return
        except Exception as caught:  # network failures differ by platform
            error = caught
            temporary.unlink(missing_ok=True)
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"Unable to retrieve verified source {url}") from error


def fetch_manifest(path: Path, root: Path, verify_only: bool = False) -> None:
    for row in load_source_manifest(path):
        target = root / row["local_path"]
        if verify_only:
            verify_source(target, int(row["bytes"]), row["sha256"])
        else:
            download_source(row["url"], target, int(row["bytes"]), row["sha256"])
