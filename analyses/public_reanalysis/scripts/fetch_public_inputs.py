#!/usr/bin/env python3
"""Download and verify the exact public inputs used by the analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import tarfile
import time
import urllib.request
from pathlib import Path


SPATIAL_MEMBERS = (
    "GSM8968967_barcodes.tsv.gz",
    "GSM8968967_features.tsv.gz",
    "GSM8968967_matrix.mtx.gz",
    "GSM8968967_scalefactors_json.json.gz",
    "GSM8968967_tissue_hires_image.png.gz",
    "GSM8968967_tissue_positions.parquet.gz",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"dataset", "url", "local_path", "bytes", "sha256"}
    if not rows or required.difference(rows[0]):
        raise ValueError(f"Manifest must contain {sorted(required)}")
    return rows


def verify(path: Path, expected_bytes: int, expected_sha256: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed_bytes = path.stat().st_size
    if observed_bytes != expected_bytes:
        raise ValueError(f"Size mismatch for {path}: {observed_bytes} != {expected_bytes}")
    observed_sha256 = digest(path)
    if observed_sha256 != expected_sha256:
        raise ValueError(f"SHA-256 mismatch for {path}: {observed_sha256} != {expected_sha256}")


def download(url: str, target: Path, expected_bytes: int, expected_sha256: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        verify(target, expected_bytes, expected_sha256)
        return
    temporary = target.with_suffix(target.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "public-data-reproduction/1.0"})
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            verify(temporary, expected_bytes, expected_sha256)
            temporary.replace(target)
            return
        except Exception as error:  # network errors vary by platform
            last_error = error
            if temporary.exists():
                temporary.unlink()
            if attempt < 3:
                time.sleep(2**attempt)
    raise RuntimeError(f"Unable to download {url}") from last_error


def extract_spatial_archive(archive: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    wanted = set(SPATIAL_MEMBERS)
    with tarfile.open(archive, "r") as handle:
        members = {member.name: member for member in handle.getmembers()}
        missing = wanted.difference(members)
        if missing:
            raise ValueError(f"Archive is missing required members: {sorted(missing)}")
        for name in SPATIAL_MEMBERS:
            member = members[name]
            if not member.isfile() or Path(name).name != name:
                raise ValueError(f"Unsafe archive member: {name}")
            source = handle.extractfile(member)
            if source is None:
                raise ValueError(f"Cannot read archive member: {name}")
            target = output_dir / name
            temporary = target.with_suffix(target.suffix + ".part")
            with temporary.open("wb") as destination:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
            temporary.replace(target)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("config/sources.tsv"))
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--list", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = load_manifest(args.root / args.manifest)
    for row in rows:
        target = args.root / row["local_path"]
        if args.list:
            print(f"{row['dataset']}\t{target}\t{row['bytes']}\t{row['sha256']}")
            continue
        if args.verify_only:
            verify(target, int(row["bytes"]), row["sha256"])
        else:
            download(row["url"], target, int(row["bytes"]), row["sha256"])
    if not args.list and not args.verify_only:
        archive = args.root / "data/raw/GSE269379/GSE269379_RAW.tar"
        extract_spatial_archive(archive, args.root / "data/raw/GSE269379/extracted")


if __name__ == "__main__":
    main()
