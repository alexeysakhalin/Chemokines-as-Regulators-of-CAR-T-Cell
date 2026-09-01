#!/usr/bin/env python3
"""Download checksum-pinned public inputs into an ignored local data directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import sys
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CHUNK_SIZE = 1024 * 1024


@dataclass(frozen=True)
class PublicFile:
    dataset_id: str
    file_id: str
    url: str
    size_bytes: int
    sha256: str
    access_class: str
    local_relative_path: Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[PublicFile]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "dataset_id",
        "file_id",
        "exact_file_url",
        "size_bytes",
        "sha256",
        "access_class",
        "local_relative_path",
    }
    if not rows:
        raise ValueError(f"No file records found in {path}")
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(sorted(missing))}")

    records: list[PublicFile] = []
    for row_number, row in enumerate(rows, start=2):
        digest = row["sha256"].strip().lower()
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError(f"Invalid SHA-256 on manifest row {row_number}")
        local_path = Path(row["local_relative_path"].strip())
        if local_path.is_absolute() or ".." in local_path.parts:
            raise ValueError(f"Unsafe local path on manifest row {row_number}: {local_path}")
        records.append(
            PublicFile(
                dataset_id=row["dataset_id"].strip(),
                file_id=row["file_id"].strip(),
                url=row["exact_file_url"].strip(),
                size_bytes=int(row["size_bytes"]),
                sha256=digest,
                access_class=row["access_class"].strip().lower(),
                local_relative_path=local_path,
            )
        )
    return records


def select_records(records: list[PublicFile], datasets: list[str]) -> list[PublicFile]:
    if not datasets:
        return records
    requested = set(datasets)
    selected = [record for record in records if record.dataset_id in requested]
    missing = sorted(requested.difference(record.dataset_id for record in selected))
    if missing:
        raise ValueError(f"Unknown dataset_id: {', '.join(missing)}")
    return selected


def resolve_destination(project_root: Path, record: PublicFile) -> Path:
    destination = (project_root / record.local_relative_path).resolve()
    allowed_root = (project_root / "data" / "raw").resolve()
    if destination != allowed_root and allowed_root not in destination.parents:
        raise ValueError(f"Destination must remain under data/raw: {destination}")
    return destination


def verify(path: Path, record: PublicFile) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing"
    observed_size = path.stat().st_size
    if observed_size != record.size_bytes:
        return False, f"size mismatch ({observed_size} != {record.size_bytes})"
    observed_digest = sha256_file(path)
    if observed_digest != record.sha256:
        return False, f"SHA-256 mismatch ({observed_digest})"
    return True, "verified"


def download(record: PublicFile, destination: Path, force: bool) -> None:
    if record.access_class != "public":
        raise PermissionError(
            f"{record.dataset_id}/{record.file_id} is {record.access_class}; "
            "controlled inputs must be retrieved through their authorized repository"
        )
    valid, message = verify(destination, record)
    if valid:
        print(f"verified\t{record.dataset_id}\t{record.file_id}\t{destination}")
        return
    if destination.exists() and not force:
        raise FileExistsError(
            f"Refusing to replace {destination}: {message}. Re-run with --force after review."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".partial",
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
        try:
            request = urllib.request.Request(
                record.url,
                headers={"User-Agent": "chemokine-cart-reproducibility/0.1"},
            )
            with urllib.request.urlopen(request) as response:
                shutil.copyfileobj(response, temporary, length=CHUNK_SIZE)
            temporary.flush()
            os.fsync(temporary.fileno())
        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise

    valid, message = verify(temporary_path, record)
    if not valid:
        temporary_path.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded file failed verification: {message}")
    os.replace(temporary_path, destination)
    print(f"downloaded\t{record.dataset_id}\t{record.file_id}\t{destination}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch public study files recorded with exact sizes and SHA-256 digests."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("config/public_datasets.tsv"),
        help="Public-file manifest (default: config/public_datasets.tsv).",
    )
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Dataset accession to select; repeat for multiple accessions.",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root used to resolve local_relative_path.",
    )
    parser.add_argument(
        "--verify-only", action="store_true", help="Verify existing files without downloading."
    )
    parser.add_argument(
        "--list", action="store_true", help="List selected files without downloading."
    )
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing file that fails verification."
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    records = select_records(read_manifest(args.manifest), args.dataset)
    project_root = args.project_root.resolve()

    failures = 0
    for record in records:
        destination = resolve_destination(project_root, record)
        if args.list:
            print(
                f"{record.dataset_id}\t{record.file_id}\t{record.size_bytes}\t"
                f"{record.sha256}\t{record.url}\t{destination}"
            )
            continue
        if args.verify_only:
            valid, message = verify(destination, record)
            print(
                f"{'verified' if valid else 'failed'}\t{record.dataset_id}\t{record.file_id}\t{message}"
            )
            failures += int(not valid)
            continue
        try:
            download(record, destination, args.force)
        except Exception as error:
            failures += 1
            print(f"failed\t{record.dataset_id}\t{record.file_id}\t{error}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
