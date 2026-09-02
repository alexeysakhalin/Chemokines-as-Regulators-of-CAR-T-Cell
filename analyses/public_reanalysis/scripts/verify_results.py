#!/usr/bin/env python3
"""Compare scientific tables and figures with the frozen reviewed outputs."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    args = parser.parse_args()
    expected_files = sorted(
        path.relative_to(args.expected)
        for path in args.expected.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    )
    observed_files = sorted(
        path.relative_to(args.observed)
        for path in args.observed.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    )
    missing = [path for path in expected_files if not (args.observed / path).is_file()]
    if missing:
        raise FileNotFoundError(f"Observed run is missing: {missing}")
    unexpected = sorted(set(observed_files).difference(expected_files))
    if unexpected:
        raise ValueError(f"Observed run contains unexpected files: {unexpected}")
    mismatches = [
        path
        for path in expected_files
        if digest(args.observed / path) != digest(args.expected / path)
    ]
    if mismatches:
        raise ValueError(f"Outputs differ from frozen results: {mismatches}")
    print(f"Verified {len(expected_files)} byte-identical outputs")


if __name__ == "__main__":
    main()
