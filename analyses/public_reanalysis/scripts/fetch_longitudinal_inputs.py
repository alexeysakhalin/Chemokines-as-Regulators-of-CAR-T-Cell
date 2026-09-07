#!/usr/bin/env python3
"""Download or verify the checksum-pinned longitudinal public inputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_ROOT))

from longitudinal_extension.sources import fetch_manifest, load_source_manifest  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("config/longitudinal_sources.tsv"))
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--list", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = args.root / args.manifest
    if args.list:
        for row in load_source_manifest(manifest):
            print(
                "\t".join(
                    [row["dataset"], row["local_path"], row["bytes"], row["sha256"], row["role"]]
                )
            )
        return
    fetch_manifest(manifest, args.root, verify_only=args.verify_only)


if __name__ == "__main__":
    main()
