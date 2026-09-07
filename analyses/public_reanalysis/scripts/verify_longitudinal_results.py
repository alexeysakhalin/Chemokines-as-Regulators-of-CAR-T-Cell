#!/usr/bin/env python3
"""Verify a longitudinal run against frozen data, code, and output digests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_ROOT))

from longitudinal_extension.verification import (  # noqa: E402
    load_manifest,
    verify_checksums,
    verify_observed_manifest,
    verify_results,
    verify_working_tree_manifest,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument(
        "--module-root",
        type=Path,
        default=MODULE_ROOT,
        help="Repository analysis-module root used to resolve frozen provenance paths",
    )
    args = parser.parse_args()
    verify_checksums(args.expected)
    expected_manifest = load_manifest(args.expected / "analysis_manifest.json")
    verify_working_tree_manifest(expected_manifest, args.module_root)
    verify_results(args.observed, args.expected)
    verify_observed_manifest(args.observed, expected_manifest)
    count = len(expected_manifest["canonical_outputs"])
    print(f"Verified current provenance and {count} byte-identical canonical outputs")


if __name__ == "__main__":
    main()
