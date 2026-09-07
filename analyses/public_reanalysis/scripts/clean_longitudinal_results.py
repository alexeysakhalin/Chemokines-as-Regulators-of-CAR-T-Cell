#!/usr/bin/env python3
"""Remove only the explicitly selected local longitudinal result directory."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=Path, required=True)
    args = parser.parse_args()
    target = args.path.resolve()
    expected_suffix = Path("results/longitudinal_extension/local")
    if not target.parts[-len(expected_suffix.parts) :] == expected_suffix.parts:
        raise ValueError(f"Refusing to remove unexpected path: {target}")
    if target.is_dir():
        shutil.rmtree(target)


if __name__ == "__main__":
    main()
