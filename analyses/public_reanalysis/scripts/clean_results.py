#!/usr/bin/env python3
"""Remove only the declared local output directory."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=Path, required=True)
    args = parser.parse_args()
    target = args.path.resolve()
    expected_suffix = Path("results/local")
    if tuple(target.parts[-2:]) != tuple(expected_suffix.parts):
        raise ValueError(f"Refusing to remove undeclared path: {target}")
    if target.exists():
        shutil.rmtree(target)


if __name__ == "__main__":
    main()
