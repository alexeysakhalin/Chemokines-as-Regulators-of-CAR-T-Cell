"""Verify byte-identical scientific outputs across two clean pipeline executions."""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"
MANIFEST_PATH = PROJECT_ROOT / "data" / "demo" / "manifest.tsv"
COMPARED_DIRECTORIES = ("tables", "figures")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _inventory(root: Path) -> dict[str, str]:
    inventory: dict[str, str] = {}
    for directory_name in COMPARED_DIRECTORIES:
        directory = root / directory_name
        if not directory.is_dir():
            raise RuntimeError(f"Expected output directory is missing: {directory}")
        for path in sorted(directory.rglob("*")):
            if not path.is_file():
                continue
            if path.stat().st_size == 0:
                raise RuntimeError(f"Output file is empty: {path}")
            inventory[path.relative_to(root).as_posix()] = _sha256(path)
    return inventory


def _run(output_dir: Path, *, random_seed: int) -> None:
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONHASHSEED": str(random_seed),
            "LANG": "C.UTF-8",
            "LC_ALL": "C.UTF-8",
            "TZ": "UTC",
            "SOURCE_DATE_EPOCH": "1788220800",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "MPLBACKEND": "Agg",
        }
    )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "chemokine_cart.pipeline",
            "--stage",
            "all",
            "--config",
            str(CONFIG_PATH),
            "--manifest",
            str(MANIFEST_PATH),
            "--output-dir",
            str(output_dir),
            "--random-seed",
            str(random_seed),
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        check=True,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare scientific outputs from two clean pipeline executions."
    )
    parser.add_argument("--random-seed", type=int, default=20260901)
    args = parser.parse_args(argv)

    if not MANIFEST_PATH.is_file():
        raise SystemExit(
            "Demonstration manifest is missing. Run "
            "'chemokine-cart make-demo --output-dir data/demo --force' first."
        )

    with tempfile.TemporaryDirectory(prefix=".determinism-", dir=PROJECT_ROOT) as temporary:
        temporary_root = Path(temporary)
        first = temporary_root / "first"
        second = temporary_root / "second"
        _run(first, random_seed=args.random_seed)
        _run(second, random_seed=args.random_seed)
        first_inventory = _inventory(first)
        second_inventory = _inventory(second)

    if first_inventory != second_inventory:
        paths = sorted(set(first_inventory) | set(second_inventory))
        changed = [
            path for path in paths if first_inventory.get(path) != second_inventory.get(path)
        ]
        raise SystemExit("Non-deterministic scientific outputs: " + ", ".join(changed))

    print(f"Verified {len(first_inventory)} byte-identical scientific output files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
