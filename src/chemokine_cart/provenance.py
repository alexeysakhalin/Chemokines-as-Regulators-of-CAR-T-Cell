"""Deterministic provenance records for every analysis stage."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from collections.abc import Iterable, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

TRACKED_PACKAGES: tuple[str, ...] = (
    "chemokine-cart",
    "anndata",
    "h5py",
    "matplotlib",
    "numpy",
    "pandas",
    "PyYAML",
    "scanpy",
    "scikit-learn",
    "scipy",
    "seaborn",
    "snakemake",
    "statsmodels",
)


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Calculate a SHA-256 digest without loading the complete file."""

    resolved = Path(path)
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path, *, root: str | Path | None = None) -> dict[str, object]:
    """Return a stable size-and-checksum record for one file."""

    resolved = Path(path).resolve()
    display = resolved
    if root is not None:
        root_path = Path(root).resolve()
        try:
            display = resolved.relative_to(root_path)
        except ValueError:
            display = resolved
    return {
        "path": display.as_posix(),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def package_versions(packages: Sequence[str] = TRACKED_PACKAGES) -> dict[str, str | None]:
    """Return installed versions; unavailable optional packages are recorded as null."""

    versions: dict[str, str | None] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def git_revision(project_root: str | Path) -> dict[str, object]:
    """Read the local revision when available without requiring a Git checkout."""

    root = Path(project_root)
    record: dict[str, object] = {"commit": None, "dirty": None}
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError):
        return record
    record["commit"] = commit
    record["dirty"] = bool(status.strip())
    return record


def runtime_record(project_root: str | Path, *, seed: int) -> dict[str, object]:
    """Collect the reproducibility-relevant runtime state."""

    return {
        "created_utc": datetime.now(UTC).isoformat(),
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "byte_order": sys.byteorder,
        "random_seed": int(seed),
        "thread_environment": {
            name: os.environ.get(name)
            for name in (
                "PYTHONHASHSEED",
                "LANG",
                "LC_ALL",
                "TZ",
                "SOURCE_DATE_EPOCH",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
        "git": git_revision(project_root),
        "packages": package_versions(),
    }


def write_json(path: str | Path, payload: Mapping[str, object]) -> Path:
    """Write canonical UTF-8 JSON through an atomic replacement."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, target)
    return target


def write_stage_marker(
    output_dir: str | Path,
    *,
    stage: str,
    outputs: Iterable[str | Path],
    inputs: Iterable[str | Path] = (),
    metadata: Mapping[str, object] | None = None,
) -> Path:
    """Create a non-empty stage marker only after all declared outputs exist."""

    root = Path(output_dir).resolve()
    output_paths = [Path(path).resolve() for path in outputs]
    missing = [path for path in output_paths if not path.is_file()]
    if missing:
        joined = ", ".join(path.as_posix() for path in missing)
        raise FileNotFoundError(f"stage {stage!r} is missing declared outputs: {joined}")

    input_paths = [Path(path).resolve() for path in inputs if Path(path).is_file()]
    payload: dict[str, object] = {
        "stage": stage,
        "status": "completed",
        "inputs": [file_record(path, root=root.parent) for path in input_paths],
        "outputs": [file_record(path, root=root) for path in output_paths],
    }
    if metadata:
        payload["metadata"] = dict(metadata)
    marker = root / ".workflow" / f"{stage.replace('-', '_')}.done"
    return write_json(marker, payload)
