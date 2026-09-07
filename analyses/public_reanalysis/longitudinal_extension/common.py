"""Small deterministic utilities shared by parsers and analysis code."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_two_sided_sign_test(increases: int, decreases: int) -> float:
    """Two-sided exact sign-test P value after discarding zero differences."""
    if increases < 0 or decreases < 0:
        raise ValueError("Direction counts must be non-negative")
    n = increases + decreases
    if n == 0:
        return float("nan")
    tail = min(increases, decreases)
    probability = 2.0 * sum(math.comb(n, index) for index in range(tail + 1)) / (2**n)
    return float(min(probability, 1.0))


def holm_adjust(values: pd.Series) -> pd.Series:
    """Holm family-wise adjustment, preserving the original index."""
    if values.empty:
        return values.astype(float)
    numeric = values.astype(float)
    if numeric.isna().any():
        raise ValueError("Holm adjustment does not accept missing P values")
    order = np.argsort(numeric.to_numpy(), kind="mergesort")
    ranked = numeric.to_numpy()[order]
    adjusted_ranked = np.maximum.accumulate(ranked * (len(ranked) - np.arange(len(ranked))))
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return pd.Series(adjusted, index=values.index, dtype=float)


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment, preserving the original index."""
    if values.empty:
        return values.astype(float)
    numeric = values.astype(float)
    if numeric.isna().any():
        raise ValueError("BH adjustment does not accept missing P values")
    order = np.argsort(numeric.to_numpy(), kind="mergesort")
    ranked = numeric.to_numpy()[order]
    adjusted_ranked = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return pd.Series(adjusted, index=values.index, dtype=float)


def quantile(values: pd.Series, probability: float) -> float:
    """Return a deterministic linear sample quantile."""
    return float(values.astype(float).quantile(probability, interpolation="linear"))


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Write a stable UTF-8 TSV with LF newlines and explicit missing values."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def write_json(payload: dict[str, Any], path: Path) -> None:
    """Write a stable JSON document."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
