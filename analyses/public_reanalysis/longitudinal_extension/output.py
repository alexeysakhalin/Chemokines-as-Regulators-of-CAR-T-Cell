"""Deterministic output and report generation for the longitudinal extension."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .common import sha256, write_json, write_tsv

OUTPUT_TABLES = (
    "sample_aggregates.tsv",
    "marker_aggregates.tsv",
    "contrast_inclusion.tsv",
    "contrast_results.tsv",
    "marker_contrast_results.tsv",
    "leave_one_out.tsv",
)
CANONICAL_OUTPUTS = (*OUTPUT_TABLES, "REPORT.md")
PYTHON_REQUIREMENT = ">=3.11,<3.12"


def _is_raw_input(relative: Path) -> bool:
    """Return whether a repository-relative path denotes downloaded source data."""
    return len(relative.parts) >= 2 and relative.parts[:2] == ("data", "raw")


def _canonical_input_key(path: Path, root: Path, digest: str) -> str:
    """Create a checkout-independent manifest key for an analysis input.

    Version-controlled inputs retain their repository-relative path. Downloaded files and
    inputs supplied from a cache outside the checkout are content-addressed, so invoking the
    workflow with ``/tmp/...`` instead of ``data/raw/...`` cannot change the manifest bytes.
    """
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(root.resolve())
    except ValueError:
        return f"source/sha256/{digest}"
    if _is_raw_input(relative):
        return f"source/sha256/{digest}"
    return relative.as_posix()


def _canonical_code_key(path: Path, root: Path) -> str:
    """Return a repository-relative code path, rejecting undeclared external code."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(f"Code path is outside the manifest root: {path}") from error


def _digest_map(paths: list[Path], root: Path, *, inputs: bool) -> dict[str, str]:
    """Hash files under stable logical keys and reject ambiguous key collisions."""
    entries: dict[str, str] = {}
    for path in sorted(paths, key=lambda item: str(item)):
        if not path.is_file():
            raise ValueError(f"Manifest input is not a file: {path}")
        digest = sha256(path)
        key = (
            _canonical_input_key(path, root, digest) if inputs else _canonical_code_key(path, root)
        )
        if key in entries and entries[key] != digest:
            raise ValueError(f"Manifest key collision for {key}")
        entries[key] = digest
    return dict(sorted(entries.items()))


def _format_p(value: object) -> str:
    if value is None or pd.isna(value):
        return "not adjusted"
    numeric = float(value)
    if numeric < 0.001:
        return f"{numeric:.2e}"
    return f"{numeric:.3f}"


def _direction_sentence(row: pd.Series) -> str:
    n = int(row["n_pairs"])
    noun = "patient" if n == 1 else "patients"
    return (
        f"{n} {noun}; median change {float(row['median_change_percentage_points']):+.2f} "
        f"percentage points; {int(row['n_increase'])} increased, "
        f"{int(row['n_decrease'])} decreased, {int(row['n_unchanged'])} unchanged; "
        f"exact sign-test P={_format_p(row['sign_test_p'])}"
    )


def build_report(results: pd.DataFrame) -> str:
    """Create manuscript-drafting text from numeric results without hidden inference."""
    lines = [
        "# Longitudinal CXCR6 public-cohort extension",
        "",
        "## Scope",
        "",
        "The patient is the independent unit. Each cohort is summarized separately; cells are "
        "never treated as biological replicates, and no pooled cross-study P value is calculated.",
        "",
        "## Locked CXCR6-detection contrasts",
        "",
    ]
    for _, row in results.iterrows():
        lines.append(f"- **{row['dataset']} — {row['contrast']}:** {_direction_sentence(row)}.")
    early = results[results["family_id"].eq("early_cxcr6_fraction_three_cohort")]
    if len(early) == 3:
        significant = int((early["holm_p"] < 0.05).sum())
        family_sentence = (
            "none of the three early tests had Holm-adjusted P<0.05"
            if significant == 0
            else f"{significant} of 3 cohort-specific tests had Holm-adjusted P<0.05"
        )
        lines.extend(
            [
                "",
                "## Cross-cohort interpretation",
                "",
                f"In the finalized exploratory three-test early comparison family, {family_sentence}. "
                "The cohort medians were "
                "positive in GSE197268 and GSE235760 and negative in GSE162975. These data are "
                "consistent with temporal and context dependence of CXCR6 transcript detection, not a universal "
                "early increase.",
            ]
        )
    lines.extend(
        [
            "",
            "## Boundaries",
            "",
            "These transcript-level observational comparisons do not establish trafficking, "
            "chemotaxis, retention, egress, a CXCL16 protein gradient, treatment efficacy, "
            "neurotoxicity, or causality. Product-stratified estimates are descriptive and were not "
            "used to claim a time-by-product interaction.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(
    output_dir: Path,
    *,
    samples: pd.DataFrame,
    markers: pd.DataFrame,
    inclusion: pd.DataFrame,
    results: pd.DataFrame,
    marker_results: pd.DataFrame,
    leave_one_out: pd.DataFrame,
    input_paths: list[Path],
    code_paths: list[Path],
    seed: int,
    manifest_root: Path | None = None,
    source_manifest_path: Path | None = None,
    audit_tables: dict[str, pd.DataFrame] | None = None,
) -> None:
    """Write tables, report, run manifest, and an internal SHA-256 manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in (
        ("sample_aggregates.tsv", samples),
        ("marker_aggregates.tsv", markers),
        ("contrast_inclusion.tsv", inclusion),
        ("contrast_results.tsv", results),
        ("marker_contrast_results.tsv", marker_results),
        ("leave_one_out.tsv", leave_one_out),
    ):
        write_tsv(frame, output_dir / name)
    for name, frame in sorted((audit_tables or {}).items()):
        if Path(name).name != name or not name.endswith(".tsv"):
            raise ValueError(f"Audit table must be a plain TSV filename: {name}")
        write_tsv(frame, output_dir / name)
    (output_dir / "REPORT.md").write_text(build_report(results), encoding="utf-8", newline="\n")
    root = (manifest_root or Path.cwd()).resolve()
    canonical_names = (*CANONICAL_OUTPUTS, *sorted((audit_tables or {}).keys()))
    canonical_outputs = {
        name: sha256(output_dir / name) for name in canonical_names if (output_dir / name).is_file()
    }
    manifest = {
        "analysis": "patient-level longitudinal CXCR6 public-cohort extension",
        "seed": seed,
        "python_requirement": PYTHON_REQUIREMENT,
        "inputs": _digest_map(input_paths, root, inputs=True),
        "code": _digest_map(code_paths, root, inputs=False),
        "canonical_outputs": canonical_outputs,
        "independent_unit": "patient",
        "cell_level_p_values": False,
        "cross_study_pooling": False,
    }
    if source_manifest_path is not None:
        manifest["source_manifest"] = {
            "path": _canonical_code_key(source_manifest_path, root),
            "sha256": sha256(source_manifest_path),
        }
    if audit_tables:
        manifest["audit_outputs"] = {
            name: sha256(output_dir / name) for name in sorted(audit_tables)
        }
    write_json(manifest, output_dir / "analysis_manifest.json")
    files = sorted(
        path for path in output_dir.iterdir() if path.is_file() and path.name != "SHA256SUMS"
    )
    checksum_lines = [f"{sha256(path)}  {path.name}" for path in files]
    (output_dir / "SHA256SUMS").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8", newline="\n"
    )


def assert_finite_results(results: pd.DataFrame) -> None:
    """Reject incomplete numeric results before they are frozen."""
    numeric = results[
        [
            "n_pairs",
            "n_increase",
            "n_decrease",
            "n_unchanged",
            "median_change_fraction",
            "sign_test_p",
        ]
    ].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("Contrast results contain non-finite required values")
