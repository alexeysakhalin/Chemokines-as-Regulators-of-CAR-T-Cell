"""Stage-oriented analysis pipeline used by the command line and Snakemake."""

from __future__ import annotations

import argparse
import io
import json
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import functional, single_cell, spatial
from .io import (
    DataValidationError,
    load_config,
    read_tsv,
    resolve_project_path,
    validate_manifest,
    validate_table,
)
from .provenance import file_record, runtime_record, sha256_file, write_json, write_stage_marker

plt.switch_backend("Agg")
matplotlib.rcParams["svg.hashsalt"] = "chemokine-cart-v0.1"

STAGES: tuple[str, ...] = (
    "validate",
    "single-cell",
    "spatial",
    "functional",
    "figures",
    "report",
)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool | np.bool_):
        return bool(value)
    if isinstance(value, int | np.integer) and value in {0, 1}:
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n"}:
        return False
    raise ValueError(f"value is not Boolean-like: {value!r}")


def _effective_seed(config: Mapping[str, Any]) -> int:
    """Return the effective seed stored in the resolved project config."""

    return int(config.get("project", {}).get("random_seed", 20260901))


def _project_context(
    config_path: str | Path,
    manifest_path: str | Path,
    output_dir: str | Path,
    random_seed: int | None = None,
) -> tuple[dict[str, Any], Path, Path, Path]:
    config_file = Path(config_path).resolve()
    config = load_config(config_file)
    if random_seed is not None:
        project_config = config.get("project")
        if not isinstance(project_config, dict):
            raise DataValidationError(["config: project must be a mapping"])
        project_config["random_seed"] = int(random_seed)
    project_root = Path(config["_project_root"]).resolve()
    manifest_file = Path(manifest_path)
    if not manifest_file.is_absolute():
        manifest_file = resolve_project_path(project_root, manifest_file.as_posix())
    manifest_file = manifest_file.resolve()
    output = Path(output_dir)
    if not output.is_absolute():
        output = resolve_project_path(project_root, output.as_posix())
    return config, project_root, manifest_file, output.resolve()


def _manifest_frame(manifest_path: Path) -> pd.DataFrame:
    _, rows = read_tsv(manifest_path)
    manifest = pd.DataFrame(rows)
    if manifest.empty:
        raise DataValidationError([f"{manifest_path}: manifest contains no records"])
    required = {"record_id", "data_origin", "evidence_eligible"}
    missing = sorted(required.difference(manifest.columns))
    if missing:
        raise DataValidationError(
            [f"{manifest_path}: required manifest columns missing: {', '.join(missing)}"]
        )
    try:
        manifest["evidence_eligible"] = manifest["evidence_eligible"].map(_as_bool)
    except ValueError as exc:
        raise DataValidationError(
            [f"{manifest_path}: invalid evidence_eligible value: {exc}"]
        ) from exc
    return manifest


def _is_synthetic_only(manifest: pd.DataFrame) -> bool:
    return set(manifest["data_origin"].astype(str).str.lower()) == {"synthetic"}


def _enforce_evidence_policy(config: Mapping[str, Any], manifest: pd.DataFrame) -> None:
    """Block non-synthetic scientific stages when any input is ineligible."""

    required = {"record_id", "data_origin", "evidence_eligible"}
    missing_columns = sorted(required.difference(manifest.columns))
    if missing_columns:
        raise DataValidationError(
            [f"manifest required evidence columns missing: {', '.join(missing_columns)}"]
        )

    origins = manifest["data_origin"].astype("string").str.strip().str.lower()
    missing_origin = manifest["data_origin"].isna() | origins.eq("").fillna(True)
    if missing_origin.any():
        records = manifest.loc[missing_origin, "record_id"].astype(str).tolist()
        raise DataValidationError(
            [
                f"scientific stages require an explicit data_origin; missing for: {', '.join(records[:10])}"
            ]
        )
    allowed_origins = {"experimental", "public", "synthetic"}
    invalid_origin = ~origins.isin(allowed_origins)
    if invalid_origin.any():
        values = sorted(set(origins.loc[invalid_origin].astype(str)))
        raise DataValidationError(
            [f"scientific stages received unsupported data_origin value(s): {', '.join(values)}"]
        )

    missing_eligibility = manifest["evidence_eligible"].isna() | manifest[
        "evidence_eligible"
    ].astype(str).str.strip().eq("")
    if missing_eligibility.any():
        records = manifest.loc[missing_eligibility, "record_id"].astype(str).tolist()
        raise DataValidationError(
            [
                "scientific stages require an explicit evidence_eligible value for every "
                f"manifest record; missing for: {', '.join(records[:10])}"
            ]
        )
    try:
        eligible = manifest["evidence_eligible"].map(_as_bool)
    except ValueError as exc:
        raise DataValidationError(
            [f"manifest has an invalid evidence_eligible value: {exc}"]
        ) from exc

    synthetic_rows = origins.eq("synthetic")
    if eligible.loc[synthetic_rows].any():
        raise DataValidationError(
            ["synthetic records must be explicitly marked evidence_eligible=false"]
        )
    if synthetic_rows.all():
        return
    if synthetic_rows.any():
        raise DataValidationError(
            [
                "synthetic records cannot be mixed into a scientific run; run the synthetic-only "
                "software demonstration separately"
            ]
        )
    evidence_config = config.get("evidence", {})
    if not isinstance(evidence_config, Mapping):
        raise DataValidationError(["config: evidence must be a mapping"])
    require_eligible = _as_bool(
        evidence_config.get("require_evidence_eligible_for_scientific_outputs", True)
    )
    if not require_eligible:
        return
    ineligible = manifest.loc[~eligible, "record_id"].astype(str).tolist()
    if ineligible:
        examples = ", ".join(ineligible[:10])
        suffix = "" if len(ineligible) <= 10 else f" (+{len(ineligible) - 10} more)"
        raise DataValidationError(
            [
                "scientific stages require evidence_eligible=true for every non-synthetic "
                f"manifest record; ineligible record(s): {examples}{suffix}"
            ]
        )


def _source_paths(
    manifest: pd.DataFrame,
    project_root: Path,
    *,
    modalities: Iterable[str],
) -> list[tuple[pd.Series, Path]]:
    selected = manifest.loc[manifest["modality"].isin(tuple(modalities))]
    if selected.empty:
        names = ", ".join(sorted(set(modalities)))
        raise DataValidationError([f"manifest contains no local records for modality: {names}"])
    records: list[tuple[pd.Series, Path]] = []
    for _, record in selected.sort_values("record_id", kind="mergesort").iterrows():
        source = str(record["source_path"]).strip()
        if not source:
            raise DataValidationError(
                [
                    f"record {record['record_id']}: accession-only inputs must be materialized "
                    "locally before analysis"
                ]
            )
        records.append((record, resolve_project_path(project_root, source, must_exist=True)))
    return records


def _require_current_stage_marker(
    output_dir: Path,
    *,
    stage: str,
    config_path: str | Path,
    manifest_path: Path,
    synthetic_only: bool,
    random_seed: int,
) -> None:
    """Reject missing, stale, or modified upstream stage outputs."""

    marker_path = output_dir / ".workflow" / f"{stage.replace('-', '_')}.done"
    if not marker_path.is_file():
        raise DataValidationError(
            [f"stage {stage!r} must complete for the current inputs before downstream execution"]
        )
    try:
        payload = json.loads(marker_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise DataValidationError([f"invalid upstream stage marker: {marker_path}"]) from exc
    if payload.get("stage") != stage or payload.get("status") != "completed":
        raise DataValidationError([f"invalid upstream stage marker state: {marker_path}"])

    marker_root = output_dir.parent.resolve()

    def resolve_record_path(record: Mapping[str, Any], *, output: bool) -> Path:
        path = Path(str(record.get("path", "")))
        if path.is_absolute():
            return path.resolve()
        return ((output_dir if output else marker_root) / path).resolve()

    input_records = payload.get("inputs", [])
    if not isinstance(input_records, list):
        raise DataValidationError([f"invalid input records in upstream marker: {marker_path}"])
    for current_input in (Path(config_path).resolve(), manifest_path.resolve()):
        matches = [
            record
            for record in input_records
            if isinstance(record, Mapping)
            and resolve_record_path(record, output=False) == current_input
        ]
        if len(matches) != 1 or matches[0].get("sha256") != sha256_file(current_input):
            raise DataValidationError(
                [
                    f"stage {stage!r} was not produced from the current "
                    f"{current_input.name}; rerun upstream stages"
                ]
            )

    output_records = payload.get("outputs", [])
    if not isinstance(output_records, list) or not output_records:
        raise DataValidationError([f"upstream marker has no declared outputs: {marker_path}"])
    for record in output_records:
        if not isinstance(record, Mapping):
            raise DataValidationError([f"invalid output record in upstream marker: {marker_path}"])
        path = resolve_record_path(record, output=True)
        if not path.is_file() or record.get("sha256") != sha256_file(path):
            raise DataValidationError(
                [f"upstream stage {stage!r} output is missing or modified: {path}"]
            )

    metadata = payload.get("metadata", {})
    if not isinstance(metadata, Mapping):
        raise DataValidationError([f"invalid metadata in upstream marker: {marker_path}"])
    if metadata.get("synthetic_only") is not synthetic_only:
        raise DataValidationError(
            [f"upstream stage {stage!r} synthetic/evidence status does not match current manifest"]
        )
    if int(metadata.get("random_seed", -1)) != int(random_seed):
        raise DataValidationError(
            [f"upstream stage {stage!r} random seed does not match the current run"]
        )


def _load_tables(
    manifest: pd.DataFrame,
    project_root: Path,
    *,
    modalities: Iterable[str],
) -> pd.DataFrame:
    selected_records = _source_paths(manifest, project_root, modalities=modalities)
    records_by_source: dict[Path, list[pd.Series]] = {}
    for record, source in selected_records:
        records_by_source.setdefault(source, []).append(record)

    frames: list[pd.DataFrame] = []
    identity_columns = (
        "sample_id",
        "patient_id",
        "section_id",
        "timepoint_id",
        "data_origin",
    )
    attached_columns = (
        "experimental_arm",
        "compartment",
        "batch_id",
        "reference_version",
        "protocol_version",
        "assay_panel_version",
        "condition",
        "comparator",
    )
    for source, records in sorted(records_by_source.items(), key=lambda item: item[0].as_posix()):
        source_frame = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
        synthetic_multi = (
            len(records) == 1
            and str(records[0]["data_origin"]).lower() == "synthetic"
            and str(records[0]["sample_id"]).upper() == "SYNTH_COHORT"
        )
        if not synthetic_multi:
            if "sample_id" not in source_frame.columns:
                raise DataValidationError([f"{source}: sample_id is required for manifest linkage"])
            manifest_samples = [str(record["sample_id"]) for record in records]
            if len(manifest_samples) != len(set(manifest_samples)):
                raise DataValidationError(
                    [f"{source}: duplicate manifest records for the same sample_id"]
                )
            observed_samples = set(source_frame["sample_id"].astype(str))
            expected_samples = set(manifest_samples)
            if observed_samples != expected_samples:
                missing = sorted(expected_samples.difference(observed_samples))
                unregistered = sorted(observed_samples.difference(expected_samples))
                details: list[str] = []
                if missing:
                    details.append("missing manifest samples: " + ", ".join(missing))
                if unregistered:
                    details.append("unregistered table samples: " + ", ".join(unregistered))
                raise DataValidationError([f"{source}: " + "; ".join(details)])

        for record in records:
            if synthetic_multi:
                frame = source_frame.copy()
            else:
                frame = source_frame.loc[
                    source_frame["sample_id"].astype(str) == str(record["sample_id"])
                ].copy()
            if frame.empty:
                raise DataValidationError(
                    [
                        f"record {record['record_id']}: no rows matched sample_id {record['sample_id']}"
                    ]
                )

            for column in identity_columns:
                if column not in frame.columns:
                    raise DataValidationError(
                        [f"{source}: required identity column {column!r} missing"]
                    )
                if synthetic_multi and column in {
                    "sample_id",
                    "patient_id",
                    "section_id",
                    "timepoint_id",
                }:
                    continue
                observed = set(frame[column].astype(str))
                expected = {str(record[column])}
                if observed != expected:
                    raise DataValidationError(
                        [
                            f"record {record['record_id']}: table {column} values "
                            f"{sorted(observed)!r} do not match manifest value "
                            f"{record[column]!r}"
                        ]
                    )

            if "modality" in frame.columns:
                observed_modalities = set(frame["modality"].astype(str))
                expected_modality = {str(record["modality"])}
                if observed_modalities != expected_modality:
                    raise DataValidationError(
                        [
                            f"record {record['record_id']}: table modality values "
                            f"{sorted(observed_modalities)!r} do not match manifest value "
                            f"{record['modality']!r}"
                        ]
                    )

            for column in attached_columns:
                manifest_value = str(record.get(column, ""))
                if column in frame.columns and not synthetic_multi:
                    observed = set(frame[column].astype(str))
                    if observed != {manifest_value}:
                        raise DataValidationError(
                            [
                                f"record {record['record_id']}: table {column} values "
                                f"{sorted(observed)!r} do not match manifest value "
                                f"{manifest_value!r}"
                            ]
                        )
                if column not in frame.columns or not synthetic_multi:
                    frame[column] = manifest_value
            frame["_manifest_record_id"] = str(record["record_id"])
            frame["_manifest_modality"] = str(record["modality"])
            frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False)


def _write_table(frame: pd.DataFrame, path: Path, *, sort_by: Sequence[str] = ()) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    output = frame.copy()
    valid_sort = [column for column in sort_by if column in output.columns]
    if valid_sort and not output.empty:
        output = output.sort_values(valid_sort, kind="mergesort").reset_index(drop=True)
    output.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        float_format="%.10g",
        na_rep="",
    )
    return path


def _bh_adjust(values: Sequence[float]) -> np.ndarray:
    p_values = np.asarray(values, dtype=float)
    if p_values.size == 0:
        return p_values
    order = np.argsort(p_values, kind="mergesort")
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def _stage_validate(
    config: Mapping[str, Any],
    project_root: Path,
    manifest_path: Path,
    output_dir: Path,
) -> list[Path]:
    manifest_config = config.get("manifest", {})
    schema_path = resolve_project_path(
        project_root,
        str(manifest_config.get("schema", "resources/manifest.schema.tsv")),
        must_exist=True,
    )
    manifest_report = validate_manifest(
        manifest_path,
        schema_path,
        project_root,
        check_files=True,
        verify_checksums=True,
    )
    manifest = _manifest_frame(manifest_path)
    schema_config = config.get("schemas", {})
    table_reports: list[dict[str, object]] = []
    for record, source in _source_paths(
        manifest,
        project_root,
        modalities=tuple(manifest["modality"].unique()),
    ):
        modality = str(record["modality"])
        schema_value = schema_config.get(modality)
        if schema_value is None and modality in {"scrna", "citeseq"}:
            schema_value = schema_config.get("single_cell")
        if schema_value is None:
            raise DataValidationError([f"config has no table schema for modality {modality!r}"])
        table_schema = resolve_project_path(project_root, str(schema_value), must_exist=True)
        report = validate_table(source, table_schema, allow_extra_columns=True)
        row = report.as_dict()
        row["record_id"] = str(record["record_id"])
        row["modality"] = modality
        row["sha256"] = sha256_file(source)
        table_reports.append(row)

    single_cell_modalities = tuple(
        modality for modality in ("scrna", "citeseq") if modality in set(manifest["modality"])
    )
    if single_cell_modalities:
        single_table = _load_tables(manifest, project_root, modalities=single_cell_modalities)
        single_config = config.get("single_cell", {})
        single_cell.validate_required_columns(
            single_table,
            (
                str(single_config.get("car_column", "car_positive")),
                str(single_config.get("origin_column", "t_cell_origin")),
                str(single_config.get("state_column", "cell_state")),
            ),
            table_name="single_cell",
        )

    payload = {
        "status": "valid",
        "manifest": manifest_report.as_dict(),
        "records": table_reports,
        "synthetic_only": _is_synthetic_only(manifest),
        "evidence_eligible_records": int(manifest["evidence_eligible"].sum()),
        "note": (
            "Synthetic inputs are eligible for software testing only."
            if _is_synthetic_only(manifest)
            else "Evidence eligibility is controlled record by record in the manifest."
        ),
    }
    output = write_json(output_dir / "tables" / "validation_report.json", payload)
    write_stage_marker(
        output_dir,
        stage="validate",
        outputs=[output],
        inputs=[
            config["_config_path"],
            manifest_path,
            *[
                path
                for _, path in _source_paths(
                    manifest,
                    project_root,
                    modalities=tuple(manifest["modality"].unique()),
                )
            ],
        ],
        metadata={
            "random_seed": _effective_seed(config),
            "synthetic_only": payload["synthetic_only"],
        },
    )
    return [output]


def _consistent_cell_metadata(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    key_columns: Sequence[str] = ("sample_id", "cell_id"),
) -> pd.DataFrame:
    conflicts: list[str] = []
    for column in columns:
        counts = table.groupby(list(key_columns), observed=True)[column].nunique(dropna=False)
        if (counts > 1).any():
            conflicts.append(column)
    if conflicts:
        raise DataValidationError(
            ["cell-level metadata conflict across modalities: " + ", ".join(conflicts)]
        )
    return table.loc[:, [*key_columns, *columns]].drop_duplicates(list(key_columns), keep="first")


def _ordered_unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values))


def _single_cell_sample_metadata(
    table: pd.DataFrame,
    *,
    sample_columns: Sequence[str],
) -> pd.DataFrame:
    conflicts: list[str] = []
    for column in sample_columns:
        counts = table.groupby("sample_id", observed=True)[column].nunique(dropna=False)
        if (counts > 1).any():
            conflicts.append(column)
    if conflicts:
        raise DataValidationError(
            ["sample-level metadata conflict across records: " + ", ".join(conflicts)]
        )
    metadata = table.loc[:, ["sample_id", *sample_columns]].drop_duplicates(
        "sample_id", keep="first"
    )
    batch_ids = (
        table.groupby("sample_id", observed=True)["batch_id"]
        .agg(lambda values: "|".join(sorted(set(map(str, values)))))
        .rename("batch_ids")
        .reset_index()
    )
    return metadata.merge(batch_ids, on="sample_id", validate="one_to_one")


def _single_cell_qc(
    cell_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    product_compartments: set[str],
    product_timepoints: set[str],
    min_cells_per_group: int,
) -> pd.DataFrame:
    sample_columns = [
        "sample_id",
        "patient_id",
        "timepoint_id",
        "experimental_arm",
        "compartment",
        "batch_ids",
    ]
    expected_rows: list[dict[str, object]] = []
    for row in sample_metadata.loc[:, sample_columns].to_dict("records"):
        is_product = (
            str(row["compartment"]).lower() in product_compartments
            or str(row["timepoint_id"]) in product_timepoints
        )
        origins = ("CAR_T", "CAR_negative_product_T") if is_product else ("CAR_T", "endogenous_T")
        for origin in origins:
            expected_rows.append({**row, "t_cell_origin": origin})
    expected = pd.DataFrame.from_records(expected_rows)
    observed = (
        cell_metadata.groupby(["sample_id", "t_cell_origin"], observed=True, dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    qc = expected.merge(
        observed,
        on=["sample_id", "t_cell_origin"],
        how="left",
        validate="one_to_one",
    )
    qc["n_cells"] = qc["n_cells"].fillna(0).astype("int64")
    qc["minimum_cells"] = int(min_cells_per_group)
    qc["analysis_eligible"] = qc["n_cells"] >= min_cells_per_group
    qc["status"] = np.select(
        [qc["n_cells"] == 0, qc["n_cells"] < min_cells_per_group],
        ["absent", "below_minimum_cells"],
        default="eligible",
    )
    return qc


def _single_cell_missingness(
    sample_metadata: pd.DataFrame,
    *,
    timepoint_order: Sequence[str],
) -> pd.DataFrame:
    patient_arms = sample_metadata.loc[:, ["patient_id", "experimental_arm"]].drop_duplicates(
        ignore_index=True
    )
    expected = patient_arms.merge(
        pd.DataFrame({"timepoint_id": list(timepoint_order)}), how="cross"
    )
    observed = (
        sample_metadata.groupby(
            ["patient_id", "experimental_arm", "timepoint_id"],
            observed=True,
            dropna=False,
        )["sample_id"]
        .nunique()
        .rename("n_samples")
        .reset_index()
    )
    result = expected.merge(
        observed,
        on=["patient_id", "experimental_arm", "timepoint_id"],
        how="left",
        validate="one_to_one",
    )
    result["n_samples"] = result["n_samples"].fillna(0).astype("int64")
    result["status"] = np.where(result["n_samples"] > 0, "present", "missing")
    return result


def _empty_pseudobulk(group_columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=[*group_columns, "n_cells"])


def _aggregate_count_type(
    table: pd.DataFrame,
    *,
    feature_type: str,
    group_columns: Sequence[str],
    metadata_columns: Sequence[str],
    min_cells: int,
    long_format: bool,
    wide_feature_columns: Sequence[str] = (),
) -> pd.DataFrame:
    if long_format:
        subset = table.loc[table["feature_type"] == feature_type].copy()
        if subset.empty:
            return _empty_pseudobulk(group_columns)
        duplicate_key = ["sample_id", "cell_id", "feature_type", "feature_id"]
        if subset.duplicated(duplicate_key, keep=False).any():
            raise DataValidationError(
                [
                    "duplicate single-cell feature rows detected for key "
                    "sample_id/cell_id/feature_type/feature_id"
                ]
            )
        subset["feature_key"] = feature_type + "__" + subset["feature_id"].astype(str)
        wide = subset.pivot(
            index=["sample_id", "cell_id", *metadata_columns],
            columns="feature_key",
            values="count",
        ).fillna(0)
        wide = wide.reset_index()
        feature_columns = [
            column for column in wide.columns if str(column).startswith(feature_type + "__")
        ]
    else:
        feature_columns = list(wide_feature_columns)
        if not feature_columns:
            return _empty_pseudobulk(group_columns)
        duplicate_key = ["sample_id", "cell_id"]
        if table.duplicated(duplicate_key, keep=False).any():
            raise DataValidationError(
                ["wide single-cell input contains duplicate sample_id/cell_id rows"]
            )
        wide = table.loc[:, ["sample_id", "cell_id", *metadata_columns, *feature_columns]].copy()

    return single_cell.aggregate_pseudobulk(
        wide,
        feature_columns,
        group_columns=group_columns,
        min_cells=min_cells,
    )


def _stage_single_cell(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    manifest_path: Path,
    output_dir: Path,
) -> list[Path]:
    module_config = config.get("single_cell", {})
    modalities = tuple(module_config.get("modalities", ("scrna", "citeseq")))
    table = _load_tables(manifest, project_root, modalities=modalities)
    car_column = str(module_config.get("car_column", "car_positive"))
    origin_column = str(module_config.get("origin_column", "t_cell_origin"))
    state_column = str(module_config.get("state_column", "cell_state"))
    subtype_column = str(module_config.get("subtype_column", "t_cell_subtype"))
    configured_groups = tuple(module_config.get("group_columns", ()))
    additional_groups = [
        str(column)
        for column in configured_groups
        if str(column)
        not in {
            "sample_id",
            "patient_id",
            "timepoint_id",
            origin_column,
            "t_cell_origin",
        }
    ]
    required = (
        "cell_id",
        "sample_id",
        "patient_id",
        "timepoint_id",
        "experimental_arm",
        "compartment",
        "batch_id",
        origin_column,
        car_column,
        state_column,
        "annotation_version",
        "car_gate_version",
        *additional_groups,
    )
    single_cell.validate_required_columns(table, required, table_name="single_cell")
    table["car_positive"] = table[car_column].map(_as_bool)
    table["t_cell_origin"] = table[origin_column].astype(str)
    table["cell_state"] = table[state_column].astype(str)
    if subtype_column in table.columns:
        subtype_values = table[subtype_column].astype(str).str.strip()
        table["t_cell_subtype"] = subtype_values.mask(subtype_values == "", "unspecified")
    else:
        table["t_cell_subtype"] = "unspecified"

    product_compartments = {
        str(value).lower() for value in module_config.get("product_compartments", ("product",))
    }
    product_timepoints = {
        str(value) for value in module_config.get("product_timepoints", ("PRODUCT",))
    }
    product_mask = table["compartment"].astype(str).str.lower().isin(product_compartments) | table[
        "timepoint_id"
    ].astype(str).isin(product_timepoints)
    expected_origin = np.select(
        [table["car_positive"], product_mask],
        ["CAR_T", "CAR_negative_product_T"],
        default="endogenous_T",
    )
    if not np.array_equal(table["t_cell_origin"].to_numpy(), expected_origin):
        mismatch = table.loc[
            table["t_cell_origin"].to_numpy() != expected_origin,
            ["sample_id", "cell_id", "timepoint_id", "compartment", "t_cell_origin"],
        ].head(5)
        raise DataValidationError(
            [
                "t_cell_origin is inconsistent with CAR gate and product/blood context; "
                f"examples: {mismatch.to_dict('records')}"
            ]
        )

    metadata_columns = (
        "patient_id",
        "timepoint_id",
        "experimental_arm",
        "compartment",
        "t_cell_origin",
        "car_positive",
        "cell_state",
        "t_cell_subtype",
        "annotation_version",
        "car_gate_version",
        *additional_groups,
    )
    cell_metadata = _consistent_cell_metadata(table, metadata_columns)
    sample_metadata = _single_cell_sample_metadata(
        table,
        sample_columns=(
            "patient_id",
            "timepoint_id",
            "experimental_arm",
            "compartment",
        ),
    )
    configured_timepoints = tuple(module_config.get("timepoint_order", ()))
    timepoint_order = configured_timepoints or tuple(
        sorted(sample_metadata["timepoint_id"].unique())
    )
    if len(timepoint_order) != len(set(timepoint_order)):
        raise DataValidationError(["single_cell.timepoint_order contains duplicates"])
    missing_timepoints = set(sample_metadata["timepoint_id"]).difference(timepoint_order)
    if missing_timepoints:
        raise DataValidationError(
            [
                "single_cell.timepoint_order omits observed timepoints: "
                + ", ".join(sorted(missing_timepoints))
            ]
        )

    configured_states = tuple(module_config.get("states", ()))
    states = configured_states or tuple(sorted(cell_metadata["cell_state"].unique()))
    other_state = str(module_config.get("other_state_label", "other"))
    cell_metadata["analysis_state"] = np.where(
        cell_metadata["cell_state"].isin(states),
        cell_metadata["cell_state"],
        other_state,
    )
    table["analysis_state"] = np.where(
        table["cell_state"].isin(states), table["cell_state"], other_state
    )
    min_cells = int(module_config.get("minimum_cells_per_group", 3))
    min_pseudobulk_cells = int(module_config.get("minimum_cells_per_pseudobulk", 3))
    min_patients = int(module_config.get("minimum_patients_per_stratum", 3))
    fraction_groups = _ordered_unique(
        (
            "sample_id",
            "patient_id",
            "timepoint_id",
            "experimental_arm",
            "compartment",
            *additional_groups,
            "t_cell_origin",
            "t_cell_subtype",
        )
    )
    single_cell.validate_required_columns(
        cell_metadata, fraction_groups, table_name="single_cell metadata"
    )
    fractions = single_cell.calculate_state_fractions(
        cell_metadata,
        group_columns=fraction_groups,
        state_col="analysis_state",
        states=states,
        denominator="all",
        other_label=other_state,
        min_cells_per_group=min_cells,
    )
    fractions = fractions.rename(columns={"analysis_state": "cell_state"})
    fractions = fractions.merge(
        sample_metadata.loc[:, ["sample_id", "batch_ids"]],
        on="sample_id",
        validate="many_to_one",
    )

    inference_groups = _ordered_unique(
        (
            "timepoint_id",
            "experimental_arm",
            "compartment",
            *additional_groups,
            "t_cell_origin",
            "t_cell_subtype",
            "cell_state",
        )
    )
    eligible_fractions = fractions.loc[fractions["passes_min_cells"]].copy()
    duplicate_inference = eligible_fractions.duplicated(
        [*inference_groups, "patient_id"], keep=False
    )
    if duplicate_inference.any():
        raise DataValidationError(
            [
                "multiple samples remain for one patient within a single-cell inference "
                "stratum; define and document a biological-replicate aggregation policy"
            ]
        )
    confidence_level = float(config.get("statistics", {}).get("confidence_level", 0.95))
    bootstrap_iterations = int(module_config.get("bootstrap_iterations", 2000))
    bootstrap_columns = [
        *inference_groups,
        "n_patients",
        "estimate",
        "ci_low",
        "ci_high",
        "confidence_level",
        "n_boot",
        "ci_method",
        "status",
    ]
    if eligible_fractions.empty:
        bootstrap = pd.DataFrame(columns=bootstrap_columns)
    else:
        bootstrap = single_cell.bootstrap_patient_ci(
            eligible_fractions,
            "fraction",
            patient_col="patient_id",
            group_columns=inference_groups,
            confidence_level=confidence_level,
            n_boot=bootstrap_iterations,
            random_state=_effective_seed(config),
            min_patients=min_patients,
            on_insufficient="record",
        )
    all_inference_strata = fractions.loc[:, list(inference_groups)].drop_duplicates(
        ignore_index=True
    )
    bootstrap = all_inference_strata.merge(
        bootstrap,
        on=list(inference_groups),
        how="left",
        validate="one_to_one",
    )
    no_eligible_patients = bootstrap["status"].isna()
    bootstrap.loc[no_eligible_patients, "n_patients"] = 0
    bootstrap.loc[no_eligible_patients, "confidence_level"] = confidence_level
    bootstrap.loc[no_eligible_patients, "n_boot"] = bootstrap_iterations
    bootstrap.loc[no_eligible_patients, "ci_method"] = "patient_percentile_bootstrap"
    bootstrap.loc[no_eligible_patients, "status"] = "no_eligible_patients"
    bootstrap["n_patients"] = bootstrap["n_patients"].astype("int64")
    bootstrap["n_boot"] = bootstrap["n_boot"].astype("int64")
    bootstrap.insert(0, "estimand", "mean_patient_fraction")

    reference_timepoint = str(module_config.get("reference_timepoint", ""))
    paired_origins = {
        str(value) for value in module_config.get("paired_change_origins", ("CAR_T",))
    }
    paired_input = eligible_fractions.loc[
        eligible_fractions["t_cell_origin"].isin(paired_origins)
    ].copy()
    paired_groups = _ordered_unique(
        (
            "experimental_arm",
            *additional_groups,
            "t_cell_origin",
            "t_cell_subtype",
            "cell_state",
        )
    )
    if (
        reference_timepoint
        and reference_timepoint in set(paired_input["timepoint_id"])
        and paired_input["timepoint_id"].nunique() > 1
    ):
        paired_comparisons = tuple(
            timepoint
            for timepoint in timepoint_order
            if timepoint != reference_timepoint and timepoint in set(paired_input["timepoint_id"])
        )
        paired_changes = single_cell.bootstrap_paired_change_ci(
            paired_input,
            "fraction",
            reference_timepoint=reference_timepoint,
            timepoint_col="timepoint_id",
            patient_col="patient_id",
            group_columns=paired_groups,
            comparison_timepoints=paired_comparisons,
            confidence_level=confidence_level,
            n_boot=bootstrap_iterations,
            random_state=_effective_seed(config),
            min_patients=min_patients,
            on_insufficient="record",
        )
        paired_changes.insert(0, "estimand", "mean_within_patient_change")
    else:
        paired_changes = pd.DataFrame(
            columns=[
                "estimand",
                *paired_groups,
                "reference_timepoint",
                "comparison_timepoint",
                "n_patients",
                "estimate",
                "ci_low",
                "ci_high",
                "confidence_level",
                "n_boot",
                "ci_method",
                "status",
            ]
        )

    qc = _single_cell_qc(
        cell_metadata,
        sample_metadata,
        product_compartments=product_compartments,
        product_timepoints=product_timepoints,
        min_cells_per_group=min_cells,
    )
    missingness = _single_cell_missingness(
        sample_metadata,
        timepoint_order=timepoint_order,
    )

    count_metadata_columns = _ordered_unique(
        (
            "patient_id",
            "timepoint_id",
            "experimental_arm",
            "compartment",
            "batch_id",
            *additional_groups,
            "t_cell_origin",
            "t_cell_subtype",
            "analysis_state",
        )
    )
    pseudobulk_groups = _ordered_unique(
        (
            "sample_id",
            "patient_id",
            "timepoint_id",
            "experimental_arm",
            "compartment",
            "batch_id",
            *additional_groups,
            "t_cell_origin",
            "t_cell_subtype",
            "analysis_state",
        )
    )
    long_fields = {"feature_id", "feature_type", "count"}
    present_long_fields = long_fields.intersection(table.columns)
    if present_long_fields and present_long_fields != long_fields:
        raise DataValidationError(
            ["long single-cell counts require feature_id, feature_type, and count together"]
        )
    if present_long_fields == long_fields:
        if (table["feature_id"].astype(str).str.len() == 0).any():
            raise DataValidationError(["single-cell feature_id may not be empty"])
        feature_types = set(table["feature_type"].astype(str))
        unsupported_types = feature_types.difference({"RNA", "ADT"})
        if unsupported_types:
            raise DataValidationError(
                ["unsupported feature_type values: " + ", ".join(sorted(unsupported_types))]
            )
        table["count"] = pd.to_numeric(table["count"], errors="raise")
        counts = table["count"].to_numpy(dtype=float)
        if not np.isfinite(counts).all() or (counts < 0).any():
            raise DataValidationError(["single-cell raw counts must be finite and non-negative"])
        if not np.allclose(counts, np.round(counts)):
            raise DataValidationError(
                ["single-cell pseudobulk requires untransformed integer counts"]
            )
        rna_pseudobulk = _aggregate_count_type(
            table,
            feature_type="RNA",
            group_columns=pseudobulk_groups,
            metadata_columns=count_metadata_columns,
            min_cells=min_pseudobulk_cells,
            long_format=True,
        )
        adt_pseudobulk = _aggregate_count_type(
            table,
            feature_type="ADT",
            group_columns=pseudobulk_groups,
            metadata_columns=count_metadata_columns,
            min_cells=min_pseudobulk_cells,
            long_format=True,
        )
    else:
        rna_suffix = str(module_config.get("rna_count_suffix", "_RNA_count"))
        adt_suffix = str(module_config.get("adt_count_suffix", "_ADT_count"))
        rna_columns = sorted(column for column in table.columns if str(column).endswith(rna_suffix))
        adt_columns = sorted(column for column in table.columns if str(column).endswith(adt_suffix))
        ambiguous_counts = sorted(
            column
            for column in table.columns
            if str(column).endswith("_count") and column not in {*rna_columns, *adt_columns}
        )
        if ambiguous_counts:
            raise DataValidationError(
                ["count columns lack an explicit RNA or ADT suffix: " + ", ".join(ambiguous_counts)]
            )
        feature_columns = [*rna_columns, *adt_columns]
        if feature_columns:
            for column in feature_columns:
                table[column] = pd.to_numeric(table[column], errors="raise")
            values = table[feature_columns].to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < 0).any():
                raise DataValidationError(
                    ["single-cell raw counts must be finite and non-negative"]
                )
            if not np.allclose(values, np.round(values)):
                raise DataValidationError(
                    ["single-cell pseudobulk requires untransformed integer counts"]
                )
            rna_pseudobulk = _aggregate_count_type(
                table,
                feature_type="RNA",
                group_columns=pseudobulk_groups,
                metadata_columns=count_metadata_columns,
                min_cells=min_pseudobulk_cells,
                long_format=False,
                wide_feature_columns=rna_columns,
            )
            adt_pseudobulk = _aggregate_count_type(
                table,
                feature_type="ADT",
                group_columns=pseudobulk_groups,
                metadata_columns=count_metadata_columns,
                min_cells=min_pseudobulk_cells,
                long_format=False,
                wide_feature_columns=adt_columns,
            )
        else:
            rna_pseudobulk = _empty_pseudobulk(pseudobulk_groups)
            adt_pseudobulk = _empty_pseudobulk(pseudobulk_groups)
    rna_pseudobulk = rna_pseudobulk.rename(columns={"analysis_state": "cell_state"})
    adt_pseudobulk = adt_pseudobulk.rename(columns={"analysis_state": "cell_state"})

    state_order = (*states, other_state)
    for frame in (fractions, bootstrap, rna_pseudobulk, adt_pseudobulk):
        if "timepoint_id" in frame.columns:
            frame["timepoint_id"] = pd.Categorical(
                frame["timepoint_id"],
                categories=timepoint_order,
                ordered=True,
            )
        if "cell_state" in frame.columns:
            frame["cell_state"] = pd.Categorical(
                frame["cell_state"],
                categories=state_order,
                ordered=True,
            )
    for frame in (qc, missingness):
        frame["timepoint_id"] = pd.Categorical(
            frame["timepoint_id"],
            categories=timepoint_order,
            ordered=True,
        )
    if not paired_changes.empty:
        for column in ("reference_timepoint", "comparison_timepoint"):
            paired_changes[column] = pd.Categorical(
                paired_changes[column],
                categories=timepoint_order,
                ordered=True,
            )
        paired_changes["cell_state"] = pd.Categorical(
            paired_changes["cell_state"],
            categories=state_order,
            ordered=True,
        )

    table_dir = output_dir / "tables"
    outputs = [
        _write_table(
            fractions,
            table_dir / "single_cell_state_fractions.tsv",
            sort_by=("patient_id", "timepoint_id", "t_cell_origin", "cell_state"),
        ),
        _write_table(
            bootstrap,
            table_dir / "single_cell_state_bootstrap.tsv",
            sort_by=("timepoint_id", "experimental_arm", "t_cell_origin", "cell_state"),
        ),
        _write_table(
            paired_changes,
            table_dir / "single_cell_state_paired_changes.tsv",
            sort_by=(
                "experimental_arm",
                "t_cell_origin",
                "cell_state",
                "comparison_timepoint",
            ),
        ),
        _write_table(
            rna_pseudobulk,
            table_dir / "single_cell_rna_pseudobulk.tsv",
            sort_by=("patient_id", "timepoint_id", "t_cell_origin", "cell_state"),
        ),
        _write_table(
            adt_pseudobulk,
            table_dir / "single_cell_adt_pseudobulk.tsv",
            sort_by=("patient_id", "timepoint_id", "t_cell_origin", "cell_state"),
        ),
        _write_table(
            qc,
            table_dir / "single_cell_sample_qc.tsv",
            sort_by=("patient_id", "timepoint_id", "sample_id", "t_cell_origin"),
        ),
        _write_table(
            missingness,
            table_dir / "single_cell_missingness.tsv",
            sort_by=("patient_id", "experimental_arm", "timepoint_id"),
        ),
    ]
    input_paths = [path for _, path in _source_paths(manifest, project_root, modalities=modalities)]
    write_stage_marker(
        output_dir,
        stage="single-cell",
        outputs=outputs,
        inputs=[config["_config_path"], manifest_path, *input_paths],
        metadata={
            "analysis_unit": "patient_id",
            "n_cells": int(len(cell_metadata)),
            "states": list(states),
            "other_state_label": other_state,
            "minimum_cells_per_group": min_cells,
            "minimum_patients_per_stratum": min_patients,
            "random_seed": _effective_seed(config),
            "timepoint_order": list(timepoint_order),
            "annotation_versions": sorted(set(table["annotation_version"].astype(str))),
            "car_gate_versions": sorted(set(table["car_gate_version"].astype(str))),
            "synthetic_only": _is_synthetic_only(manifest),
        },
    )
    return outputs


def _stage_spatial(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    manifest_path: Path,
    output_dir: Path,
) -> list[Path]:
    table = _load_tables(manifest, project_root, modalities=("spatial",))
    table["x"] = pd.to_numeric(table["x"], errors="raise")
    table["y"] = pd.to_numeric(table["y"], errors="raise")
    for column in ("is_ligand_source", "is_car_t"):
        table[column] = table[column].map(_as_bool)
    validated = spatial.validate_spatial_table(
        table,
        boolean_columns=("is_ligand_source", "is_car_t"),
    )
    module_config = config.get("spatial", {})
    distances = spatial.distance_to_landmarks(
        validated,
        query_columns={
            "ligand_source": str(module_config.get("ligand_source_column", "is_ligand_source")),
            "car_t": str(module_config.get("car_t_column", "is_car_t")),
        },
        structure_col=str(module_config.get("landmark_column", "structure")),
        coordinate_cols=tuple(module_config.get("coordinate_columns", ("x", "y"))),
        coordinate_scale_um=float(module_config.get("coordinate_scale_um", 1.0)),
    )
    patient_summary = spatial.summarize_patient_distances(distances)

    observed_types = sorted(validated["cell_type"].astype(str).unique())
    focal_types = [name for name in ("CAR_T", "endogenous_T") if name in observed_types]
    if not focal_types:
        focal_types = [name for name in observed_types if "T" in name][:1]
    neighbor_types = [name for name in observed_types if name not in set(focal_types)]
    if not focal_types or not neighbor_types:
        raise DataValidationError(
            ["spatial neighborhood analysis requires at least one T-cell and one neighbor type"]
        )
    neighborhood = spatial.neighborhood_enrichment(
        validated,
        radius=float(module_config.get("neighborhood_radius_um", 50.0)),
        focal_types=focal_types,
        neighbor_types=neighbor_types,
        coordinate_scale_um=float(module_config.get("coordinate_scale_um", 1.0)),
        n_permutations=int(module_config.get("permutation_iterations", 999)),
        seed=_effective_seed(config),
    )

    table_dir = output_dir / "tables"
    outputs = [
        _write_table(
            distances,
            table_dir / "spatial_distances.tsv",
            sort_by=("patient_id", "section_id", "query_role", "cell_id", "landmark"),
        ),
        _write_table(
            patient_summary,
            table_dir / "spatial_patient_summary.tsv",
            sort_by=("patient_id", "query_role", "landmark"),
        ),
        _write_table(
            neighborhood,
            table_dir / "spatial_neighborhood_enrichment.tsv",
            sort_by=("focal_type", "neighbor_type"),
        ),
    ]
    inputs = [path for _, path in _source_paths(manifest, project_root, modalities=("spatial",))]
    write_stage_marker(
        output_dir,
        stage="spatial",
        outputs=outputs,
        inputs=[config["_config_path"], manifest_path, *inputs],
        metadata={
            "analysis_unit": "patient_id",
            "permutation_strata": ["patient_id", "section_id"],
            "random_seed": _effective_seed(config),
            "synthetic_only": _is_synthetic_only(manifest),
        },
    )
    return outputs


def _protein_summary(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = table.copy()
    table["response"] = pd.to_numeric(table["response"], errors="raise")
    biological_keys = [
        "biological_replicate",
        "patient_id",
        "timepoint_id",
        "analyte",
        "molecular_form",
        "unit",
        "assay",
    ]
    patient_level = (
        table.groupby(biological_keys, observed=True, dropna=False)["response"]
        .mean()
        .rename("response")
        .reset_index()
    )
    cohort_keys = ["timepoint_id", "analyte", "molecular_form", "unit", "assay"]
    cohort = (
        patient_level.groupby(cohort_keys, observed=True, dropna=False)["response"]
        .agg(n_biological_replicates="count", mean="mean", sd="std", median="median")
        .reset_index()
    )
    return patient_level, cohort


def _functional_design_columns(table: pd.DataFrame) -> list[str]:
    candidates = (
        "patient_id",
        "timepoint_id",
        "assay",
        "condition",
        "product",
        "unit",
        "chemokine",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "cxcl16_form",
        "pair_id",
    )
    return [column for column in candidates if column in table.columns]


def _functional_comparisons(
    table: pd.DataFrame,
    *,
    seed: int,
    iterations: int,
    confidence_level: float,
) -> pd.DataFrame:
    context_columns = [
        column
        for column in (
            "assay",
            "condition",
            "unit",
            "chemokine",
            "dose",
            "dose_unit",
            "time",
            "time_unit",
            "cxcl16_form",
            "timepoint_id",
        )
        if column in table.columns
    ]
    records: list[dict[str, object]] = []
    iterator = table.groupby(context_columns, observed=True, dropna=False, sort=True)
    for context_key, subset in iterator:
        if not isinstance(context_key, tuple):
            context_key = (context_key,)
        if set(subset["product"].dropna()) != {"identical_control", "receptor_matched"}:
            continue
        pair_column = "pair_id" if "pair_id" in subset.columns else "biological_replicate"
        result = functional.compare_receptor_matched_product(
            subset,
            paired=True,
            pair_col=pair_column,
            biological_replicate_col="biological_replicate",
            response_col="response",
            technical_replicate_col="technical_replicate",
            design_cols=["product", pair_column],
            bootstrap_iterations=iterations,
            confidence_level=confidence_level,
            seed=seed,
        )
        record = dict(zip(context_columns, context_key, strict=True))
        record.update(asdict(result))
        record["notes"] = "; ".join(result.notes)
        records.append(record)
    comparisons = pd.DataFrame.from_records(records)
    if not comparisons.empty:
        comparisons["q_value_bh"] = _bh_adjust(comparisons["p_value"].to_numpy())
    return comparisons


def _stage_functional(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    manifest_path: Path,
    output_dir: Path,
) -> list[Path]:
    functional_table = _load_tables(manifest, project_root, modalities=("functional",))
    for optional_column in (
        "technical_replicate",
        "pair_id",
        "chemokine",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "cxcl16_form",
    ):
        if optional_column in functional_table.columns:
            functional_table[optional_column] = functional_table[optional_column].replace("", pd.NA)
    functional_table["response"] = pd.to_numeric(functional_table["response"], errors="raise")
    design = _functional_design_columns(functional_table)
    qc = functional.validate_functional_data(
        functional_table,
        biological_replicate_col="biological_replicate",
        response_col="response",
        assay_col="assay",
        technical_replicate_col="technical_replicate",
        design_cols=design,
    )
    if not qc.passed:
        raise DataValidationError(list(qc.errors))
    biological = functional.aggregate_biological_replicates(
        functional_table,
        biological_replicate_col="biological_replicate",
        response_col="response",
        technical_replicate_col="technical_replicate",
        design_cols=design,
        method="mean",
    )
    module_config = config.get("functional", {})
    comparisons = _functional_comparisons(
        functional_table,
        seed=_effective_seed(config),
        iterations=int(module_config.get("bootstrap_iterations", 2000)),
        confidence_level=float(config.get("statistics", {}).get("confidence_level", 0.95)),
    )

    protein_table = _load_tables(manifest, project_root, modalities=("protein",))
    protein_table["response"] = pd.to_numeric(protein_table["response"], errors="raise")
    protein_qc = functional.validate_protein_data(
        protein_table,
        biological_replicate_col="biological_replicate",
        response_col="response",
        analyte_col="analyte",
        form_col="molecular_form",
        unit_col="unit",
        assay_col="assay",
        technical_replicate_col="technical_replicate",
    )
    if not protein_qc.passed:
        raise DataValidationError(list(protein_qc.errors))
    protein_patient, protein_cohort = _protein_summary(protein_table)
    table_dir = output_dir / "tables"
    outputs = [
        _write_table(
            biological,
            table_dir / "functional_biological_summary.tsv",
            sort_by=("assay", "condition", "biological_replicate", "product"),
        ),
        _write_table(
            comparisons,
            table_dir / "functional_product_comparisons.tsv",
            sort_by=("assay", "condition"),
        ),
        _write_table(
            protein_patient,
            table_dir / "cxcl16_patient_measurements.tsv",
            sort_by=("analyte", "molecular_form", "patient_id", "timepoint_id"),
        ),
        _write_table(
            protein_cohort,
            table_dir / "cxcl16_cohort_summary.tsv",
            sort_by=("analyte", "molecular_form", "timepoint_id"),
        ),
        write_json(table_dir / "functional_qc.json", qc.to_dict()),
        write_json(table_dir / "protein_qc.json", protein_qc.to_dict()),
    ]
    inputs = [
        *[path for _, path in _source_paths(manifest, project_root, modalities=("functional",))],
        *[path for _, path in _source_paths(manifest, project_root, modalities=("protein",))],
    ]
    write_stage_marker(
        output_dir,
        stage="functional",
        outputs=outputs,
        inputs=[config["_config_path"], manifest_path, *inputs],
        metadata={
            "analysis_unit": "biological_replicate",
            "paired_comparator": "identical_control",
            "modified_product": "receptor_matched",
            "cxcl16_cross_scale_test": False,
            "random_seed": _effective_seed(config),
            "synthetic_only": _is_synthetic_only(manifest),
        },
    )
    return outputs


def _mark_synthetic(figure: plt.Figure, synthetic: bool) -> None:
    if synthetic:
        figure.text(
            0.5,
            0.008,
            "SYNTHETIC SOFTWARE TEST — NOT SCIENTIFIC EVIDENCE",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#A61B1B",
            weight="bold",
        )


def _save_figure(
    figure: plt.Figure,
    base_path: Path,
    *,
    formats: Sequence[str],
    dpi: int,
    title: str,
) -> list[Path]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    fixed_timestamp = datetime(2026, 9, 1, tzinfo=UTC)
    for extension in formats:
        path = base_path.with_suffix(f".{extension}")
        if extension.lower() == "pdf":
            metadata = {
                "Title": title,
                "Subject": "Chemokine-CAR-T reproducibility analysis",
                "CreationDate": fixed_timestamp,
                "ModDate": fixed_timestamp,
            }
        elif extension.lower() == "svg":
            metadata = {
                "Title": title,
                "Description": "Chemokine-CAR-T reproducibility analysis",
                "Date": "2026-09-01",
            }
        else:
            metadata = {
                "Title": title,
                "Description": "Chemokine-CAR-T reproducibility analysis",
            }
        buffer = io.BytesIO()
        figure.savefig(
            buffer,
            format=extension,
            dpi=dpi if extension.lower() == "png" else None,
            bbox_inches="tight",
            metadata=metadata,
        )
        payload = buffer.getvalue()
        if not payload:
            raise RuntimeError(f"Figure writer produced an empty file: {path}")
        path.write_bytes(payload)
        if path.stat().st_size != len(payload):
            raise RuntimeError(f"Figure writer produced an incomplete file: {path}")
        outputs.append(path)
    plt.close(figure)
    return outputs


def _figure_single_cell(config: Mapping[str, Any], output_dir: Path, synthetic: bool) -> list[Path]:
    fractions = pd.read_csv(output_dir / "tables" / "single_cell_state_fractions.tsv", sep="\t")
    if "passes_min_cells" in fractions.columns:
        fractions = fractions.loc[fractions["passes_min_cells"].map(_as_bool)].copy()
    if fractions.empty:
        raise DataValidationError(
            ["no single-cell state fractions pass the configured minimum cell threshold"]
        )
    module_config = config.get("single_cell", {})
    states = tuple(module_config.get("states", ())) or tuple(
        sorted(fractions["cell_state"].unique())
    )
    configured_timepoints = tuple(module_config.get("timepoint_order", ()))
    timepoint_order = configured_timepoints or None
    fractions["plot_series"] = (
        fractions["experimental_arm"].astype(str)
        + " | "
        + fractions["t_cell_origin"].astype(str)
        + " | "
        + fractions["t_cell_subtype"].astype(str)
    )
    figure, _ = single_cell.plot_state_fractions(
        fractions,
        states=states,
        origin_col="plot_series",
        timepoint_order=timepoint_order,
    )
    figure.suptitle("Patient-level CAR-T and endogenous T-cell states", y=1.02)
    _mark_synthetic(figure, synthetic)
    return _save_figure(
        figure,
        output_dir / "figures" / "figure6a_single_cell_states",
        formats=tuple(config.get("output", {}).get("figure_formats", ("pdf", "svg", "png"))),
        dpi=int(config.get("output", {}).get("png_dpi", 300)),
        title="Patient-level CAR-T and endogenous T-cell states",
    )


def _figure_spatial(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    output_dir: Path,
    synthetic: bool,
) -> list[Path]:
    table = _load_tables(manifest, project_root, modalities=("spatial",))
    table["x"] = pd.to_numeric(table["x"], errors="raise")
    table["y"] = pd.to_numeric(table["y"], errors="raise")
    for column in ("is_ligand_source", "is_car_t"):
        table[column] = table[column].map(_as_bool)
    figure, _ = spatial.plot_spatial_panel(table)
    figure.suptitle("Chemokine sources, CAR-T cells, and tissue landmarks", y=1.02)
    _mark_synthetic(figure, synthetic)
    return _save_figure(
        figure,
        output_dir / "figures" / "figure6b_spatial_map",
        formats=tuple(config.get("output", {}).get("figure_formats", ("pdf", "svg", "png"))),
        dpi=int(config.get("output", {}).get("png_dpi", 300)),
        title="Chemokine sources, CAR-T cells, and tissue landmarks",
    )


def _figure_cxcl16(config: Mapping[str, Any], output_dir: Path, synthetic: bool) -> list[Path]:
    data = pd.read_csv(output_dir / "tables" / "cxcl16_patient_measurements.tsv", sep="\t")
    data = data.loc[data["analyte"].eq("CXCL16")].copy()
    if data.empty:
        raise DataValidationError(["protein table contains no CXCL16 measurements"])
    panels = list(data.groupby(["molecular_form", "unit", "assay"], observed=True, sort=True))
    figure, axes = plt.subplots(
        1,
        len(panels),
        figsize=(5.2 * len(panels), 4.2),
        squeeze=False,
        constrained_layout=True,
    )
    for axis, ((form, unit, assay), subset) in zip(axes.ravel(), panels, strict=True):
        timepoints = list(pd.unique(subset["timepoint_id"]))
        x = {timepoint: index for index, timepoint in enumerate(timepoints)}
        for _, patient in subset.groupby("biological_replicate", observed=True, sort=True):
            patient = patient.sort_values(
                "timepoint_id",
                key=lambda values: values.map(x),
                kind="mergesort",
            )
            axis.plot(
                patient["timepoint_id"].map(x),
                patient["response"],
                marker="o",
                linewidth=1.0,
                alpha=0.65,
            )
        axis.set_xticks(range(len(timepoints)), timepoints)
        axis.set_title(f"{form} CXCL16 | {assay}")
        axis.set_xlabel("Time point")
        axis.set_ylabel(str(unit))
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle("CXCL16 protein forms shown on assay-specific scales", y=1.04)
    _mark_synthetic(figure, synthetic)
    return _save_figure(
        figure,
        output_dir / "figures" / "figure6c_cxcl16_protein_validation",
        formats=tuple(config.get("output", {}).get("figure_formats", ("pdf", "svg", "png"))),
        dpi=int(config.get("output", {}).get("png_dpi", 300)),
        title="CXCL16 protein validation",
    )


def _figure_functional(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    output_dir: Path,
    synthetic: bool,
) -> list[Path]:
    data = _load_tables(manifest, project_root, modalities=("functional",))
    data["response"] = pd.to_numeric(data["response"], errors="raise")
    assays = [name for name in ("chemotaxis", "retention", "egress") if name in set(data["assay"])]
    if not assays:
        raise DataValidationError(["functional table has no chemotaxis, retention, or egress rows"])
    figure, axes = plt.subplots(
        1,
        len(assays),
        figsize=(5.0 * len(assays), 4.2),
        squeeze=False,
        constrained_layout=True,
    )
    for axis, assay in zip(axes.ravel(), assays, strict=True):
        subset = data.loc[data["assay"].eq(assay)].copy()
        units = sorted(subset["unit"].unique())
        if len(units) != 1:
            raise DataValidationError(
                [f"assay {assay!r} contains multiple units; split it into prespecified endpoints"]
            )
        functional.plot_migration(
            subset,
            condition_col="condition",
            group_col="product",
            biological_replicate_col="biological_replicate",
            response_col="response",
            technical_replicate_col="technical_replicate",
            response_label=str(units[0]),
            ax=axis,
        )
        axis.set_title(assay.capitalize())
    figure.suptitle("Chemotaxis, retention, and egress assays", y=1.04)
    _mark_synthetic(figure, synthetic)
    return _save_figure(
        figure,
        output_dir / "figures" / "figure6d_functional_assays",
        formats=tuple(config.get("output", {}).get("figure_formats", ("pdf", "svg", "png"))),
        dpi=int(config.get("output", {}).get("png_dpi", 300)),
        title="Chemotaxis, retention, and egress assays",
    )


def _figure_effects(config: Mapping[str, Any], output_dir: Path, synthetic: bool) -> list[Path]:
    comparisons = pd.read_csv(
        output_dir / "tables" / "functional_product_comparisons.tsv", sep="\t"
    )
    if comparisons.empty:
        raise DataValidationError(["no complete modified-versus-control comparisons are available"])
    labels = comparisons["assay"].astype(str)
    if "condition" in comparisons.columns:
        labels = labels + " | " + comparisons["condition"].astype(str)
    y = np.arange(len(comparisons))
    estimate = comparisons["mean_difference"].to_numpy(dtype=float)
    low = comparisons["mean_difference_ci_low"].to_numpy(dtype=float)
    high = comparisons["mean_difference_ci_high"].to_numpy(dtype=float)
    figure, axis = plt.subplots(
        figsize=(7.2, max(3.4, 0.7 * len(comparisons) + 1.8)), constrained_layout=True
    )
    axis.errorbar(
        estimate,
        y,
        xerr=np.vstack((estimate - low, high - estimate)),
        fmt="o",
        color="#24557A",
        ecolor="#7EA6C2",
        capsize=3,
    )
    axis.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    axis.set_yticks(y, labels)
    axis.set_xlabel("Mean difference: receptor-matched minus identical control")
    axis.set_title("Paired biological-unit effect estimates")
    axis.spines[["top", "right"]].set_visible(False)
    _mark_synthetic(figure, synthetic)
    return _save_figure(
        figure,
        output_dir / "figures" / "figure6e_matched_product_effects",
        formats=tuple(config.get("output", {}).get("figure_formats", ("pdf", "svg", "png"))),
        dpi=int(config.get("output", {}).get("png_dpi", 300)),
        title="Matched product effect estimates",
    )


def _stage_figures(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    manifest_path: Path,
    output_dir: Path,
) -> list[Path]:
    synthetic = _is_synthetic_only(manifest)
    for upstream_stage in ("single-cell", "spatial", "functional"):
        _require_current_stage_marker(
            output_dir,
            stage=upstream_stage,
            config_path=config["_config_path"],
            manifest_path=manifest_path,
            synthetic_only=synthetic,
            random_seed=_effective_seed(config),
        )
    outputs = [
        *_figure_single_cell(config, output_dir, synthetic),
        *_figure_spatial(config, project_root, manifest, output_dir, synthetic),
        *_figure_cxcl16(config, output_dir, synthetic),
        *_figure_functional(config, project_root, manifest, output_dir, synthetic),
        *_figure_effects(config, output_dir, synthetic),
    ]
    manifest_rows = [file_record(path, root=output_dir) for path in outputs]
    figure_manifest = _write_table(
        pd.DataFrame(manifest_rows),
        output_dir / "figures" / "figure_manifest.tsv",
        sort_by=("path",),
    )
    outputs.append(figure_manifest)
    write_stage_marker(
        output_dir,
        stage="figures",
        outputs=outputs,
        inputs=[
            config["_config_path"],
            manifest_path,
            output_dir / "tables" / "single_cell_state_fractions.tsv",
            output_dir / "tables" / "spatial_patient_summary.tsv",
            output_dir / "tables" / "cxcl16_patient_measurements.tsv",
            output_dir / "tables" / "functional_product_comparisons.tsv",
        ],
        metadata={
            "figure_count": 5,
            "random_seed": _effective_seed(config),
            "synthetic_only": synthetic,
        },
    )
    return outputs


def _stage_report(
    config: Mapping[str, Any],
    project_root: Path,
    manifest: pd.DataFrame,
    manifest_path: Path,
    output_dir: Path,
) -> list[Path]:
    seed = _effective_seed(config)
    synthetic = _is_synthetic_only(manifest)
    for upstream_stage in ("single-cell", "spatial", "functional", "figures"):
        _require_current_stage_marker(
            output_dir,
            stage=upstream_stage,
            config_path=config["_config_path"],
            manifest_path=manifest_path,
            synthetic_only=synthetic,
            random_seed=seed,
        )
    input_paths = [
        path
        for _, path in _source_paths(
            manifest,
            project_root,
            modalities=tuple(manifest["modality"].unique()),
        )
    ]
    if synthetic:
        conclusion = (
            "This run used deterministic synthetic inputs. It verifies software execution only "
            "and provides no biological or clinical evidence."
        )
    else:
        conclusion = (
            "This report records the completed analytical workflow. Scientific interpretation "
            "requires confirmation that every included record is evidence eligible and that all "
            "prespecified quality-control criteria were met."
        )
    report_text = "\n".join(
        [
            "# Reproducibility report",
            "",
            conclusion,
            "",
            "## Completed components",
            "",
            "- patient-level CAR-T and endogenous T-cell state summaries;",
            "- spatial distances and conditional patient/section-blocked neighborhood analysis "
            "(not a population-level patient test; not yet tissue-mask constrained);",
            "- assay-specific soluble and membrane CXCL16 summaries;",
            "- chemotaxis, retention, and egress analyses;",
            "- paired comparison with the otherwise identical CAR-T control;",
            "- vector and raster figure exports with checksums.",
            "",
            "## Statistical unit",
            "",
            "Patients or independent donors are the inferential units. Cells, spots, fields, "
            "wells, and technical replicates are nested measurements and are aggregated or "
            "handled within specimen-preserving null models.",
            "",
            "## Provenance",
            "",
            "See `provenance/reproducibility_record.json` and `provenance/output_checksums.tsv`.",
            "",
        ]
    )
    report_path = output_dir / "analysis_report.md"
    report_path.write_text(report_text, encoding="utf-8", newline="\n")
    generated = sorted(
        path
        for path in output_dir.rglob("*")
        if path.is_file()
        and ".workflow" not in path.parts
        and "provenance" not in path.parts
        and "logs" not in path.parts
    )
    software_environment = [
        file_record(project_root / filename, root=project_root)
        for filename in ("requirements.lock.txt", "environment.yml")
    ]
    try:
        output_display = output_dir.relative_to(project_root).as_posix()
    except ValueError:
        output_display = output_dir.as_posix()
    provenance_payload = {
        "project": str(config.get("project", {}).get("name", "chemokine-cart")),
        "project_version": str(config.get("project", {}).get("version", "0.1.0.dev0")),
        "runtime": runtime_record(project_root, seed=seed),
        "synthetic_only": synthetic,
        "evidence_eligible": bool(manifest["evidence_eligible"].all()) and not synthetic,
        "configuration": file_record(config["_config_path"], root=project_root),
        "manifest": file_record(manifest_path, root=project_root),
        "invocation": {
            "entrypoint": "python -m chemokine_cart.pipeline",
            "stage": "report",
            "workflow_stages": list(STAGES),
            "config": file_record(config["_config_path"], root=project_root)["path"],
            "manifest": file_record(manifest_path, root=project_root)["path"],
            "output_dir": output_display,
            "random_seed": seed,
            "process_argv": [sys.executable, *sys.argv],
            "reproduction_command": [
                "python",
                "-m",
                "chemokine_cart.pipeline",
                "--stage",
                "all",
                "--config",
                file_record(config["_config_path"], root=project_root)["path"],
                "--manifest",
                file_record(manifest_path, root=project_root)["path"],
                "--output-dir",
                output_display,
                "--random-seed",
                str(seed),
            ],
        },
        "software_environment": software_environment,
        "inputs": [file_record(path, root=project_root) for path in input_paths],
        "outputs": [file_record(path, root=output_dir) for path in generated],
    }
    provenance_path = write_json(
        output_dir / "provenance" / "reproducibility_record.json",
        provenance_payload,
    )
    checksum_table = _write_table(
        pd.DataFrame(provenance_payload["outputs"]),
        output_dir / "provenance" / "output_checksums.tsv",
        sort_by=("path",),
    )
    outputs = [report_path, provenance_path, checksum_table]
    write_stage_marker(
        output_dir,
        stage="report",
        outputs=outputs,
        inputs=[
            config["_config_path"],
            manifest_path,
            project_root / "requirements.lock.txt",
            project_root / "environment.yml",
            *input_paths,
        ],
        metadata={
            "synthetic_only": synthetic,
            "evidence_eligible": provenance_payload["evidence_eligible"],
            "random_seed": seed,
        },
    )
    return outputs


def run_stage(
    stage: str,
    *,
    config_path: str | Path,
    manifest_path: str | Path,
    output_dir: str | Path,
    random_seed: int | None = None,
) -> list[Path]:
    """Execute one validated stage and return its declared outputs."""

    if stage not in STAGES:
        raise ValueError(f"stage must be one of: {', '.join(STAGES)}")
    config, project_root, manifest_file, output = _project_context(
        config_path, manifest_path, output_dir, random_seed
    )
    output.mkdir(parents=True, exist_ok=True)
    manifest = _manifest_frame(manifest_file)
    if stage == "validate":
        return _stage_validate(config, project_root, manifest_file, output)
    _enforce_evidence_policy(config, manifest)
    if stage == "single-cell":
        return _stage_single_cell(config, project_root, manifest, manifest_file, output)
    if stage == "spatial":
        return _stage_spatial(config, project_root, manifest, manifest_file, output)
    if stage == "functional":
        return _stage_functional(config, project_root, manifest, manifest_file, output)
    if stage == "figures":
        return _stage_figures(config, project_root, manifest, manifest_file, output)
    return _stage_report(config, project_root, manifest, manifest_file, output)


def run_pipeline(
    *,
    config_path: str | Path,
    manifest_path: str | Path,
    output_dir: str | Path,
    dry_run: bool = False,
    random_seed: int | None = None,
) -> list[Path]:
    """Run the complete deterministic workflow or print its execution plan."""

    if dry_run:
        config, _, manifest_file, output = _project_context(
            config_path,
            manifest_path,
            output_dir,
            random_seed,
        )
        plan = {
            "config": str(Path(config["_config_path"])),
            "manifest": str(manifest_file),
            "output_dir": str(output),
            "random_seed": _effective_seed(config),
            "stages": list(STAGES),
        }
        print(json.dumps(plan, indent=2, sort_keys=True))
        return []
    outputs: list[Path] = []
    for stage in STAGES:
        outputs.extend(
            run_stage(
                stage,
                config_path=config_path,
                manifest_path=manifest_path,
                output_dir=output_dir,
                random_seed=random_seed,
            )
        )
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run one reproducible analysis stage.")
    parser.add_argument("--stage", choices=(*STAGES, "all"), required=True)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--random-seed", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.stage == "all":
            run_pipeline(
                config_path=args.config,
                manifest_path=args.manifest,
                output_dir=args.output_dir,
                dry_run=args.dry_run,
                random_seed=args.random_seed,
            )
        elif args.dry_run:
            print(json.dumps({"stage": args.stage, "status": "planned"}, indent=2))
        else:
            run_stage(
                args.stage,
                config_path=args.config,
                manifest_path=args.manifest,
                output_dir=args.output_dir,
                random_seed=args.random_seed,
            )
    except (DataValidationError, FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
