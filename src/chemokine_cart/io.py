"""Input contracts, manifest validation, and deterministic demonstration data."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
FALSE_VALUES = {"0", "false", "no", "n"}
TRUE_VALUES = {"1", "true", "yes", "y"}


class DataValidationError(ValueError):
    """Raised when an input violates a declared data contract."""

    def __init__(self, errors: Sequence[str]):
        self.errors = tuple(errors)
        super().__init__("\n".join(self.errors))


@dataclass(frozen=True)
class ColumnRule:
    """One column definition read from a tab-separated schema file."""

    name: str
    data_type: str
    required: bool
    nullable: bool
    unique: bool
    pattern: str | None
    allowed: tuple[str, ...]
    description: str


@dataclass(frozen=True)
class ValidationReport:
    """Machine-readable result of validating one table."""

    path: Path
    rows: int
    columns: tuple[str, ...]
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path.as_posix(),
            "rows": self.rows,
            "columns": list(self.columns),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class ProjectValidationReport:
    """Validation result for a manifest and its canonical local tables."""

    manifest: ValidationReport
    tables: tuple[ValidationReport, ...] = ()
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        return {
            "manifest": self.manifest.as_dict(),
            "tables": [report.as_dict() for report in self.tables],
            "warnings": list(self.warnings),
        }


def _parse_bool(value: str, *, field: str) -> bool:
    normalized = value.strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ValueError(f"{field} must be a boolean value")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a UTF-8 TSV while preserving all values as strings."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise DataValidationError([f"{path}: missing header row"])
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise DataValidationError([f"{path}: duplicate column names"])
        rows = [dict(row) for row in reader]
    return list(reader.fieldnames), rows


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    """Write a deterministic UTF-8 TSV with Unix line endings."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def read_schema(path: Path) -> tuple[ColumnRule, ...]:
    """Read and validate the compact TSV schema format used by the project."""

    required_schema_columns = {
        "column",
        "type",
        "required",
        "nullable",
        "unique",
        "pattern",
        "allowed",
        "description",
    }
    header, rows = read_tsv(path)
    missing = sorted(required_schema_columns.difference(header))
    if missing:
        raise DataValidationError([f"{path}: schema columns missing: {', '.join(missing)}"])

    errors: list[str] = []
    rules: list[ColumnRule] = []
    names: set[str] = set()
    for index, row in enumerate(rows, start=2):
        name = row["column"].strip()
        data_type = row["type"].strip().lower()
        if not name:
            errors.append(f"{path}:{index}: empty schema column name")
            continue
        if name in names:
            errors.append(f"{path}:{index}: duplicate rule for {name!r}")
        names.add(name)
        if data_type not in {"string", "integer", "number", "boolean"}:
            errors.append(f"{path}:{index}: unsupported type {data_type!r}")
        try:
            required = _parse_bool(row["required"], field="required")
            nullable = _parse_bool(row["nullable"], field="nullable")
            unique = _parse_bool(row["unique"], field="unique")
        except ValueError as exc:
            errors.append(f"{path}:{index}: {exc}")
            continue
        pattern = row["pattern"].strip() or None
        if pattern:
            try:
                re.compile(pattern)
            except re.error as exc:
                errors.append(f"{path}:{index}: invalid regular expression: {exc}")
        allowed = tuple(part for part in row["allowed"].split("|") if part)
        rules.append(
            ColumnRule(
                name=name,
                data_type=data_type,
                required=required,
                nullable=nullable,
                unique=unique,
                pattern=pattern,
                allowed=allowed,
                description=row["description"].strip(),
            )
        )
    if errors:
        raise DataValidationError(errors)
    return tuple(rules)


def _validate_scalar(value: str, rule: ColumnRule, location: str) -> list[str]:
    errors: list[str] = []
    if value == "":
        if not rule.nullable:
            errors.append(f"{location}: {rule.name} may not be empty")
        return errors

    if rule.data_type == "integer":
        try:
            int(value)
        except ValueError:
            errors.append(f"{location}: {rule.name} must be an integer")
    elif rule.data_type == "number":
        try:
            number = float(value)
            if not math.isfinite(number):
                raise ValueError
        except ValueError:
            errors.append(f"{location}: {rule.name} must be a finite number")
    elif rule.data_type == "boolean":
        try:
            _parse_bool(value, field=rule.name)
        except ValueError as exc:
            errors.append(f"{location}: {exc}")

    if rule.pattern and re.fullmatch(rule.pattern, value) is None:
        errors.append(f"{location}: {rule.name} does not match {rule.pattern!r}")
    if rule.allowed and value not in rule.allowed:
        errors.append(f"{location}: {rule.name} must be one of {', '.join(rule.allowed)}")
    return errors


def validate_table(
    path: Path,
    schema_path: Path,
    *,
    allow_extra_columns: bool = True,
) -> ValidationReport:
    """Validate a TSV against a schema and report all detected errors."""

    header, rows = read_tsv(path)
    rules = read_schema(schema_path)
    errors: list[str] = []
    required = {rule.name for rule in rules if rule.required}
    missing = sorted(required.difference(header))
    if missing:
        errors.append(f"{path}: required columns missing: {', '.join(missing)}")
    known = {rule.name for rule in rules}
    extra = sorted(set(header).difference(known))
    if extra and not allow_extra_columns:
        errors.append(f"{path}: undeclared columns present: {', '.join(extra)}")

    unique_values: dict[str, dict[str, int]] = {
        rule.name: {} for rule in rules if rule.unique and rule.name in header
    }
    for row_number, row in enumerate(rows, start=2):
        location = f"{path}:{row_number}"
        for rule in rules:
            if rule.name not in header:
                continue
            value = (row.get(rule.name) or "").strip()
            errors.extend(_validate_scalar(value, rule, location))
            if rule.unique and value:
                previous = unique_values[rule.name].get(value)
                if previous is not None:
                    errors.append(
                        f"{location}: duplicate {rule.name} {value!r}; first seen on row {previous}"
                    )
                else:
                    unique_values[rule.name][value] = row_number

    if errors:
        raise DataValidationError(errors)
    warnings = tuple(f"undeclared column retained: {name}" for name in extra)
    return ValidationReport(path=path, rows=len(rows), columns=tuple(header), warnings=warnings)


def resolve_project_path(project_root: Path, value: str, *, must_exist: bool = False) -> Path:
    """Resolve a project-relative path and reject absolute or escaping paths."""

    candidate = Path(value)
    if candidate.is_absolute():
        raise DataValidationError([f"path must be relative to the project root: {value}"])
    root = project_root.resolve()
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise DataValidationError([f"path escapes the project root: {value}"]) from exc
    if must_exist and not resolved.exists():
        raise DataValidationError([f"path does not exist: {value}"])
    return resolved


def load_config(config_path: Path) -> dict[str, Any]:
    """Load YAML configuration and attach its resolved project root."""

    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - dependency error is environment-specific
        raise RuntimeError("PyYAML is required to read config files") from exc

    with config_path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise DataValidationError([f"{config_path}: top-level YAML value must be a mapping"])
    project = loaded.get("project")
    if not isinstance(project, dict):
        raise DataValidationError([f"{config_path}: missing project mapping"])
    root_value = str(project.get("root", ".."))
    root_candidate = Path(root_value)
    if root_candidate.is_absolute():
        raise DataValidationError([f"{config_path}: project.root must be relative"])
    project_root = (config_path.parent / root_candidate).resolve()
    loaded["_config_path"] = config_path.resolve()
    loaded["_project_root"] = project_root
    return loaded


def validate_manifest(
    manifest_path: Path,
    schema_path: Path,
    project_root: Path,
    *,
    check_files: bool = False,
    verify_checksums: bool = False,
) -> ValidationReport:
    """Validate identifiers, provenance, relative paths, and file checksums."""

    base_report = validate_table(manifest_path, schema_path, allow_extra_columns=False)
    _, rows = read_tsv(manifest_path)
    errors: list[str] = []
    warnings = list(base_report.warnings)

    identity_columns = ("record_id", "sample_id", "patient_id", "section_id", "timepoint_id")
    for row_number, row in enumerate(rows, start=2):
        location = f"{manifest_path}:{row_number}"
        for column in identity_columns:
            value = row[column].strip()
            if ID_PATTERN.fullmatch(value) is None:
                errors.append(f"{location}: invalid {column} {value!r}")

        source_path = row["source_path"].strip()
        accession = row["accession"].strip()
        checksum = row["sha256"].strip().lower()
        origin = row["data_origin"].strip()
        evidence_eligible = row["evidence_eligible"].strip().lower() in TRUE_VALUES
        if not source_path and not accession:
            errors.append(f"{location}: source_path or accession is required")
        if origin == "synthetic" and evidence_eligible:
            errors.append(f"{location}: synthetic data must have evidence_eligible=false")
        if source_path:
            try:
                resolved = resolve_project_path(project_root, source_path)
            except DataValidationError as exc:
                errors.extend(f"{location}: {message}" for message in exc.errors)
                continue
            if not checksum:
                errors.append(f"{location}: sha256 is required for a local source_path")
            elif SHA256_PATTERN.fullmatch(checksum) is None:
                errors.append(
                    f"{location}: sha256 must contain 64 lowercase hexadecimal characters"
                )
            if check_files and not resolved.is_file():
                errors.append(f"{location}: source file not found: {source_path}")
            elif verify_checksums and resolved.is_file() and checksum:
                observed = sha256_file(resolved)
                if observed != checksum:
                    errors.append(
                        f"{location}: checksum mismatch for {source_path}; "
                        f"manifest={checksum}, observed={observed}"
                    )
        elif checksum:
            warnings.append(f"{location}: checksum supplied without a local source_path")

    if errors:
        raise DataValidationError(errors)
    return ValidationReport(
        path=base_report.path,
        rows=base_report.rows,
        columns=base_report.columns,
        warnings=tuple(warnings),
    )


def validate_configured_manifest(
    config: Mapping[str, Any],
    *,
    manifest_override: str | None = None,
    check_files: bool | None = None,
) -> ValidationReport:
    """Resolve and validate the manifest declared by a loaded config."""

    root = Path(config["_project_root"])
    manifest_config = config.get("manifest")
    if not isinstance(manifest_config, Mapping):
        raise DataValidationError(["config: missing manifest mapping"])
    manifest_value = manifest_override or str(manifest_config.get("path", ""))
    schema_value = str(manifest_config.get("schema", ""))
    if not manifest_value or not schema_value:
        raise DataValidationError(["config: manifest.path and manifest.schema are required"])
    manifest_path = resolve_project_path(root, manifest_value, must_exist=True)
    schema_path = resolve_project_path(root, schema_value, must_exist=True)
    should_check = (
        bool(manifest_config.get("check_files", False)) if check_files is None else check_files
    )
    verify = bool(manifest_config.get("verify_checksums", True)) and should_check
    return validate_manifest(
        manifest_path,
        schema_path,
        root,
        check_files=should_check,
        verify_checksums=verify,
    )


def validate_configured_inputs(
    config: Mapping[str, Any],
    *,
    manifest_override: str | None = None,
    check_files: bool | None = None,
) -> ProjectValidationReport:
    """Validate the manifest and any local canonical TSV inputs it declares."""

    manifest_report = validate_configured_manifest(
        config,
        manifest_override=manifest_override,
        check_files=check_files,
    )
    manifest_config = config["manifest"]
    should_check = (
        bool(manifest_config.get("check_files", False)) if check_files is None else check_files
    )
    if not should_check:
        return ProjectValidationReport(manifest=manifest_report)

    root = Path(config["_project_root"])
    schema_config = config.get("schemas")
    if not isinstance(schema_config, Mapping):
        raise DataValidationError(["config: missing schemas mapping"])
    schema_key_by_modality = {
        "scrna": "single_cell",
        "citeseq": "single_cell",
        "spatial": "spatial",
        "protein": "protein",
        "functional": "functional",
    }
    _, rows = read_tsv(manifest_report.path)
    reports: list[ValidationReport] = []
    warnings: list[str] = []
    errors: list[str] = []
    records_by_source: dict[Path, list[tuple[int, dict[str, str]]]] = {}
    for row_number, row in enumerate(rows, start=2):
        source_value = row["source_path"].strip()
        if not source_value:
            continue
        source_path = resolve_project_path(root, source_value, must_exist=True)
        records_by_source.setdefault(source_path, []).append((row_number, row))

    for source_path, source_records in records_by_source.items():
        first_row_number, first_record = source_records[0]
        source_value = first_record["source_path"].strip()
        if source_path.suffix.lower() != ".tsv":
            warnings.append(
                f"{manifest_report.path}:{first_row_number}: no tabular schema applied to "
                f"{source_value}; checksum validation completed"
            )
            continue

        modalities = {record["modality"].strip() for _, record in source_records}
        schema_keys = {schema_key_by_modality.get(modality) for modality in modalities}
        if None in schema_keys or any(key not in schema_config for key in schema_keys):
            errors.append(
                f"{manifest_report.path}:{first_row_number}: no configured schema for "
                f"modalities {sorted(modalities)!r}"
            )
            continue
        if len(schema_keys) != 1:
            errors.append(
                f"{source_path}: one source file may not mix modalities governed by "
                "different table schemas"
            )
            continue
        schema_key = next(iter(schema_keys))
        assert schema_key is not None
        schema_path = resolve_project_path(root, str(schema_config[schema_key]), must_exist=True)
        try:
            table_report = validate_table(source_path, schema_path, allow_extra_columns=True)
        except DataValidationError as exc:
            errors.extend(exc.errors)
            continue
        reports.append(table_report)
        warnings.extend(table_report.warnings)

        header, table_rows = read_tsv(source_path)
        synthetic_multi = (
            len(source_records) == 1
            and first_record["data_origin"].strip().lower() == "synthetic"
            and first_record["sample_id"].strip().upper() == "SYNTH_COHORT"
        )
        if synthetic_multi:
            for column in ("data_origin", "modality"):
                if column not in header:
                    continue
                observed = {item.get(column, "").strip() for item in table_rows}
                expected = {first_record[column].strip()}
                if observed != expected:
                    errors.append(
                        f"{source_path}: {column} values {sorted(observed)!r} do not match "
                        f"manifest value {first_record[column]!r}"
                    )
        else:
            manifest_samples = [record["sample_id"].strip() for _, record in source_records]
            if len(manifest_samples) != len(set(manifest_samples)):
                errors.append(f"{source_path}: duplicate manifest records for the same sample_id")
            if "sample_id" not in header:
                errors.append(f"{source_path}: sample_id is required for manifest linkage")
            else:
                observed_samples = {item["sample_id"].strip() for item in table_rows}
                expected_samples = set(manifest_samples)
                missing = sorted(expected_samples.difference(observed_samples))
                unregistered = sorted(observed_samples.difference(expected_samples))
                if missing:
                    errors.append(
                        f"{source_path}: manifest samples absent from table: {', '.join(missing)}"
                    )
                if unregistered:
                    errors.append(
                        f"{source_path}: table contains unregistered samples: "
                        f"{', '.join(unregistered)}"
                    )
                for _, record in source_records:
                    sample_rows = [
                        item
                        for item in table_rows
                        if item["sample_id"].strip() == record["sample_id"].strip()
                    ]
                    for column in (
                        "patient_id",
                        "section_id",
                        "timepoint_id",
                        "data_origin",
                        "modality",
                    ):
                        if column not in header or not sample_rows:
                            continue
                        observed = {item[column].strip() for item in sample_rows}
                        expected = {record[column].strip()}
                        if observed != expected:
                            errors.append(
                                f"{source_path}: sample {record['sample_id']!r} has {column} "
                                f"values {sorted(observed)!r}, expected {record[column]!r}"
                            )

        if modalities.intersection({"scrna", "citeseq"}):
            single_cell_config = config.get("single_cell", {})
            dynamic_columns = (
                str(single_cell_config.get("car_column", "car_positive")),
                str(single_cell_config.get("origin_column", "t_cell_origin")),
                str(single_cell_config.get("state_column", "cell_state")),
            )
            missing_dynamic = sorted(set(dynamic_columns).difference(header))
            if missing_dynamic:
                errors.append(
                    f"{source_path}: configured single-cell columns missing: "
                    f"{', '.join(missing_dynamic)}"
                )

            long_columns = {"feature_id", "feature_type", "count"}
            present_long = long_columns.intersection(header)
            if present_long and present_long != long_columns:
                missing_long = sorted(long_columns.difference(header))
                errors.append(
                    f"{source_path}: incomplete long-format count contract; missing "
                    f"{', '.join(missing_long)}"
                )
            elif present_long == long_columns:
                first_seen: dict[tuple[str, str, str, str], int] = {}
                for table_row_number, item in enumerate(table_rows, start=2):
                    key = (
                        item["sample_id"].strip(),
                        item["cell_id"].strip(),
                        item["feature_type"].strip(),
                        item["feature_id"].strip(),
                    )
                    previous = first_seen.get(key)
                    if previous is not None:
                        errors.append(
                            f"{source_path}:{table_row_number}: duplicate single-cell "
                            "sample_id/cell_id/feature_type/feature_id row; first seen on "
                            f"row {previous}"
                        )
                    else:
                        first_seen[key] = table_row_number

    if errors:
        raise DataValidationError(errors)
    return ProjectValidationReport(
        manifest=manifest_report,
        tables=tuple(reports),
        warnings=tuple(warnings),
    )


def _demo_single_cell_rows(rng: random.Random) -> list[dict[str, Any]]:
    states = ("memory", "effector", "dysfunctional")
    rows: list[dict[str, Any]] = []
    for patient_index in range(1, 4):
        patient_id = f"SYNTH_P{patient_index:02d}"
        for time_index, timepoint in enumerate(("PRODUCT", "T0", "T1")):
            compartment = "PRODUCT" if timepoint == "PRODUCT" else "PB"
            sample_id = f"{patient_id}_{timepoint}_{compartment}"
            for cell_index in range(1, 13):
                state = states[(cell_index + patient_index + time_index) % len(states)]
                rows.append(
                    {
                        "cell_id": f"{sample_id}_C{cell_index:03d}",
                        "sample_id": sample_id,
                        "patient_id": patient_id,
                        "section_id": "NA",
                        "timepoint_id": timepoint,
                        "compartment": "product" if timepoint == "PRODUCT" else "blood",
                        "modality": "citeseq",
                        "t_cell_origin": (
                            "CAR_T"
                            if cell_index <= 8
                            else (
                                "CAR_negative_product_T"
                                if timepoint == "PRODUCT"
                                else "endogenous_T"
                            )
                        ),
                        "car_positive": "true" if cell_index <= 8 else "false",
                        "cell_state": state,
                        "t_cell_subtype": "CD4" if cell_index % 2 else "CD8",
                        "annotation_version": "SYNTH_STATE_v1",
                        "car_gate_version": "SYNTH_CAR_GATE_v1",
                        "CCL5_RNA_count": int(rng.random() * 25),
                        "CXCR3_RNA_count": int(rng.random() * 18),
                        "CXCR3_ADT_count": int(rng.random() * 40),
                        "data_origin": "synthetic",
                    }
                )
    return rows


def _demo_spatial_rows(rng: random.Random) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    compartments = ("vessel", "stroma", "tumor_nest")
    cell_types = ("endothelial", "fibroblast", "tumor", "CAR_T")
    chemokines = ("CXCL9", "CXCL10", "CXCL12", "CXCL16")
    for index in range(1, 49):
        cell_type = cell_types[(index - 1) % len(cell_types)]
        compartment = compartments[(index - 1) % len(compartments)]
        rows.append(
            {
                "cell_id": f"SYNTH_SP_CELL{index:03d}",
                "sample_id": "SYNTH_P01_T1_TUMOR",
                "patient_id": "SYNTH_P01",
                "section_id": "SYNTH_SEC01",
                "timepoint_id": "T1",
                "x": f"{rng.uniform(0, 800):.6f}",
                "y": f"{rng.uniform(0, 800):.6f}",
                "cell_type": cell_type,
                "is_ligand_source": "true"
                if cell_type in {"endothelial", "fibroblast"}
                else "false",
                "is_car_t": "true" if cell_type == "CAR_T" else "false",
                "structure": compartment if cell_type != "CAR_T" else "",
                "chemokine": chemokines[(index - 1) % len(chemokines)],
                "chemokine_signal": f"{rng.uniform(0, 10):.6f}",
                "data_origin": "synthetic",
            }
        )
    return rows


def _demo_protein_rows(rng: random.Random) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    analytes = (
        ("CXCL9", "soluble"),
        ("CXCL10", "soluble"),
        ("CXCL16", "soluble"),
        ("CXCL16", "membrane"),
    )
    for patient_index in range(1, 4):
        for timepoint in ("T0", "T1"):
            for analyte, form in analytes:
                rows.append(
                    {
                        "measurement_id": f"SYNTH_PROT_{patient_index}_{timepoint}_{analyte}_{form}",
                        "sample_id": f"SYNTH_P{patient_index:02d}_{timepoint}_PB",
                        "patient_id": f"SYNTH_P{patient_index:02d}",
                        "section_id": "NA",
                        "timepoint_id": timepoint,
                        "analyte": analyte,
                        "molecular_form": form,
                        "response": f"{rng.uniform(5, 200):.6f}",
                        "unit": "pg_per_mL" if form == "soluble" else "MFI",
                        "assay": "immunoassay" if form == "soluble" else "flow_cytometry",
                        "biological_replicate": patient_index,
                        "technical_replicate": 1,
                        "data_origin": "synthetic",
                    }
                )
    return rows


def _demo_functional_rows(rng: random.Random) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    assays = ("chemotaxis", "retention", "egress", "cxcl16")
    products = ("identical_control", "receptor_matched")
    for patient_index in range(1, 4):
        for assay in assays:
            forms = ("soluble", "membrane") if assay == "cxcl16" else ("",)
            for product in products:
                for form in forms:
                    form_id = form or "NA"
                    condition = f"CXCL16_{form}" if assay == "cxcl16" else "CXCL10_gradient"
                    rows.append(
                        {
                            "measurement_id": (
                                f"SYNTH_FUNC_{patient_index}_{assay}_{product}_{form_id}"
                            ),
                            "sample_id": f"SYNTH_P{patient_index:02d}_T1_PRODUCT",
                            "patient_id": f"SYNTH_P{patient_index:02d}",
                            "section_id": "NA",
                            "timepoint_id": "T1",
                            "assay": assay,
                            "product": product,
                            "condition": condition,
                            "pair_id": f"SYNTH_PAIR_{patient_index}_{assay}_{form_id}",
                            "chemokine": "CXCL10" if assay != "cxcl16" else "CXCL16",
                            "dose": "100",
                            "dose_unit": "ng_per_mL",
                            "time": "4",
                            "time_unit": "hour",
                            "cxcl16_form": form,
                            "response": f"{rng.uniform(0.1, 1.0):.6f}",
                            "unit": "fraction_of_input",
                            "biological_replicate": patient_index,
                            "technical_replicate": 1,
                            "data_origin": "synthetic",
                        }
                    )
    return rows


def make_demo_dataset(
    project_root: Path,
    output_dir: str = "data/demo",
    *,
    seed: int = 20260901,
    force: bool = False,
) -> Path:
    """Create deterministic synthetic inputs for testing, never for evidence."""

    target = resolve_project_path(project_root, output_dir)
    if target.exists() and not target.is_dir():
        raise FileExistsError(f"{output_dir} exists and is not a directory")
    if target.exists() and any(target.iterdir()) and not force:
        raise FileExistsError(f"{output_dir} is not empty; pass --force to replace demo files")
    target.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)

    tables: tuple[tuple[str, Sequence[str], list[dict[str, Any]], str, str], ...] = (
        (
            "single_cell.tsv",
            (
                "cell_id",
                "sample_id",
                "patient_id",
                "section_id",
                "timepoint_id",
                "compartment",
                "modality",
                "t_cell_origin",
                "car_positive",
                "cell_state",
                "t_cell_subtype",
                "annotation_version",
                "car_gate_version",
                "CCL5_RNA_count",
                "CXCR3_RNA_count",
                "CXCR3_ADT_count",
                "data_origin",
            ),
            _demo_single_cell_rows(rng),
            "citeseq",
            "single_cell_long",
        ),
        (
            "spatial.tsv",
            (
                "cell_id",
                "sample_id",
                "patient_id",
                "section_id",
                "timepoint_id",
                "x",
                "y",
                "cell_type",
                "is_ligand_source",
                "is_car_t",
                "structure",
                "chemokine",
                "chemokine_signal",
                "data_origin",
            ),
            _demo_spatial_rows(rng),
            "spatial",
            "spatial_objects",
        ),
        (
            "protein.tsv",
            (
                "measurement_id",
                "sample_id",
                "patient_id",
                "section_id",
                "timepoint_id",
                "analyte",
                "molecular_form",
                "response",
                "unit",
                "assay",
                "biological_replicate",
                "technical_replicate",
                "data_origin",
            ),
            _demo_protein_rows(rng),
            "protein",
            "protein_gradient",
        ),
        (
            "functional.tsv",
            (
                "measurement_id",
                "sample_id",
                "patient_id",
                "section_id",
                "timepoint_id",
                "assay",
                "product",
                "condition",
                "pair_id",
                "chemokine",
                "dose",
                "dose_unit",
                "time",
                "time_unit",
                "cxcl16_form",
                "response",
                "unit",
                "biological_replicate",
                "technical_replicate",
                "data_origin",
            ),
            _demo_functional_rows(rng),
            "functional",
            "functional_assay",
        ),
    )

    manifest_rows: list[dict[str, Any]] = []
    relative_target = target.relative_to(project_root.resolve())
    for record_index, (filename, columns, rows, modality, assay) in enumerate(tables, start=1):
        table_path = target / filename
        write_tsv(table_path, columns, rows)
        manifest_rows.append(
            {
                "record_id": f"SYNTH_RECORD_{record_index:02d}",
                "sample_id": "SYNTH_COHORT",
                "patient_id": "SYNTH_MULTI",
                "section_id": "SYNTH_MIXED",
                "timepoint_id": "SYNTH_MIXED",
                "modality": modality,
                "assay": assay,
                "experimental_arm": "MULTI",
                "compartment": "MULTI",
                "batch_id": "SYNTH_BATCH_01",
                "reference_version": "SYNTH_REFERENCE_v1",
                "protocol_version": "SYNTH_PROTOCOL_v1",
                "assay_panel_version": "SYNTH_PANEL_v1",
                "condition": "synthetic_software_test",
                "comparator": "synthetic_software_test",
                "source_path": (relative_target / filename).as_posix(),
                "accession": "",
                "sha256": sha256_file(table_path),
                "data_origin": "synthetic",
                "evidence_eligible": "false",
                "notes": "Synthetic software-test data; never use as scientific evidence.",
            }
        )

    manifest_path = target / "manifest.tsv"
    write_tsv(
        manifest_path,
        (
            "record_id",
            "sample_id",
            "patient_id",
            "section_id",
            "timepoint_id",
            "modality",
            "assay",
            "experimental_arm",
            "compartment",
            "batch_id",
            "reference_version",
            "protocol_version",
            "assay_panel_version",
            "condition",
            "comparator",
            "source_path",
            "accession",
            "sha256",
            "data_origin",
            "evidence_eligible",
            "notes",
        ),
        manifest_rows,
    )
    notice = {
        "data_origin": "synthetic",
        "evidence_eligible": False,
        "purpose": "software testing only",
        "seed": seed,
    }
    (target / "SYNTHETIC_DATA_NOTICE.json").write_text(
        json.dumps(notice, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest_path
