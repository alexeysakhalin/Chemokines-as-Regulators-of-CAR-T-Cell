"""Entry point for the complete reproducible analysis workflow."""

import csv
import sys
from pathlib import Path


PROJECT_ROOT = Path(workflow.basedir).resolve()
DEFAULT_CONFIG_PATH = (PROJECT_ROOT / "config" / "config.yaml").resolve()
_override_configfiles = [Path(path).resolve() for path in workflow.overwrite_configfiles]
if len(_override_configfiles) > 1:
    raise ValueError("Use exactly one complete analysis config file.")
CONFIG_PATH = str(_override_configfiles[0] if _override_configfiles else DEFAULT_CONFIG_PATH)

configfile: CONFIG_PATH


_output_config = config.get("output", {})
_default_output = _output_config.get("root", "results") if isinstance(_output_config, dict) else "results"
_output_value = config.get("output_dir", _default_output)

_manifest_config = config.get("manifest", {})
_default_manifest = (
    _manifest_config.get("path", "config/manifest.tsv")
    if isinstance(_manifest_config, dict)
    else _manifest_config
)
_manifest_value = config.get("manifest_path", _default_manifest)

_project_config = config.get("project", {})
_default_seed = (
    _project_config.get("random_seed", 20260901)
    if isinstance(_project_config, dict)
    else 20260901
)

OUTPUT_DIR = str(PROJECT_ROOT / _output_value)
MANIFEST_PATH = str(PROJECT_ROOT / _manifest_value)
RANDOM_SEED = int(config.get("random_seed", _default_seed))
STATUS_DIR = str(Path(OUTPUT_DIR) / ".workflow")
PYTHON_EXECUTABLE = sys.executable


def _project_file(value):
    path = (PROJECT_ROOT / str(value)).resolve()
    try:
        path.relative_to(PROJECT_ROOT)
    except ValueError as exc:
        raise ValueError(f"Workflow input escapes the project root: {value}") from exc
    return str(path)


def _local_manifest_inputs(manifest_path):
    with Path(manifest_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "source_path" not in reader.fieldnames:
            raise ValueError(f"Manifest has no source_path column: {manifest_path}")
        values = {
            _project_file(row["source_path"].strip())
            for row in reader
            if row.get("source_path", "").strip()
        }
    return sorted(values)


LOCAL_INPUTS = _local_manifest_inputs(MANIFEST_PATH)
_schema_config = config.get("schemas", {})
_manifest_schema = config.get("manifest", {}).get(
    "schema", "resources/manifest.schema.tsv"
)
SCHEMA_INPUTS = sorted(
    {
        _project_file(_manifest_schema),
        *(_project_file(value) for value in _schema_config.values()),
    }
)
ANALYSIS_SOURCES = sorted(
    str(path.resolve()) for path in (PROJECT_ROOT / "src" / "chemokine_cart").glob("*.py")
)
WORKFLOW_SOURCES = [
    str((PROJECT_ROOT / "Snakefile").resolve()),
    str((PROJECT_ROOT / "workflow" / "rules" / "analysis.smk").resolve()),
]
ENVIRONMENT_INPUTS = [
    str((PROJECT_ROOT / "pyproject.toml").resolve()),
    str((PROJECT_ROOT / "requirements.lock.txt").resolve()),
    str((PROJECT_ROOT / "environment.yml").resolve()),
]

TABLE_DIR = Path(OUTPUT_DIR) / "tables"
FIGURE_DIR = Path(OUTPUT_DIR) / "figures"
PROVENANCE_DIR = Path(OUTPUT_DIR) / "provenance"

VALIDATION_OUTPUT = str(TABLE_DIR / "validation_report.json")
SINGLE_CELL_OUTPUTS = [
    str(TABLE_DIR / "single_cell_state_fractions.tsv"),
    str(TABLE_DIR / "single_cell_state_bootstrap.tsv"),
    str(TABLE_DIR / "single_cell_state_paired_changes.tsv"),
    str(TABLE_DIR / "single_cell_rna_pseudobulk.tsv"),
    str(TABLE_DIR / "single_cell_adt_pseudobulk.tsv"),
    str(TABLE_DIR / "single_cell_sample_qc.tsv"),
    str(TABLE_DIR / "single_cell_missingness.tsv"),
]
SPATIAL_OUTPUTS = [
    str(TABLE_DIR / "spatial_distances.tsv"),
    str(TABLE_DIR / "spatial_patient_summary.tsv"),
    str(TABLE_DIR / "spatial_neighborhood_enrichment.tsv"),
]
FUNCTIONAL_OUTPUTS = [
    str(TABLE_DIR / "functional_biological_summary.tsv"),
    str(TABLE_DIR / "functional_product_comparisons.tsv"),
    str(TABLE_DIR / "cxcl16_patient_measurements.tsv"),
    str(TABLE_DIR / "cxcl16_cohort_summary.tsv"),
    str(TABLE_DIR / "functional_qc.json"),
    str(TABLE_DIR / "protein_qc.json"),
]
_figure_formats = tuple(config.get("output", {}).get("figure_formats", ("pdf", "svg", "png")))
_figure_bases = (
    "figure6a_single_cell_states",
    "figure6b_spatial_map",
    "figure6c_cxcl16_protein_validation",
    "figure6d_functional_assays",
    "figure6e_matched_product_effects",
)
FIGURE_OUTPUTS = [
    str(FIGURE_DIR / f"{base}.{extension}")
    for base in _figure_bases
    for extension in _figure_formats
]
FIGURE_OUTPUTS.append(str(FIGURE_DIR / "figure_manifest.tsv"))
REPORT_OUTPUTS = [
    str(Path(OUTPUT_DIR) / "analysis_report.md"),
    str(PROVENANCE_DIR / "reproducibility_record.json"),
    str(PROVENANCE_DIR / "output_checksums.tsv"),
]
STATUS_OUTPUTS = [
    str(Path(STATUS_DIR) / f"{stage}.done")
    for stage in ("validate", "single_cell", "spatial", "functional", "figures", "report")
]
ALL_OUTPUTS = [
    *STATUS_OUTPUTS,
    VALIDATION_OUTPUT,
    *SINGLE_CELL_OUTPUTS,
    *SPATIAL_OUTPUTS,
    *FUNCTIONAL_OUTPUTS,
    *FIGURE_OUTPUTS,
    *REPORT_OUTPUTS,
]


rule all:
    input:
        ALL_OUTPUTS


include: "workflow/rules/analysis.smk"
