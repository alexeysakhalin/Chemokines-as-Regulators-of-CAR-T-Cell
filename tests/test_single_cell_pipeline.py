from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import chemokine_cart.pipeline as pipeline_module
from chemokine_cart.io import (
    DataValidationError,
    load_config,
    validate_configured_inputs,
    validate_table,
)
from chemokine_cart.pipeline import (
    _aggregate_count_type,
    _enforce_evidence_policy,
    _load_tables,
    _manifest_frame,
    _require_current_stage_marker,
    _stage_single_cell,
)
from chemokine_cart.provenance import write_stage_marker


def _manifest_record(
    *,
    record_id: str,
    sample_id: str,
    patient_id: str,
    source_path: str,
) -> dict[str, object]:
    return {
        "record_id": record_id,
        "sample_id": sample_id,
        "patient_id": patient_id,
        "section_id": "NA",
        "timepoint_id": "D7",
        "modality": "citeseq",
        "assay": "single_cell_long",
        "experimental_arm": "control",
        "compartment": "blood",
        "batch_id": "B01",
        "reference_version": "REF_v1",
        "protocol_version": "PROTOCOL_v1",
        "assay_panel_version": "PANEL_v1",
        "condition": "baseline",
        "comparator": "none",
        "source_path": source_path,
        "accession": "",
        "sha256": "0" * 64,
        "data_origin": "experimental",
        "evidence_eligible": True,
        "notes": "",
    }


def test_evidence_policy_blocks_ineligible_non_synthetic_records() -> None:
    manifest = pd.DataFrame(
        {
            "record_id": ["R_ELIGIBLE", "R_QC_HOLD"],
            "data_origin": ["public", "public"],
            "evidence_eligible": [True, False],
        }
    )
    config = {"evidence": {"require_evidence_eligible_for_scientific_outputs": True}}

    with pytest.raises(DataValidationError, match="R_QC_HOLD"):
        _enforce_evidence_policy(config, manifest)


def test_evidence_policy_allows_explicitly_non_evidentiary_synthetic_demo() -> None:
    manifest = pd.DataFrame(
        {
            "record_id": ["SYNTH_RECORD"],
            "data_origin": ["synthetic"],
            "evidence_eligible": [False],
        }
    )
    config = {"evidence": {"require_evidence_eligible_for_scientific_outputs": True}}

    _enforce_evidence_policy(config, manifest)


def test_evidence_policy_rejects_missing_eligibility() -> None:
    manifest = pd.DataFrame(
        {
            "record_id": ["R_MISSING"],
            "data_origin": ["public"],
            "evidence_eligible": [None],
        }
    )

    with pytest.raises(DataValidationError, match="explicit evidence_eligible.*R_MISSING"):
        _enforce_evidence_policy({}, manifest)


def test_evidence_policy_requires_synthetic_demo_to_be_explicitly_non_evidentiary() -> None:
    manifest = pd.DataFrame(
        {
            "record_id": ["SYNTH_RECORD"],
            "data_origin": ["synthetic"],
            "evidence_eligible": [True],
        }
    )

    with pytest.raises(DataValidationError, match="explicitly marked evidence_eligible=false"):
        _enforce_evidence_policy({}, manifest)


def test_evidence_policy_rejects_synthetic_records_mixed_into_scientific_run() -> None:
    manifest = pd.DataFrame(
        {
            "record_id": ["SYNTH_RECORD", "PUBLIC_RECORD"],
            "data_origin": ["synthetic", "public"],
            "evidence_eligible": [False, True],
        }
    )

    with pytest.raises(DataValidationError, match="cannot be mixed into a scientific run"):
        _enforce_evidence_policy({}, manifest)


def test_evidence_policy_rejects_blank_data_origin() -> None:
    manifest = pd.DataFrame(
        {
            "record_id": ["R_BLANK_ORIGIN"],
            "data_origin": [""],
            "evidence_eligible": [True],
        }
    )

    with pytest.raises(DataValidationError, match="explicit data_origin.*R_BLANK_ORIGIN"):
        _enforce_evidence_policy({}, manifest)


def test_scientific_stage_enforces_evidence_policy_before_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest_path = tmp_path / "manifest.tsv"
    pd.DataFrame(
        {
            "record_id": ["R_QC_HOLD"],
            "data_origin": ["public"],
            "evidence_eligible": [False],
        }
    ).to_csv(manifest_path, sep="\t", index=False)
    config = {"evidence": {"require_evidence_eligible_for_scientific_outputs": True}}
    monkeypatch.setattr(
        pipeline_module,
        "_project_context",
        lambda *_args, **_kwargs: (config, tmp_path, manifest_path, tmp_path / "results"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_stage_single_cell",
        lambda *_args, **_kwargs: pytest.fail("scientific stage ran despite evidence hold"),
    )

    with pytest.raises(DataValidationError, match="R_QC_HOLD"):
        pipeline_module.run_stage(
            "single-cell",
            config_path=tmp_path / "config.yaml",
            manifest_path=manifest_path,
            output_dir=tmp_path / "results",
        )


def test_downstream_stage_rejects_marker_from_different_manifest(tmp_path: Path) -> None:
    output_dir = tmp_path / "results"
    output = output_dir / "tables" / "upstream.tsv"
    output.parent.mkdir(parents=True)
    output.write_text("value\n1\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("project: {}\n", encoding="utf-8")
    old_manifest = tmp_path / "old_manifest.tsv"
    old_manifest.write_text("record_id\nOLD\n", encoding="utf-8")
    current_manifest = tmp_path / "current_manifest.tsv"
    current_manifest.write_text("record_id\nCURRENT\n", encoding="utf-8")
    write_stage_marker(
        output_dir,
        stage="single-cell",
        outputs=[output],
        inputs=[config, old_manifest],
        metadata={"synthetic_only": False, "random_seed": 20260901},
    )

    with pytest.raises(DataValidationError, match="current current_manifest.tsv"):
        _require_current_stage_marker(
            output_dir,
            stage="single-cell",
            config_path=config,
            manifest_path=current_manifest,
            synthetic_only=False,
            random_seed=20260901,
        )


def test_long_schema_allows_repeated_cell_across_distinct_features(tmp_path: Path) -> None:
    table = pd.DataFrame(
        {
            "cell_id": ["C1", "C1"],
            "sample_id": ["S1", "S1"],
            "patient_id": ["P1", "P1"],
            "section_id": ["NA", "NA"],
            "timepoint_id": ["D7", "D7"],
            "modality": ["citeseq", "citeseq"],
            "t_cell_origin": ["CAR_T", "CAR_T"],
            "car_positive": ["true", "true"],
            "cell_state": ["memory", "memory"],
            "annotation_version": ["STATE_v1", "STATE_v1"],
            "car_gate_version": ["GATE_v1", "GATE_v1"],
            "feature_id": ["CD3D", "CXCR3"],
            "feature_type": ["RNA", "ADT"],
            "count": [2, 5],
            "data_origin": ["experimental", "experimental"],
        }
    )
    path = tmp_path / "long.tsv"
    table.to_csv(path, sep="\t", index=False)
    schema = Path(__file__).parents[1] / "resources" / "single_cell.schema.tsv"

    report = validate_table(path, schema, allow_extra_columns=True)

    assert report.rows == 2


def test_configured_single_cell_column_names_are_validated() -> None:
    project_root = Path(__file__).parents[1]
    config = load_config(project_root / "config" / "config.yaml")
    config["single_cell"]["state_column"] = "missing_configured_state"

    with pytest.raises(DataValidationError, match="configured single-cell columns missing"):
        validate_configured_inputs(
            config,
            manifest_override="data/demo/manifest.tsv",
            check_files=True,
        )


def test_manifest_linkage_filters_shared_source_once_and_validates_identity(
    tmp_path: Path,
) -> None:
    source = tmp_path / "cells.tsv"
    pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "patient_id": ["P1", "P2"],
            "section_id": ["NA", "NA"],
            "timepoint_id": ["D7", "D7"],
            "modality": ["citeseq", "citeseq"],
            "data_origin": ["experimental", "experimental"],
        }
    ).to_csv(source, sep="\t", index=False)
    manifest = pd.DataFrame(
        [
            _manifest_record(
                record_id="R1",
                sample_id="S1",
                patient_id="P1",
                source_path="cells.tsv",
            ),
            _manifest_record(
                record_id="R2",
                sample_id="S2",
                patient_id="P2",
                source_path="cells.tsv",
            ),
        ]
    )

    loaded = _load_tables(manifest, tmp_path, modalities=("citeseq",))

    assert len(loaded) == 2
    assert set(loaded["sample_id"]) == {"S1", "S2"}
    assert set(loaded["experimental_arm"]) == {"control"}

    manifest.loc[manifest["sample_id"] == "S2", "patient_id"] = "WRONG"
    with pytest.raises(DataValidationError, match="do not match manifest"):
        _load_tables(manifest, tmp_path, modalities=("citeseq",))


def test_duplicate_cell_feature_rows_are_rejected() -> None:
    table = pd.DataFrame(
        {
            "sample_id": ["S1", "S1"],
            "cell_id": ["C1", "C1"],
            "patient_id": ["P1", "P1"],
            "timepoint_id": ["D7", "D7"],
            "t_cell_origin": ["CAR_T", "CAR_T"],
            "cell_state": ["memory", "memory"],
            "feature_id": ["CXCR3", "CXCR3"],
            "feature_type": ["RNA", "RNA"],
            "count": [2, 2],
        }
    )

    with pytest.raises(DataValidationError, match="duplicate single-cell feature"):
        _aggregate_count_type(
            table,
            feature_type="RNA",
            group_columns=(
                "sample_id",
                "patient_id",
                "timepoint_id",
                "t_cell_origin",
                "cell_state",
            ),
            metadata_columns=(
                "patient_id",
                "timepoint_id",
                "t_cell_origin",
                "cell_state",
            ),
            min_cells=1,
            long_format=True,
        )


def test_demo_stage_separates_origins_modalities_and_reports_qc(tmp_path: Path) -> None:
    project_root = Path(__file__).parents[1]
    config = load_config(project_root / "config" / "config.yaml")
    manifest_path = project_root / "data" / "demo" / "manifest.tsv"
    manifest = _manifest_frame(manifest_path)

    outputs = _stage_single_cell(
        config,
        project_root,
        manifest,
        manifest_path,
        tmp_path,
    )

    output_names = {path.name for path in outputs}
    assert "single_cell_rna_pseudobulk.tsv" in output_names
    assert "single_cell_adt_pseudobulk.tsv" in output_names
    assert "single_cell_pseudobulk.tsv" not in output_names

    qc = pd.read_csv(tmp_path / "tables" / "single_cell_sample_qc.tsv", sep="\t")
    product_origins = set(qc.loc[qc["compartment"] == "product", "t_cell_origin"])
    assert product_origins == {"CAR_T", "CAR_negative_product_T"}
    blood_origins = set(qc.loc[qc["compartment"] == "blood", "t_cell_origin"])
    assert blood_origins == {"CAR_T", "endogenous_T"}

    fractions = pd.read_csv(tmp_path / "tables" / "single_cell_state_fractions.tsv", sep="\t")
    totals = fractions.groupby(["sample_id", "t_cell_origin", "t_cell_subtype"], observed=True)[
        "fraction"
    ].sum()
    assert (totals.round(12) == 1.0).all()
    assert "other" in set(fractions["cell_state"])

    missingness = pd.read_csv(tmp_path / "tables" / "single_cell_missingness.tsv", sep="\t")
    first_patient = missingness.loc[
        missingness["patient_id"] == missingness["patient_id"].iloc[0],
        "timepoint_id",
    ]
    assert first_patient.tolist() == ["PRODUCT", "T0", "T1"]

    paired = pd.read_csv(
        tmp_path / "tables" / "single_cell_state_paired_changes.tsv",
        sep="\t",
    )
    assert paired["comparison_timepoint"].drop_duplicates().tolist() == ["T0", "T1"]


def test_stage_records_strata_when_all_samples_fail_cell_threshold(tmp_path: Path) -> None:
    project_root = Path(__file__).parents[1]
    config = load_config(project_root / "config" / "config.yaml")
    config["single_cell"]["minimum_cells_per_group"] = 1_000
    manifest_path = project_root / "data" / "demo" / "manifest.tsv"
    manifest = _manifest_frame(manifest_path)

    _stage_single_cell(
        config,
        project_root,
        manifest,
        manifest_path,
        tmp_path,
    )

    bootstrap = pd.read_csv(
        tmp_path / "tables" / "single_cell_state_bootstrap.tsv",
        sep="\t",
    )
    assert not bootstrap.empty
    assert set(bootstrap["status"]) == {"no_eligible_patients"}
    assert set(bootstrap["n_patients"]) == {0}
    assert bootstrap[["estimate", "ci_low", "ci_high"]].isna().all().all()
