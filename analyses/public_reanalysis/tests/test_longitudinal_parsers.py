from __future__ import annotations

import gzip
import hashlib
import io
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from longitudinal_extension.parsers import (
    _csr_target_counts,
    parse_dense_gene_by_cell,
    parse_gse162975_matrix,
    parse_gse197268_archives,
    read_h5ad_column,
)


def test_h5ad_categorical_reader(tmp_path) -> None:
    path = tmp_path / "categorical.h5"
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        field = obs.create_group("field")
        field.create_dataset("categories", data=np.asarray([b"A", b"B"]))
        field.create_dataset("codes", data=np.asarray([1, 0, 1]))
    with h5py.File(path, "r") as handle:
        assert read_h5ad_column(handle["obs"], "field").tolist() == ["B", "A", "B"]


def test_csr_target_count_reader_handles_empty_rows(tmp_path) -> None:
    path = tmp_path / "matrix.h5"
    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.attrs["shape"] = [3, 4]
        matrix.create_dataset("indptr", data=[0, 0, 2, 3])
        matrix.create_dataset("indices", data=[1, 3, 1])
        matrix.create_dataset("data", data=[2, 4, 5])
    with h5py.File(path, "r") as handle:
        totals, targets = _csr_target_counts(handle["matrix"], [1, 3])
    assert totals.tolist() == [0, 6, 5]
    assert targets[1].tolist() == [0, 2, 5]
    assert targets[3].tolist() == [0, 4, 0]


def test_dense_parser_uses_car_positive_cells_only(monkeypatch) -> None:
    import longitudinal_extension.parsers as parsers

    monkeypatch.setattr(parsers, "MARKER_GENES", ())
    content = "gene\tc1\tc2\tc3\nCAR\t1\t0\t2\nCXCR6\t0\t5\t1\nOTHER\t4\t4\t4\n"
    cells = parse_dense_gene_by_cell(
        io.StringIO(content),
        dataset="TEST",
        patient_id="P1",
        product="axi-cel",
        stage="D7",
        biological_sample_id="sample",
        eligibility_gene="CAR",
    )
    assert cells["cell_id"].tolist() == ["c1", "c3"]
    assert cells["cxcr6_umi"].tolist() == [0, 1]
    assert cells["total_umi"].tolist() == [5, 7]


def test_dense_parser_reads_r_write_table_whitespace_layout(monkeypatch) -> None:
    import longitudinal_extension.parsers as parsers

    monkeypatch.setattr(parsers, "MARKER_GENES", ())
    content = '"c1" "c2" "c3"\n"CAR" 1 0 2\n"CXCR6" 0 5 1\n"OTHER" 4 4 4\n'
    cells = parse_dense_gene_by_cell(
        io.StringIO(content),
        dataset="TEST",
        patient_id="P1",
        product="axi-cel",
        stage="D7",
        biological_sample_id="sample",
        eligibility_gene="CAR",
    )
    assert cells["cell_id"].tolist() == ["c1", "c3"]
    assert cells["cxcr6_umi"].tolist() == [0, 1]
    assert cells["total_umi"].tolist() == [5, 7]


def test_dense_parser_rejects_noninteger_r_table_counts(monkeypatch) -> None:
    import longitudinal_extension.parsers as parsers

    monkeypatch.setattr(parsers, "MARKER_GENES", ())
    content = '"c1"\n"CAR" 1\n"CXCR6" 0.5\n'
    with pytest.raises(ValueError, match="Non-integer counts"):
        parse_dense_gene_by_cell(
            io.StringIO(content),
            dataset="TEST",
            patient_id="P1",
            product="axi-cel",
            stage="D7",
            biological_sample_id="sample",
            eligibility_gene="CAR",
        )


def test_gse162975_parser_rejects_missing_target_gene(tmp_path: Path, monkeypatch) -> None:
    import longitudinal_extension.parsers as parsers

    monkeypatch.setattr(parsers, "MARKER_GENES", ())
    header = '"gene",' + ",".join(f'"sample_sc{i}"' for i in range(7_578))
    content = f'{header}\n"OTHER",0\n'
    monkeypatch.setattr(parsers.gzip, "open", lambda *args, **kwargs: io.StringIO(content))
    monkeypatch.setattr(
        parsers.np, "fromstring", lambda *args, **kwargs: np.zeros(7_578, dtype=np.int64)
    )
    crosswalk = tmp_path / "crosswalk.tsv"
    crosswalk.write_text("technical_record\tpatient_id\tdeposited_stage\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing target genes.*CXCR6"):
        parse_gse162975_matrix(tmp_path / "matrix.csv.gz", crosswalk)


def test_standardized_cell_schema_is_ordered(synthetic_cells: pd.DataFrame) -> None:
    from longitudinal_extension.parsers import write_standard_cells

    target = io.BytesIO()
    assert list(synthetic_cells.columns)[:8] == [
        "dataset",
        "patient_id",
        "product",
        "stage",
        "biological_sample_id",
        "cell_id",
        "total_umi",
        "cxcr6_umi",
    ]
    assert callable(write_standard_cells)
    assert target.getvalue() == b""


def test_gse197268_adapter_joins_author_qc_car_and_cd8_annotations(
    tmp_path: Path, monkeypatch
) -> None:
    import longitudinal_extension.parsers as parsers

    monkeypatch.setattr(parsers, "MARKER_GENES", ())
    metadata_rows = []
    subtype_rows = []
    crosswalk_rows = []
    for patient_number in range(1, 10):
        patient_id = f"Axi-R-{patient_number:02d}"
        for canonical_stage, author_timepoint in (("IP", "Infusion"), ("D7-CART", "D7-CART")):
            eligible_id = f"eligible-{patient_number}-{canonical_stage}"
            excluded_id = f"excluded-{patient_number}-{canonical_stage}"
            metadata_rows.extend(
                [
                    {
                        "cell_id": eligible_id,
                        "barcode": patient_id,
                        "timepoint": author_timepoint,
                        "generic": "Axi-cel",
                        "CAR": True,
                        "cell_type": "T",
                    },
                    {
                        "cell_id": excluded_id,
                        "barcode": patient_id,
                        "timepoint": author_timepoint,
                        "generic": "Axi-cel",
                        "CAR": False,
                        "cell_type": "T",
                    },
                ]
            )
            subtype_rows.extend(
                [
                    {"cell_id": eligible_id, "subtype": "CD8 T"},
                    {"cell_id": excluded_id, "subtype": "CD8 T"},
                ]
            )
            source_file = f"GSM{patient_number:07d}_{canonical_stage}.tar.gz"
            archive_path = tmp_path / source_file
            members = {
                "sample/features.tsv.gz": gzip.compress(b"id\tCXCR6\tGene Expression\n"),
                "sample/barcodes.tsv.gz": gzip.compress(f"{eligible_id}\n{excluded_id}\n".encode()),
                "sample/matrix.mtx.gz": gzip.compress(
                    b"%%MatrixMarket matrix coordinate integer general\n%\n1 2 1\n1 1 1\n"
                ),
            }
            with tarfile.open(archive_path, "w:gz") as archive:
                for name, content in members.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(content)
                    archive.addfile(info, io.BytesIO(content))
            content = archive_path.read_bytes()
            crosswalk_rows.append(
                {
                    "source_file": source_file,
                    "geo_sample": f"GSM{patient_number:07d}",
                    "patient_number": patient_number,
                    "patient_id": patient_id,
                    "author_patient_label": patient_id,
                    "product": "Axi-cel",
                    "canonical_stage": canonical_stage,
                    "author_timepoint": author_timepoint,
                    "raw_matrix_cells": 2,
                    "author_qc_matched_cells": 2,
                    "expected_eligible_cells": 1,
                    "expected_cxcr6_positive_cells": 1,
                    "bytes": len(content),
                    "sha256": hashlib.sha256(content).hexdigest(),
                }
            )
    obs_path = tmp_path / "obs.csv"
    subtype_path = tmp_path / "subtype.tsv"
    crosswalk_path = tmp_path / "crosswalk.tsv"
    pd.DataFrame(metadata_rows).set_index("cell_id").to_csv(obs_path)
    pd.DataFrame(subtype_rows).set_index("cell_id").to_csv(subtype_path, sep="\t")
    pd.DataFrame(crosswalk_rows).to_csv(crosswalk_path, sep="\t", index=False)

    cells = parse_gse197268_archives(tmp_path, obs_path, subtype_path, crosswalk_path)

    assert len(cells) == 18
    assert cells["cxcr6_umi"].eq(1).all()
    assert cells["cell_id"].str.startswith("eligible-").all()
