"""Auditable readers for the public processed single-cell inputs."""

from __future__ import annotations

import gzip
import io
import re
import tarfile
from collections.abc import Iterable
from pathlib import Path, PurePosixPath
from typing import BinaryIO

import h5py
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix

from .analysis import validate_cells
from .common import sha256
from .constants import CELL_COLUMNS, MARKER_GENES


def _decode(values: np.ndarray) -> np.ndarray:
    """Decode an HDF5 string array without changing non-string values."""
    array = np.asarray(values)
    if array.dtype.kind in {"O", "S"}:
        return np.asarray(
            [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in array]
        )
    return array


def read_h5ad_column(group: h5py.Group, column: str) -> np.ndarray:
    """Read an AnnData dataframe column encoded as an array or categorical."""
    node = group[column]
    if isinstance(node, h5py.Dataset):
        return _decode(node[...])
    if not isinstance(node, h5py.Group) or "codes" not in node or "categories" not in node:
        raise ValueError(f"Unsupported H5AD encoding for {group.name}/{column}")
    codes = np.asarray(node["codes"][...], dtype=np.int64)
    categories = _decode(node["categories"][...])
    if (codes < 0).any():
        raise ValueError(f"Missing categorical values in {group.name}/{column}")
    return categories[codes]


def _csr_target_counts(
    matrix: h5py.Group, target_indices: Iterable[int]
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """Read row totals and selected columns from an H5AD CSR matrix."""
    shape = tuple(int(value) for value in matrix.attrs["shape"])
    n_rows, _ = shape
    indptr = np.asarray(matrix["indptr"][...], dtype=np.int64)
    if len(indptr) != n_rows + 1 or indptr[0] != 0 or (np.diff(indptr) < 0).any():
        raise ValueError("Invalid CSR index pointer")
    totals = np.zeros(n_rows, dtype=np.int64)
    targets = {int(index): np.zeros(n_rows, dtype=np.int64) for index in target_indices}
    chunk_rows = 4096
    for first_row in range(0, n_rows, chunk_rows):
        last_row = min(first_row + chunk_rows, n_rows)
        first_value = int(indptr[first_row])
        last_value = int(indptr[last_row])
        indices = np.asarray(matrix["indices"][first_value:last_value], dtype=np.int64)
        raw_data = np.asarray(matrix["data"][first_value:last_value])
        if not np.isfinite(raw_data).all() or (raw_data < 0).any():
            raise ValueError("H5AD matrix contains invalid counts")
        rounded = np.rint(raw_data)
        if not np.array_equal(raw_data, rounded):
            raise ValueError("Expected raw integer UMI counts in H5AD raw.X")
        data = rounded.astype(np.int64, copy=False)
        local_indptr = indptr[first_row : last_row + 1] - first_value
        nonempty = np.flatnonzero(np.diff(local_indptr))
        if len(nonempty):
            totals[first_row + nonempty] = np.add.reduceat(data, local_indptr[nonempty])
        if len(indices):
            row_for_value = np.searchsorted(local_indptr[1:], np.arange(len(indices)), side="right")
            for target_index, target_values in targets.items():
                positions = np.flatnonzero(indices == target_index)
                if len(positions):
                    np.add.at(
                        target_values[first_row:last_row], row_for_value[positions], data[positions]
                    )
    return totals, targets


def parse_gse235760_h5ad(path: Path) -> pd.DataFrame:
    """Extract author-QC CAR-positive cells and raw target counts from GSE235760."""
    with h5py.File(path, "r") as handle:
        required = ("obs", "raw/X", "raw/var")
        missing = [name for name in required if name not in handle]
        if missing:
            raise ValueError(f"GSE235760 H5AD is missing: {missing}")
        obs = handle["obs"]
        var = handle["raw/var"]
        cell_ids = read_h5ad_column(obs, "_index")
        patients = read_h5ad_column(obs, "donor_id")
        stages = read_h5ad_column(obs, "Timepoint")
        transduction = read_h5ad_column(obs, "Transduction")
        sample_ids = read_h5ad_column(obs, "Sample_id")
        gene_symbols = read_h5ad_column(var, "feature_name")
        target_genes = ("CXCR6", *MARKER_GENES)
        gene_indices: dict[str, int] = {}
        for gene in target_genes:
            matches = np.flatnonzero(gene_symbols == gene)
            if len(matches) != 1:
                raise ValueError(f"Expected one raw feature for {gene}; observed {len(matches)}")
            gene_indices[gene] = int(matches[0])
        totals, target_counts = _csr_target_counts(handle["raw/X"], gene_indices.values())
    eligible = transduction == "CAR+"
    payload: dict[str, object] = {
        "dataset": "GSE235760",
        "patient_id": patients[eligible],
        "product": "varni-cel",
        "stage": stages[eligible],
        "biological_sample_id": sample_ids[eligible],
        "cell_id": cell_ids[eligible],
        "total_umi": totals[eligible],
        "cxcr6_umi": target_counts[gene_indices["CXCR6"]][eligible],
    }
    for gene in MARKER_GENES:
        payload[f"gene_{gene}_umi"] = target_counts[gene_indices[gene]][eligible]
    cells = pd.DataFrame(payload)
    if len(cells) != 18_978:
        raise ValueError(f"Expected 18,978 deposited CAR-positive cells; observed {len(cells):,}")
    if set(cells["stage"]) != {"IP", "Peak"}:
        raise ValueError(f"Unexpected GSE235760 time points: {sorted(set(cells['stage']))}")
    return validate_cells(cells)


def _open_dense_member(archive: tarfile.TarFile, member_name: str) -> io.TextIOWrapper:
    member = archive.getmember(member_name)
    if not member.isfile() or Path(member.name).name != member.name:
        raise ValueError(f"Unsafe or non-file archive member: {member_name}")
    source = archive.extractfile(member)
    if source is None:
        raise ValueError(f"Cannot read archive member: {member_name}")
    compressed = gzip.GzipFile(fileobj=source, mode="rb")
    return io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def parse_dense_gene_by_cell(
    handle: io.TextIOBase,
    *,
    dataset: str,
    patient_id: str,
    product: str,
    stage: str,
    biological_sample_id: str,
    eligibility_gene: str | None = None,
) -> pd.DataFrame:
    """Parse a tab- or whitespace-delimited dense integer gene-by-cell table.

    GEO processed matrices are commonly emitted by R ``write.table``.  In that
    representation the header contains only quoted cell barcodes, whereas each
    data row starts with a quoted gene symbol.  Small standardized fixtures may
    instead include an explicit row-name header (for example, ``gene``).  The
    first data row is therefore used to determine which of the two layouts was
    supplied without relying on a dataset-specific barcode pattern.
    """
    header = handle.readline().rstrip("\r\n").split()
    if not header:
        raise ValueError("Dense matrix header has no cell columns")

    first_data_line = handle.readline()
    if not first_data_line:
        raise ValueError("Dense matrix contains no gene rows")

    def parse_row(
        line: str, line_number: int, expected_values: int | None = None
    ) -> tuple[str, np.ndarray]:
        fields = line.rstrip("\r\n").split(maxsplit=1)
        if len(fields) != 2:
            raise ValueError(f"Malformed dense matrix row {line_number}")
        gene = fields[0].strip('"')
        count_text = fields[1].replace('"', "")
        if re.search(r"[^0-9 \t]", count_text):
            raise ValueError(f"Non-integer counts on dense matrix row {line_number}")
        values = np.fromstring(count_text, dtype=np.int64, sep=" ")
        if expected_values is not None and len(values) != expected_values:
            raise ValueError(f"Invalid counts on dense matrix row {line_number}")
        return gene, values

    first_gene, first_values = parse_row(first_data_line, 2)
    if len(header) == len(first_values):
        barcode_fields = header
    elif len(header) == len(first_values) + 1:
        barcode_fields = header[1:]
    else:
        raise ValueError("Dense matrix header and first gene row have incompatible column counts")
    cell_ids = np.asarray([value.strip('"') for value in barcode_fields], dtype=object)
    if not len(cell_ids):
        raise ValueError("Dense matrix header has no cell columns")
    if len(set(cell_ids)) != len(cell_ids):
        raise ValueError("Dense matrix cell barcodes are not unique")
    requested = {"CXCR6", *MARKER_GENES}
    if eligibility_gene:
        requested.add(eligibility_gene)
    target_counts = {gene: np.zeros(len(cell_ids), dtype=np.int64) for gene in requested}
    total_umi = np.zeros(len(cell_ids), dtype=np.int64)
    seen: set[str] = set()

    def consume_row(gene: str, values: np.ndarray) -> None:
        total_umi[:] += values
        if gene in requested:
            target_counts[gene] += values
            seen.add(gene)

    consume_row(first_gene, first_values)
    for line_number, line in enumerate(handle, start=3):
        gene, values = parse_row(line, line_number, len(cell_ids))
        consume_row(gene, values)
    missing = requested.difference(seen)
    if missing:
        raise ValueError(f"Dense matrix is missing target genes: {sorted(missing)}")
    eligible = (
        np.ones(len(cell_ids), dtype=bool)
        if eligibility_gene is None
        else target_counts[eligibility_gene] > 0
    )
    payload: dict[str, object] = {
        "dataset": dataset,
        "patient_id": patient_id,
        "product": product,
        "stage": stage,
        "biological_sample_id": biological_sample_id,
        "cell_id": cell_ids[eligible],
        "total_umi": total_umi[eligible],
        "cxcr6_umi": target_counts["CXCR6"][eligible],
    }
    for gene in MARKER_GENES:
        payload[f"gene_{gene}_umi"] = target_counts[gene][eligible]
    return validate_cells(pd.DataFrame(payload))


def parse_gse273170_archive(archive_path: Path, crosswalk_path: Path) -> pd.DataFrame:
    """Extract CAR-transcript-positive cells from the GSE273170 dense matrices."""
    crosswalk = pd.read_csv(crosswalk_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_file", "patient_id", "stage", "product", "sha256"}
    missing = required.difference(crosswalk.columns)
    if missing:
        raise ValueError(f"GSE273170 crosswalk is missing columns: {sorted(missing)}")
    frames: list[pd.DataFrame] = []
    with tarfile.open(archive_path, "r") as archive:
        members = {member.name: member for member in archive.getmembers()}
        for row in crosswalk.to_dict("records"):
            if row["source_file"] not in members:
                raise ValueError(f"Archive is missing {row['source_file']}")
            source = archive.extractfile(members[row["source_file"]])
            if source is None:
                raise ValueError(f"Cannot read {row['source_file']}")
            compressed_bytes = source.read()
            import hashlib

            observed_sha = hashlib.sha256(compressed_bytes).hexdigest()
            if observed_sha != row["sha256"]:
                raise ValueError(f"SHA-256 mismatch for {row['source_file']}")
            with io.TextIOWrapper(
                gzip.GzipFile(fileobj=io.BytesIO(compressed_bytes)), encoding="utf-8", newline=""
            ) as text_handle:
                frames.append(
                    parse_dense_gene_by_cell(
                        text_handle,
                        dataset="GSE273170",
                        patient_id=row["patient_id"],
                        product=row["product"],
                        stage=row["stage"],
                        biological_sample_id=row["source_file"],
                        eligibility_gene="CAR",
                    )
                )
    return validate_cells(pd.concat(frames, ignore_index=True))


def parse_gse162975_matrix(matrix_path: Path, crosswalk_path: Path) -> pd.DataFrame:
    """Parse the sequence-validated GSE162975 dense UMI matrix and GEO crosswalk."""
    crosswalk = pd.read_csv(crosswalk_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"technical_record", "patient_id", "deposited_stage"}
    missing = required.difference(crosswalk.columns)
    if missing:
        raise ValueError(f"GSE162975 crosswalk is missing columns: {sorted(missing)}")
    if crosswalk["technical_record"].duplicated().any():
        raise ValueError("GSE162975 technical records must be unique")
    metadata = crosswalk.set_index("technical_record")
    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        if len(header) != 7_579:
            raise ValueError(f"Expected 7,578 GSE162975 cells; observed {len(header) - 1:,}")
        cell_ids = np.asarray([value.strip('"') for value in header[1:]], dtype=object)
        requested = {"CXCR6", *MARKER_GENES}
        target_counts = {gene: np.zeros(len(cell_ids), dtype=np.int64) for gene in requested}
        total_umi = np.zeros(len(cell_ids), dtype=np.int64)
        seen: set[str] = set()
        gene_rows = 0
        for line_number, line in enumerate(handle, start=2):
            gene, separator, remainder = line.rstrip("\r\n").partition(",")
            if not separator:
                raise ValueError(f"Malformed GSE162975 matrix row {line_number}")
            gene = gene.strip('"')
            values = np.fromstring(remainder.replace('"', ""), dtype=np.int64, sep=",")
            if len(values) != len(cell_ids) or (values < 0).any():
                raise ValueError(f"Invalid GSE162975 counts on row {line_number}")
            total_umi += values
            gene_rows += 1
            if gene in requested:
                target_counts[gene] += values
                seen.add(gene)
    missing_genes = requested.difference(seen)
    if missing_genes:
        raise ValueError(f"GSE162975 matrix is missing target genes: {sorted(missing_genes)}")
    if gene_rows != 22_397:
        raise ValueError(f"Expected 22,397 GSE162975 genes; observed {gene_rows:,}")
    if int(total_umi.sum()) != 577_215_764:
        raise ValueError("GSE162975 total UMI count does not match the audited matrix")
    technical_records = np.asarray([cell.rsplit("_sc", 1)[0] for cell in cell_ids])
    missing_records = sorted(set(technical_records).difference(metadata.index))
    if missing_records:
        raise ValueError(f"Unmapped GSE162975 technical records: {missing_records[:5]}")
    mapped = metadata.loc[technical_records]
    payload: dict[str, object] = {
        "dataset": "GSE162975",
        "patient_id": mapped["patient_id"].to_numpy(),
        "product": "CD19/CD22 CAR-T cocktail",
        "stage": mapped["deposited_stage"].to_numpy(),
        "biological_sample_id": technical_records,
        "cell_id": cell_ids,
        "total_umi": total_umi,
        "cxcr6_umi": target_counts["CXCR6"],
    }
    for gene in MARKER_GENES:
        payload[f"gene_{gene}_umi"] = target_counts[gene]
    return validate_cells(pd.DataFrame(payload))


def parse_10x_mtx(
    matrix: BinaryIO,
    features: BinaryIO,
    barcodes: BinaryIO,
    *,
    dataset: str,
    patient_id: str,
    product: str,
    stage: str,
    biological_sample_id: str,
    allowed_barcodes: set[str] | None = None,
    eligibility_gene: str | None = None,
    expected_matrix_cells: int | None = None,
) -> pd.DataFrame:
    """Read a 10x Matrix Market bundle into the locked cell interchange schema."""
    feature_table = pd.read_csv(features, sep="\t", header=None, compression="gzip", dtype=str)
    if feature_table.shape[1] < 2:
        raise ValueError("10x feature table must contain feature ID and symbol")
    symbols = feature_table.iloc[:, 1].astype(str).to_numpy()
    barcode_values = pd.read_csv(
        barcodes, sep="\t", header=None, compression="gzip", dtype=str
    ).iloc[:, 0]
    if expected_matrix_cells is not None and len(barcode_values) != expected_matrix_cells:
        raise ValueError(
            f"Expected {expected_matrix_cells:,} 10x barcodes; observed {len(barcode_values):,}"
        )
    counts = csc_matrix(mmread(gzip.GzipFile(fileobj=matrix)))
    if counts.shape != (len(symbols), len(barcode_values)):
        raise ValueError("10x matrix dimensions do not match features and barcodes")
    if counts.data.size and (
        (counts.data < 0).any() or not np.array_equal(counts.data, np.rint(counts.data))
    ):
        raise ValueError("10x matrix must contain non-negative integer counts")
    select = np.ones(len(barcode_values), dtype=bool)
    if allowed_barcodes is not None:
        select &= barcode_values.isin(allowed_barcodes).to_numpy()

    def gene_counts(gene: str) -> np.ndarray:
        matches = np.flatnonzero(symbols == gene)
        if not len(matches):
            raise ValueError(f"10x features do not contain {gene}")
        return np.asarray(counts[matches, :].sum(axis=0)).ravel().astype(np.int64)

    if eligibility_gene:
        select &= gene_counts(eligibility_gene) > 0
    payload: dict[str, object] = {
        "dataset": dataset,
        "patient_id": patient_id,
        "product": product,
        "stage": stage,
        "biological_sample_id": biological_sample_id,
        "cell_id": barcode_values.to_numpy()[select],
        "total_umi": np.asarray(counts.sum(axis=0)).ravel().astype(np.int64)[select],
        "cxcr6_umi": gene_counts("CXCR6")[select],
    }
    for gene in MARKER_GENES:
        payload[f"gene_{gene}_umi"] = gene_counts(gene)[select]
    return validate_cells(pd.DataFrame(payload))


def _load_gse197268_annotations(obs_path: Path, subtype_path: Path) -> pd.DataFrame:
    """Join the two author-provided cell annotations without inventing missing labels."""
    obs = pd.read_csv(obs_path, index_col=0, low_memory=False)
    subtype = pd.read_csv(subtype_path, sep="\t", index_col=0, low_memory=False)
    required_obs = {"barcode", "timepoint", "generic", "CAR", "cell_type"}
    missing_obs = required_obs.difference(obs.columns)
    if missing_obs:
        raise ValueError(f"GSE197268 global metadata is missing columns: {sorted(missing_obs)}")
    if "subtype" not in subtype.columns:
        raise ValueError("GSE197268 subtype metadata is missing the subtype column")
    if not obs.index.is_unique or not subtype.index.is_unique:
        raise ValueError("GSE197268 author annotation cell identifiers must be unique")
    car_strings = obs["CAR"].astype(str).str.lower()
    if not car_strings.isin({"true", "false"}).all():
        raise ValueError("GSE197268 CAR annotations must be Boolean")
    obs = obs.copy()
    obs["CAR"] = car_strings.eq("true")
    annotations = obs.join(subtype[["subtype"]], how="left", validate="one_to_one")
    annotations.index = annotations.index.astype(str)
    annotations.index.name = "cell_id"
    return annotations


def _summarize_gse197268_screening(annotations: pd.DataFrame, min_cells: int = 25) -> pd.DataFrame:
    """Count the outcome-blind GSE197268 denominator screen for every paired patient."""
    paired_stages = ("Infusion", "D7-CART")
    sample_counts = (
        annotations.loc[annotations["timepoint"].isin(paired_stages)]
        .groupby(["barcode", "generic", "timepoint"], sort=True, observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=paired_stages, fill_value=0)
    )
    sampling_frame = sample_counts.loc[
        sample_counts["Infusion"].gt(0) & sample_counts["D7-CART"].gt(0)
    ].index
    eligible_counts = (
        annotations.loc[
            annotations["CAR"]
            & annotations["subtype"].eq("CD8 T")
            & annotations["timepoint"].isin(paired_stages)
        ]
        .groupby(["barcode", "generic", "timepoint"], sort=True, observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=sampling_frame, columns=paired_stages, fill_value=0)
        .reset_index()
    )
    eligible_counts["patient_id"] = "P" + eligible_counts["barcode"].str.rsplit("-", n=1).str[-1]
    included = eligible_counts[list(paired_stages)].ge(min_cells).all(axis=1)
    result = pd.DataFrame(
        {
            "patient_id": eligible_counts["patient_id"],
            "author_patient_label": eligible_counts["barcode"],
            "product": eligible_counts["generic"],
            "ip_eligible_car_positive_cd8_t_cells": eligible_counts["Infusion"].astype(int),
            "d7_cart_eligible_car_positive_cd8_t_cells": eligible_counts["D7-CART"].astype(int),
            "included_at_min25": included.astype(int),
            "reason": np.where(included, "included", "below_minimum_eligible_cells"),
        }
    )
    return result.sort_values("patient_id", kind="mergesort").reset_index(drop=True)


def build_gse197268_screening(
    obs_path: Path, subtype_path: Path, min_cells: int = 25
) -> pd.DataFrame:
    """Rebuild the complete 21-patient GSE197268 denominator-screen audit."""
    screening = _summarize_gse197268_screening(
        _load_gse197268_annotations(obs_path, subtype_path), min_cells=min_cells
    )
    if min_cells == 25 and (len(screening) != 21 or screening["included_at_min25"].sum() != 9):
        raise ValueError("Expected 9 of 21 paired GSE197268 patients to pass the 25-cell rule")
    return screening


def _gse197268_archive_members(
    archive: tarfile.TarFile,
) -> tuple[tarfile.ExFileObject, tarfile.ExFileObject, tarfile.ExFileObject]:
    """Return matrix, features, and barcodes streams from one safe 10x archive."""

    def select(suffix: str) -> tarfile.ExFileObject:
        matches = []
        for member in archive.getmembers():
            parts = PurePosixPath(member.name).parts
            if member.name.startswith("/") or ".." in parts:
                raise ValueError(f"Unsafe GSE197268 archive member: {member.name}")
            if member.isfile() and member.name.endswith(suffix):
                matches.append(member)
        if len(matches) != 1:
            raise ValueError(f"Expected one {suffix} member; observed {len(matches)}")
        stream = archive.extractfile(matches[0])
        if stream is None:
            raise ValueError(f"Cannot read GSE197268 archive member: {matches[0].name}")
        return stream

    return select("matrix.mtx.gz"), select("features.tsv.gz"), select("barcodes.tsv.gz")


def parse_gse197268_archives(
    archive_dir: Path,
    obs_path: Path,
    subtype_path: Path,
    crosswalk_path: Path,
) -> pd.DataFrame:
    """Recover author-QC CAR-positive CD8 T cells from the paired GSE197268 bundles."""
    crosswalk = pd.read_csv(crosswalk_path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "source_file",
        "geo_sample",
        "patient_number",
        "patient_id",
        "author_patient_label",
        "product",
        "canonical_stage",
        "author_timepoint",
        "raw_matrix_cells",
        "author_qc_matched_cells",
        "expected_eligible_cells",
        "expected_cxcr6_positive_cells",
        "bytes",
        "sha256",
    }
    missing = required.difference(crosswalk.columns)
    if missing:
        raise ValueError(f"GSE197268 sample crosswalk is missing columns: {sorted(missing)}")
    if len(crosswalk) != 18 or crosswalk["source_file"].duplicated().any():
        raise ValueError("GSE197268 sample crosswalk must identify exactly 18 unique bundles")
    stages_per_patient = crosswalk.groupby("patient_id")["canonical_stage"].agg(set)
    if (
        len(stages_per_patient) != 9
        or not stages_per_patient.map(lambda values: values == {"IP", "D7-CART"}).all()
    ):
        raise ValueError("GSE197268 crosswalk must contain IP and D7-CART for nine patients")
    annotations = _load_gse197268_annotations(obs_path, subtype_path)
    frames: list[pd.DataFrame] = []
    for row in crosswalk.to_dict("records"):
        archive_path = archive_dir / row["source_file"]
        expected_bytes = int(row["bytes"])
        if not archive_path.is_file() or archive_path.stat().st_size != expected_bytes:
            raise ValueError(f"GSE197268 source size mismatch: {row['source_file']}")
        if sha256(archive_path) != row["sha256"]:
            raise ValueError(f"GSE197268 source checksum mismatch: {row['source_file']}")
        sample_annotations = annotations.loc[
            annotations["barcode"].eq(row["author_patient_label"])
            & annotations["timepoint"].eq(row["author_timepoint"])
            & annotations["generic"].eq(row["product"])
        ]
        expected_author_qc = int(row["author_qc_matched_cells"])
        if len(sample_annotations) != expected_author_qc:
            raise ValueError(
                f"GSE197268 author-QC sample count mismatch for {row['source_file']}: "
                f"{len(sample_annotations)} != {expected_author_qc}"
            )
        eligible_annotations = sample_annotations.loc[
            sample_annotations["CAR"] & sample_annotations["subtype"].eq("CD8 T")
        ]
        allowed_barcodes = set(eligible_annotations.index)
        expected_cells = int(row["expected_eligible_cells"])
        if len(allowed_barcodes) != expected_cells:
            raise ValueError(
                f"GSE197268 author annotation count mismatch for {row['source_file']}: "
                f"{len(allowed_barcodes)} != {expected_cells}"
            )
        with tarfile.open(archive_path, "r:gz") as archive:
            matrix, features, barcodes = _gse197268_archive_members(archive)
            cells = parse_10x_mtx(
                matrix,
                features,
                barcodes,
                dataset="GSE197268",
                patient_id=row["patient_id"],
                product=row["product"],
                stage=row["canonical_stage"],
                biological_sample_id=row["geo_sample"],
                allowed_barcodes=allowed_barcodes,
                expected_matrix_cells=int(row["raw_matrix_cells"]),
            )
        observed_positive = int(cells["cxcr6_umi"].gt(0).sum())
        expected_positive = int(row["expected_cxcr6_positive_cells"])
        if len(cells) != expected_cells or observed_positive != expected_positive:
            raise ValueError(
                f"GSE197268 raw-count audit mismatch for {row['source_file']}: "
                f"cells={len(cells)}, CXCR6-positive={observed_positive}"
            )
        frames.append(cells)
    result = validate_cells(pd.concat(frames, ignore_index=True, sort=False))
    expected_total = crosswalk["expected_eligible_cells"].astype(int).sum()
    if len(result) != expected_total:
        raise ValueError(
            f"Expected {expected_total:,} eligible GSE197268 cells; observed {len(result):,}"
        )
    return result


def read_standard_cells(path: Path) -> pd.DataFrame:
    """Read the documented lossless cell interchange table."""
    if path.suffix == ".parquet":
        frame = pd.read_parquet(path)
    else:
        frame = pd.read_csv(path, sep="\t")
    return validate_cells(frame)


def write_standard_cells(cells: pd.DataFrame, path: Path) -> None:
    """Write a deterministic gzip TSV interchange file."""
    cells = validate_cells(cells)
    path.parent.mkdir(parents=True, exist_ok=True)
    compression: str | dict[str, object] | None = None
    if path.suffix == ".gz":
        compression = {"method": "gzip", "mtime": 0}
    cells.to_csv(
        path,
        sep="\t",
        index=False,
        columns=[*CELL_COLUMNS, *[column for column in cells if column.startswith("gene_")]],
        compression=compression,
        lineterminator="\n",
    )
