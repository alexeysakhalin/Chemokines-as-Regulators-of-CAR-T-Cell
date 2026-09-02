#!/usr/bin/env python3
"""Reproduce the public CAR-T single-cell and Visium HD summaries."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon


T_CELL_STATES = ("GZMAlo CD4", "GZMAhi CD4", "CD8", "Proliferating")
STATE_LABELS = {
    "GZMAlo CD4": "GZMA-low CD4",
    "GZMAhi CD4": "GZMA-high CD4",
    "CD8": "CD8",
    "Proliferating": "Proliferating",
}
STATE_COLORS = {
    "GZMAlo CD4": "#8DB8C7",
    "GZMAhi CD4": "#D9795E",
    "CD8": "#6B6AAE",
    "Proliferating": "#E8B34B",
}
MYELOID_MARKERS = ("CD14", "CD68", "LST1", "AIF1", "TYROBP", "FCER1G")
CHEMOKINE_RECEPTORS = ("CCR5", "CCR7", "CX3CR1", "CXCR3", "CXCR4", "CXCR6")
STATE_MARKER_SETS = {
    "Memory-associated": ("CCR7", "IL7R", "LEF1", "LTB", "MAL", "SELL", "TCF7"),
    "Effector-associated": ("CCL5", "GNLY", "GZMB", "IFNG", "NKG7", "PRF1"),
    "Dysfunction-associated": ("ENTPD1", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TOX"),
}
SPATIAL_INPUT_FILES = (
    "GSM8968967_barcodes.tsv.gz",
    "GSM8968967_features.tsv.gz",
    "GSM8968967_matrix.mtx.gz",
    "GSM8968967_scalefactors_json.json.gz",
    "GSM8968967_tissue_hires_image.png.gz",
    "GSM8968967_tissue_positions.parquet.gz",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    order = np.argsort(values.to_numpy())
    ranked = values.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=values.index)


def analyze_single_cell(metadata_path: Path) -> dict[str, pd.DataFrame]:
    cells = pd.read_csv(metadata_path)
    expected = {"donor_id", "disease", "CAR", "cell_type", "CXCR6"}
    missing = expected.difference(cells.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")

    cells = cells.rename(columns={cells.columns[0]: "cell_barcode"})
    cells["analysis_group"] = np.select(
        [cells["disease"].str.startswith("ICANS_Grade"), cells["disease"].eq("IIH")],
        ["ICANS onset", "IIH comparator"],
        default="Longitudinal ICANS1",
    )
    cells["is_t_cell"] = cells["cell_type"].isin(T_CELL_STATES)
    cells["cxcr6_positive"] = cells["CXCR6"].eq("CXCR6_Pos")
    cells["car_positive"] = cells["CAR"].eq("CAR_Pos")
    cells["proliferating"] = cells["cell_type"].eq("Proliferating")

    cross = cells[cells["analysis_group"].isin(("ICANS onset", "IIH comparator"))]
    donor = (
        cross.groupby(["analysis_group", "donor_id"], observed=True)
        .agg(
            n_cells=("cell_barcode", "size"),
            n_t_cells=("is_t_cell", "sum"),
            cxcr6_positive_all=("cxcr6_positive", "sum"),
            proliferating_cells=("proliferating", "sum"),
        )
        .reset_index()
    )
    donor["cxcr6_fraction_all"] = donor["cxcr6_positive_all"] / donor["n_cells"]
    donor["proliferating_fraction_all"] = donor["proliferating_cells"] / donor["n_cells"]

    t_donor = (
        cross[cross["is_t_cell"]]
        .groupby(["analysis_group", "donor_id"], observed=True)
        .agg(cxcr6_positive_t=("cxcr6_positive", "sum"), n_t_cells_check=("cell_barcode", "size"))
        .reset_index()
    )
    t_donor["cxcr6_fraction_t"] = t_donor["cxcr6_positive_t"] / t_donor["n_t_cells_check"]
    donor = donor.merge(t_donor, on=["analysis_group", "donor_id"], validate="one_to_one")
    if not donor["n_t_cells"].eq(donor["n_t_cells_check"]).all():
        raise AssertionError("T-cell denominators do not agree")
    donor = donor.drop(columns="n_t_cells_check")

    state = (
        cross[cross["is_t_cell"]]
        .groupby(["analysis_group", "donor_id", "cell_type"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    state = state.merge(
        donor[["analysis_group", "donor_id", "n_t_cells"]],
        on=["analysis_group", "donor_id"],
        validate="many_to_one",
    )
    state["fraction_of_t_cells"] = state["n"] / state["n_t_cells"]

    endpoints = []
    for endpoint, column in (
        ("CXCR6-positive fraction among T cells", "cxcr6_fraction_t"),
        ("Proliferating-cell fraction among all cells", "proliferating_fraction_all"),
    ):
        onset = donor.loc[donor["analysis_group"].eq("ICANS onset"), column].to_numpy()
        control = donor.loc[donor["analysis_group"].eq("IIH comparator"), column].to_numpy()
        test = mannwhitneyu(onset, control, alternative="two-sided", method="exact")
        endpoints.append(
            {
                "endpoint": endpoint,
                "n_icans_donors": len(onset),
                "n_comparator_donors": len(control),
                "icans_median": float(np.median(onset)),
                "icans_min": float(np.min(onset)),
                "icans_max": float(np.max(onset)),
                "comparator_median": float(np.median(control)),
                "comparator_min": float(np.min(control)),
                "comparator_max": float(np.max(control)),
                "median_difference": float(np.median(onset) - np.median(control)),
                "mann_whitney_u": float(test.statistic),
                "p_value_exact": float(test.pvalue),
            }
        )
    tests = pd.DataFrame(endpoints)
    tests["q_value_bh"] = benjamini_hochberg(tests["p_value_exact"])

    onset_t = cross[cross["analysis_group"].eq("ICANS onset") & cross["is_t_cell"]]
    paired = (
        onset_t.groupby(["donor_id", "CAR"], observed=True)["cxcr6_positive"]
        .agg([("cxcr6_positive", "sum"), ("n_t_cells", "size"), ("fraction", "mean")])
        .reset_index()
    )
    paired_wide = paired.pivot(index="donor_id", columns="CAR", values="fraction")
    paired_test = wilcoxon(
        paired_wide["CAR_Pos"], paired_wide["CAR_Neg"], alternative="two-sided", method="exact"
    )
    car_test = pd.DataFrame(
        [
            {
                "endpoint": "Within-donor CXCR6-positive T-cell fraction: CAR-positive versus CAR-negative",
                "n_paired_donors": len(paired_wide),
                "median_car_positive": float(paired_wide["CAR_Pos"].median()),
                "median_car_negative": float(paired_wide["CAR_Neg"].median()),
                "median_paired_difference": float(
                    (paired_wide["CAR_Pos"] - paired_wide["CAR_Neg"]).median()
                ),
                "wilcoxon_statistic": float(paired_test.statistic),
                "p_value_exact": float(paired_test.pvalue),
            }
        ]
    )

    timeline_order = ["pre-ICANS", "ICANS_Grade4", "post-ICANS"]
    timeline = []
    for phase in timeline_order:
        sample = cells[(cells["donor_id"].eq("ICANS1")) & cells["disease"].eq(phase)]
        t_sample = sample[sample["is_t_cell"]]
        timeline.append(
            {
                "phase": phase,
                "n_cells": len(sample),
                "n_t_cells": len(t_sample),
                "cxcr6_fraction_t": float(t_sample["cxcr6_positive"].mean()),
                "car_fraction_t": float(t_sample["car_positive"].mean()),
                "proliferating_fraction_t": float(t_sample["proliferating"].mean()),
            }
        )
    timeline = pd.DataFrame(timeline)

    return {
        "cells": cells,
        "donor_summary": donor,
        "state_summary": state,
        "primary_tests": tests,
        "car_summary": paired,
        "car_test": car_test,
        "timeline": timeline,
    }


def analyze_gse125881(matrix_path: Path, metadata_path: Path) -> dict[str, pd.DataFrame]:
    """Summarize CXCR6 counts across the four published longitudinal phases."""
    metadata = pd.read_csv(metadata_path)
    expected = {"Cell", "Patient", "Group", "nUMI"}
    missing = expected.difference(metadata.columns)
    if missing:
        raise ValueError(f"GSE125881 metadata is missing columns: {sorted(missing)}")
    if metadata["Cell"].duplicated().any():
        raise ValueError("GSE125881 cell identifiers are not unique")

    target_genes = set(CHEMOKINE_RECEPTORS)
    for genes in STATE_MARKER_SETS.values():
        target_genes.update(genes)
    with gzip.open(matrix_path, "rt", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        target_counts = {
            row[0].upper(): np.asarray(row[1:], dtype=np.int64)
            for row in reader
            if row and row[0].upper() in target_genes
        }
    missing_genes = target_genes.difference(target_counts)
    if missing_genes:
        raise ValueError(f"GSE125881 matrix is missing target genes: {sorted(missing_genes)}")
    cell_ids = header[1:]
    cxcr6_counts = target_counts["CXCR6"]
    if len(cell_ids) != len(cxcr6_counts) or len(cell_ids) != len(metadata):
        raise ValueError("GSE125881 matrix and metadata cell counts do not match")
    count_series = pd.Series(cxcr6_counts, index=cell_ids, name="cxcr6_count")
    if set(metadata["Cell"]) != set(count_series.index):
        raise ValueError("GSE125881 matrix and metadata cell identifiers do not match")

    cells = metadata.set_index("Cell").loc[cell_ids].join(count_series, how="left")
    cells["cxcr6_detected"] = cells["cxcr6_count"].gt(0)
    summary = (
        cells.groupby(["Patient", "Group"], observed=True)
        .agg(
            n_cells=("cxcr6_detected", "size"),
            cxcr6_positive_cells=("cxcr6_detected", "sum"),
            cxcr6_umis=("cxcr6_count", "sum"),
            total_umis=("nUMI", "sum"),
        )
        .reset_index()
    )
    phase_order = ["IP", "Early", "Late", "Very late"]
    summary["Group"] = pd.Categorical(summary["Group"], phase_order, ordered=True)
    summary = summary.sort_values(["Patient", "Group"]).reset_index(drop=True)
    summary["cxcr6_detected_fraction"] = summary["cxcr6_positive_cells"] / summary["n_cells"]
    summary["cxcr6_counts_per_million_umis"] = summary["cxcr6_umis"] / summary["total_umis"] * 1e6

    detected = summary.pivot(index="Patient", columns="Group", values="cxcr6_detected_fraction")
    cpm = summary.pivot(index="Patient", columns="Group", values="cxcr6_counts_per_million_umis")
    transitions = pd.DataFrame(
        {
            "patient": detected.index,
            "early_minus_ip_detected_percentage_points": (detected["Early"] - detected["IP"]) * 100,
            "very_late_minus_early_detected_percentage_points":
                (detected["Very late"] - detected["Early"]) * 100,
            "early_minus_ip_cpm": cpm["Early"] - cpm["IP"],
            "very_late_minus_early_cpm": cpm["Very late"] - cpm["Early"],
        }
    ).reset_index(drop=True)
    aggregate = pd.DataFrame(
        [
            {
                "n_patients": len(transitions),
                "median_early_minus_ip_detected_percentage_points": float(
                    transitions["early_minus_ip_detected_percentage_points"].median()
                ),
                "min_early_minus_ip_detected_percentage_points": float(
                    transitions["early_minus_ip_detected_percentage_points"].min()
                ),
                "max_early_minus_ip_detected_percentage_points": float(
                    transitions["early_minus_ip_detected_percentage_points"].max()
                ),
                "median_very_late_minus_early_detected_percentage_points": float(
                    transitions["very_late_minus_early_detected_percentage_points"].median()
                ),
                "min_very_late_minus_early_detected_percentage_points": float(
                    transitions["very_late_minus_early_detected_percentage_points"].min()
                ),
                "max_very_late_minus_early_detected_percentage_points": float(
                    transitions["very_late_minus_early_detected_percentage_points"].max()
                ),
                "direction_concordant_early_increase": bool(
                    transitions["early_minus_ip_detected_percentage_points"].gt(0).all()
                ),
                "direction_concordant_very_late_decrease": bool(
                    transitions["very_late_minus_early_detected_percentage_points"].lt(0).all()
                ),
            }
        ]
    )

    receptor_rows = []
    module_rows = []
    for (patient, phase), indices in cells.groupby(["Patient", "Group"], observed=True).indices.items():
        positions = np.asarray(indices, dtype=int)
        for gene in CHEMOKINE_RECEPTORS:
            values = target_counts[gene][positions]
            receptor_rows.append(
                {
                    "patient": patient,
                    "phase": phase,
                    "n_cells": len(positions),
                    "gene": gene,
                    "positive_cells": int(np.sum(values > 0)),
                    "fraction_positive": float(np.mean(values > 0)),
                    "total_umis": int(values.sum()),
                }
            )
        for module, genes in STATE_MARKER_SETS.items():
            detected = np.vstack([target_counts[gene][positions] > 0 for gene in genes])
            module_rows.append(
                {
                    "patient": patient,
                    "phase": phase,
                    "n_cells": len(positions),
                    "module": module,
                    "genes": ";".join(genes),
                    "detected_gene_cell_pairs": int(detected.sum()),
                    "possible_gene_cell_pairs": int(detected.size),
                    "mean_marker_detection_fraction": float(detected.mean()),
                }
            )
    receptor_summary = pd.DataFrame(receptor_rows)
    module_summary = pd.DataFrame(module_rows)
    return {
        "summary": summary,
        "transitions": transitions,
        "aggregate": aggregate,
        "chemokine_receptors": receptor_summary,
        "state_marker_detection": module_summary,
    }


def validate_gse125881_crosswalk(summary: pd.DataFrame, crosswalk_path: Path) -> None:
    crosswalk = pd.read_csv(crosswalk_path, sep="\t")
    required = {"patient", "authors_metadata_group", "day", "n_cells"}
    missing = required.difference(crosswalk.columns)
    if missing:
        raise ValueError(f"GSE125881 crosswalk is missing columns: {sorted(missing)}")
    expected = (
        crosswalk.groupby(["patient", "authors_metadata_group"], as_index=False)["n_cells"]
        .sum()
        .rename(columns={"patient": "Patient", "authors_metadata_group": "Group"})
        .sort_values(["Patient", "Group"])
        .reset_index(drop=True)
    )
    observed = (
        summary[["Patient", "Group", "n_cells"]]
        .assign(Group=lambda frame: frame["Group"].astype("object"))
        .sort_values(["Patient", "Group"])
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(observed, expected, check_dtype=False)


def _read_gzip_parquet(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rb") as handle:
        return pd.read_parquet(io.BytesIO(handle.read()))


def _matched_null(
    distances: np.ndarray,
    observed_mask: np.ndarray,
    strata: np.ndarray,
    radius_um: float,
    n_permutations: int,
    seed: int,
) -> tuple[float, np.ndarray, float]:
    observed = float(np.mean(distances[observed_mask] <= radius_um))
    observed_indices = np.flatnonzero(observed_mask)
    counts = np.bincount(strata[observed_indices], minlength=10)
    candidates = [np.flatnonzero((strata == stratum) & ~observed_mask) for stratum in range(10)]
    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for permutation in range(n_permutations):
        selected = np.concatenate(
            [
                rng.choice(candidates[stratum], size=int(count), replace=False)
                for stratum, count in enumerate(counts)
                if count
            ]
        )
        null[permutation] = np.mean(distances[selected] <= radius_um)
    p_value = float((1 + np.sum(null >= observed)) / (n_permutations + 1))
    return observed, null, p_value


def analyze_spatial(raw_dir: Path, seed: int, n_permutations: int) -> dict[str, object]:
    prefix = "GSM8968967_"
    features = pd.read_csv(
        raw_dir / f"{prefix}features.tsv.gz",
        sep="\t",
        header=None,
        names=["gene_id", "gene", "feature_type"],
    )
    barcodes = pd.read_csv(
        raw_dir / f"{prefix}barcodes.tsv.gz", sep="\t", header=None, names=["barcode"]
    )
    with gzip.open(raw_dir / f"{prefix}matrix.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsr()
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError("Spatial matrix dimensions do not match features and barcodes")

    positions = _read_gzip_parquet(raw_dir / f"{prefix}tissue_positions.parquet.gz")
    positions = barcodes.merge(positions, on="barcode", how="left", validate="one_to_one")
    if positions.isna().any(axis=None) or not positions["in_tissue"].eq(1).all():
        raise ValueError("Spatial matrix contains unmapped or out-of-tissue barcodes")

    with gzip.open(raw_dir / f"{prefix}scalefactors_json.json.gz", "rt") as handle:
        scalefactors = json.load(handle)
    microns_per_pixel = float(scalefactors["microns_per_pixel"])
    coordinates_um = (
        positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy() * microns_per_pixel
    )

    gene_index: dict[str, int] = {}
    for gene in ("CXCL16", "CXCR6", *MYELOID_MARKERS):
        matches = features.index[features["gene"].eq(gene)].tolist()
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one feature for {gene}; found {len(matches)}")
        gene_index[gene] = matches[0]

    counts = {
        gene: np.asarray(matrix[row, :].todense()).ravel() for gene, row in gene_index.items()
    }
    myeloid_score = np.asarray(
        matrix[[gene_index[gene] for gene in MYELOID_MARKERS], :].sum(axis=0)
    ).ravel()
    library_size = np.asarray(matrix.sum(axis=0)).ravel()
    percent_rank = pd.Series(library_size).rank(method="average", pct=True).to_numpy()
    depth_stratum = np.minimum((percent_rank * 10).astype(int), 9)

    cxcl16_mask = counts["CXCL16"] > 0
    cxcr6_mask = counts["CXCR6"] > 0
    myeloid_mask = myeloid_score > 0
    myeloid_tree = cKDTree(coordinates_um[myeloid_mask])
    distance_to_myeloid, _ = myeloid_tree.query(coordinates_um)

    cxcl16_observed, cxcl16_null, cxcl16_p = _matched_null(
        distance_to_myeloid,
        cxcl16_mask,
        depth_stratum,
        radius_um=24.0,
        n_permutations=n_permutations,
        seed=seed,
    )
    marker_rows = []
    for gene, values in counts.items():
        marker_rows.append(
            {
                "marker": gene,
                "total_umis": int(values.sum()),
                "positive_bins": int(np.sum(values > 0)),
                "maximum_bin_count": int(values.max()),
            }
        )
    marker_rows.append(
        {
            "marker": "Myeloid marker union",
            "total_umis": int(myeloid_score.sum()),
            "positive_bins": int(np.sum(myeloid_mask)),
            "maximum_bin_count": int(myeloid_score.max()),
        }
    )
    marker_summary = pd.DataFrame(marker_rows)

    proximity = pd.DataFrame(
        [
            {
                "query": "CXCL16-positive bins",
                "landmark": "myeloid-marker-positive bins",
                "radius_um": 24.0,
                "n_query_bins": int(np.sum(cxcl16_mask)),
                "observed_fraction_within_radius": cxcl16_observed,
                "null_median": float(np.median(cxcl16_null)),
                "null_q025": float(np.quantile(cxcl16_null, 0.025)),
                "null_q975": float(np.quantile(cxcl16_null, 0.975)),
                "conditional_p_value": cxcl16_p,
                "median_nearest_distance_um": float(np.median(distance_to_myeloid[cxcl16_mask])),
            },
        ]
    )

    image_path = raw_dir / f"{prefix}tissue_hires_image.png"
    if not image_path.exists():
        with gzip.open(raw_dir / f"{prefix}tissue_hires_image.png.gz", "rb") as source:
            image = mpimg.imread(io.BytesIO(source.read()), format="png")
    else:
        image = mpimg.imread(image_path)

    return {
        "matrix_shape": matrix.shape,
        "matrix_nnz": int(matrix.nnz),
        "positions": positions,
        "counts": counts,
        "myeloid_score": myeloid_score,
        "marker_summary": marker_summary,
        "proximity": proximity,
        "image": image,
        "scalefactors": scalefactors,
    }


def _panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.10,
        1.10,
        label,
        transform=axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def create_figure(
    single: dict[str, pd.DataFrame],
    longitudinal: dict[str, pd.DataFrame],
    spatial: dict[str, object],
    output: Path,
) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.5,
            "svg.hashsalt": "public-reanalysis-figure6",
        }
    )
    fig = plt.figure(figsize=(7.15, 10.15))
    fig.subplots_adjust(left=0.075, right=0.98, top=0.91, bottom=0.075, wspace=1.45, hspace=0.95)
    grid = fig.add_gridspec(3, 6, height_ratios=(1.0, 1.0, 1.45))

    # A: analysis population
    ax = fig.add_subplot(grid[0, :2])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    boxes = [
        (
            0.05,
            0.68,
            0.90,
            0.25,
            "GSE125881 longitudinal scRNA-seq\n62,167 CD8+ CAR-T cells\n4 patients × 4 phases",
            "#FFF2DF",
        ),
        (
            0.05,
            0.38,
            0.90,
            0.24,
            "GSE269379 CSF scRNA-seq\n8,661 author-curated cells\n5 ICANS donors + 4 IIH donors",
            "#E8F4F2",
        ),
        (
            0.05,
            0.08,
            0.90,
            0.24,
            "GSE269379 Visium HD\n302,297 filtered 8-µm bins\n1 fatal ICANS case",
            "#F3EAF7",
        ),
    ]
    for x, y, width, height, text_value, color in boxes:
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (x, y), width, height, boxstyle="round,pad=0.015", facecolor=color, edgecolor="#55606A"
            )
        )
        ax.text(x + width / 2, y + height / 2, text_value, ha="center", va="center", fontsize=7.2)
    ax.set_title("Public inputs", pad=4)
    _panel_label(ax, "A")

    # B: longitudinal marker-detection summaries
    ax = fig.add_subplot(grid[0, 2:])
    phase_order = ["IP", "Early", "Late", "Very late"]
    module_order = list(STATE_MARKER_SETS)
    module = longitudinal["state_marker_detection"].copy()
    heatmap = (
        module.groupby(["module", "phase"], observed=True)["mean_marker_detection_fraction"]
        .median()
        .unstack()
        .reindex(index=module_order, columns=phase_order)
    )
    sns.heatmap(
        heatmap * 100,
        annot=True,
        fmt=".1f",
        cmap="YlGnBu",
        vmin=0,
        vmax=100,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Median marker detection (%)", "shrink": 0.75},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels(["Product", "Early\nD12–21", "Late\nD28–38", "Very late\nD83–112"], rotation=0)
    ax.set_yticklabels(["Memory", "Effector", "Dysfunction"], rotation=0)
    ax.set_title("GSE125881 state-marker detection (median of 4 patients)")
    ax.text(
        0.5,
        -0.24,
        "Descriptive marker detection, not cell-state classification",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=6.6,
    )
    _panel_label(ax, "B")

    # C: donor-level CXCR6 fraction
    ax = fig.add_subplot(grid[1, :2])
    donor = single["donor_summary"]
    palette = {"IIH comparator": "#7C8791", "ICANS onset": "#D76C4A"}
    for index, group in enumerate(("IIH comparator", "ICANS onset")):
        group_values = (
            donor.loc[donor["analysis_group"].eq(group)]
            .sort_values("donor_id")["cxcr6_fraction_t"]
            .to_numpy()
        )
        offsets = np.linspace(-0.07, 0.07, len(group_values))
        ax.scatter(
            index + offsets,
            group_values,
            s=28,
            color=palette[group],
            edgecolors="none",
            zorder=2,
        )
    medians = donor.groupby("analysis_group")["cxcr6_fraction_t"].median()
    for index, group in enumerate(("IIH comparator", "ICANS onset")):
        ax.plot([index - 0.24, index + 0.24], [medians[group]] * 2, color="black", lw=1.4)
    p_value = single["primary_tests"].iloc[0]["p_value_exact"]
    ax.text(0.5, 0.96, f"Exact Mann–Whitney P = {p_value:.3f}", transform=ax.transAxes, ha="center", va="top")
    ax.set_xlabel("")
    ax.set_xticks([0, 1], ["IIH comparator", "ICANS onset"])
    ax.set_ylabel("CXCR6-positive T cells")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.set_ylim(0, 0.65)
    ax.set_title("T-cell CXCR6 by donor", pad=7)
    _panel_label(ax, "C")

    # D: paired CAR status comparison
    ax = fig.add_subplot(grid[1, 2:4])
    paired = single["car_summary"].pivot(index="donor_id", columns="CAR", values="fraction")
    for donor_id, row in paired.iterrows():
        ax.plot([0, 1], [row["CAR_Neg"], row["CAR_Pos"]], color="#A1A9AF", lw=0.8, zorder=1)
        ax.scatter([0, 1], [row["CAR_Neg"], row["CAR_Pos"]], color=["#66727B", "#2D8C84"], s=22, zorder=2)
    test = single["car_test"].iloc[0]
    ax.text(0.5, 0.96, f"Exact paired P = {test['p_value_exact']:.3f}", transform=ax.transAxes, ha="center", va="top")
    ax.set_xticks([0, 1], ["CAR-negative", "CAR-positive"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0, 0.65)
    ax.set_ylabel("CXCR6-positive T cells")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.set_title("Within-donor CAR status (n = 5)", pad=7)
    _panel_label(ax, "D")

    # E: longitudinal CAR-T CXCR6 detection
    ax = fig.add_subplot(grid[1, 4:])
    phases = ["IP", "Early", "Late", "Very late"]
    trajectories = longitudinal["summary"].copy()
    x = np.arange(len(phases))
    patient_colors = dict(
        zip(sorted(trajectories["Patient"].unique()), ("#0072B2", "#E69F00", "#009E73", "#CC79A7"))
    )
    for patient, patient_data in trajectories.groupby("Patient", observed=True):
        patient_data = patient_data.set_index("Group").reindex(phases)
        ax.plot(
            x,
            patient_data["cxcr6_detected_fraction"],
            marker="o",
            markersize=3.5,
            lw=1.1,
            label=patient,
            color=patient_colors[patient],
        )
    ax.set_xticks(x, ["Product", "Early\nD12–21", "Late\nD28–38", "Very late\nD83–112"])
    ax.set_ylim(0, 0.34)
    ax.set_ylabel("CXCR6-detected cells")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.set_title("GSE125881 time course (n = 4)", pad=7)
    ax.legend(frameon=False, loc="upper right", ncol=2, handlelength=1.3, columnspacing=0.8)
    _panel_label(ax, "E")

    # F: spatial expression and conditional summaries
    ax = fig.add_subplot(grid[2, :4])
    positions = spatial["positions"]
    counts = spatial["counts"]
    image = spatial["image"]
    ax.imshow(image, origin="upper")
    cxcl16 = counts["CXCL16"] > 0
    cxcr6 = counts["CXCR6"] > 0
    ax.scatter(
        positions.loc[cxcl16, "pxl_col_in_fullres"],
        positions.loc[cxcl16, "pxl_row_in_fullres"],
        s=4,
        c="#F08A24",
        alpha=0.60,
        linewidths=0,
        label=f"CXCL16+ bins (n={int(cxcl16.sum()):,})",
    )
    ax.scatter(
        positions.loc[cxcr6, "pxl_col_in_fullres"],
        positions.loc[cxcr6, "pxl_row_in_fullres"],
        s=28,
        marker="^",
        facecolors="#007C91",
        edgecolors="white",
        linewidths=0.5,
        label=f"CXCR6+ bins (n={int(cxcr6.sum())})",
    )
    ax.set_xlim(400, 1900)
    ax.set_ylim(2280, 760)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Single fatal-ICANS Visium HD section (8-µm bins)")
    ax.legend(frameon=True, facecolor="white", loc="lower right")
    _panel_label(ax, "F")

    ax = fig.add_subplot(grid[2, 4:])
    proximity = spatial["proximity"]
    labels = ["CXCL16 →\nmyeloid markers\n≤24 µm"]
    y = np.array([0.78])
    observed = proximity["observed_fraction_within_radius"].to_numpy()
    null = proximity["null_median"].to_numpy()
    lower = proximity["null_q025"].to_numpy()
    upper = proximity["null_q975"].to_numpy()
    ax.errorbar(
        null,
        y,
        xerr=np.vstack((null - lower, upper - null)),
        fmt="o",
        color="#7C8791",
        ecolor="#A7AFB5",
        capsize=3,
        label="Depth-matched null (95% interval)",
    )
    ax.scatter(observed, y, color="#D76C4A", marker="D", s=28, label="Observed")
    for row_index, ypos in enumerate(y):
        p_value = proximity.iloc[row_index]["conditional_p_value"]
        p_text = "<0.001" if p_value < 0.001 else f"= {p_value:.3f}"
        ax.text(0.03, ypos - 0.10, f"Conditional P {p_text}", fontsize=6.6)
    ax.set_yticks(y, labels)
    ax.set_ylim(-0.05, 1.0)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Fraction within radius")
    ax.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.set_title("Depth-matched spatial checks", pad=7)
    ax.text(
        0.03,
        0.46,
        "◆ observed\n● depth-matched null\n(horizontal bar: 95% interval)\n\n"
        "CXCL16: 976 UMIs / 960 bins\n"
        "CXCR6: 14 UMIs / 14 bins\n\n"
        "CXCR6 detection is too sparse\n"
        "for an independent transcript-level\n"
        "proximity inference.",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=6.5,
    )

    fig.suptitle(
        "Exploratory patient-aware reanalysis of public CAR-T single-cell and spatial data",
        fontsize=11,
        fontweight="bold",
        y=0.975,
    )
    fig.text(
        0.5,
        0.018,
        "Comparator CSF was not obtained from CAR-T-treated patients without ICANS. Spatial results are from one case and are not population-level inference.",
        ha="center",
        va="bottom",
        fontsize=6.6,
    )

    output.mkdir(parents=True, exist_ok=True)
    creator = "CAR-T public-data reproducibility workflow"
    fixed_date = dt.datetime(2026, 9, 2, tzinfo=dt.timezone.utc)
    fig.savefig(
        output / "figure6_public_reanalysis.png",
        dpi=300,
        metadata={"Creator": creator},
    )
    fig.savefig(
        output / "figure6_public_reanalysis.pdf",
        metadata={"Creator": creator, "CreationDate": fixed_date, "ModDate": fixed_date},
    )
    fig.savefig(
        output / "figure6_public_reanalysis.svg",
        metadata={"Creator": creator, "Date": "2026-09-02"},
    )
    plt.close(fig)


def write_outputs(
    single: dict[str, pd.DataFrame],
    longitudinal: dict[str, pd.DataFrame],
    spatial: dict[str, object],
    output_dir: Path,
    input_paths: list[Path],
    seed: int,
    n_permutations: int,
) -> None:
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    table_keys = (
        "donor_summary",
        "state_summary",
        "primary_tests",
        "car_summary",
        "car_test",
        "timeline",
    )
    for key in table_keys:
        single[key].to_csv(table_dir / f"scRNA_{key}.tsv", sep="\t", index=False, float_format="%.8g")
    for key, table in longitudinal.items():
        table.to_csv(
            table_dir / f"gse125881_{key}.tsv", sep="\t", index=False, float_format="%.8g"
        )
    spatial["marker_summary"].to_csv(
        table_dir / "spatial_marker_summary.tsv", sep="\t", index=False, float_format="%.8g"
    )
    spatial["proximity"].to_csv(
        table_dir / "spatial_proximity_summary.tsv", sep="\t", index=False, float_format="%.8g"
    )
    create_figure(single, longitudinal, spatial, output_dir)

    manifest = {
        "datasets": ["GSE125881", "GSE269379"],
        "analysis_scope": "Exploratory secondary analysis",
        "random_seed": seed,
        "spatial_permutations": n_permutations,
        "gse125881_cells_total": int(longitudinal["summary"]["n_cells"].sum()),
        "gse269379_metadata_cells_total": int(len(single["cells"])),
        "spatial_features": int(spatial["matrix_shape"][0]),
        "spatial_bins": int(spatial["matrix_shape"][1]),
        "spatial_nonzero_entries": int(spatial["matrix_nnz"]),
        "inputs": [
            {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in input_paths
        ],
    }
    (output_dir / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    outputs = sorted(
        path
        for path in output_dir.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    )
    checksum_lines = [f"{sha256(path)}  {path.relative_to(output_dir)}" for path in outputs]
    (output_dir / "SHA256SUMS").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--gse125881-matrix", type=Path, required=True)
    parser.add_argument("--gse125881-metadata", type=Path, required=True)
    parser.add_argument("--gse125881-crosswalk", type=Path, required=True)
    parser.add_argument("--gse125881-soft", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument("--spatial-permutations", type=int, default=10000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    single = analyze_single_cell(args.metadata)
    longitudinal = analyze_gse125881(args.gse125881_matrix, args.gse125881_metadata)
    validate_gse125881_crosswalk(longitudinal["summary"], args.gse125881_crosswalk)
    spatial = analyze_spatial(args.raw_dir, args.seed, args.spatial_permutations)
    required_inputs = [
        args.metadata,
        args.gse125881_matrix,
        args.gse125881_metadata,
        args.gse125881_crosswalk,
        args.gse125881_soft,
    ]
    required_inputs.extend(args.raw_dir / name for name in SPATIAL_INPUT_FILES)
    missing_inputs = [path for path in required_inputs if not path.is_file()]
    if missing_inputs:
        raise FileNotFoundError(f"Required input files are missing: {missing_inputs}")
    write_outputs(
        single,
        longitudinal,
        spatial,
        args.output_dir,
        required_inputs,
        args.seed,
        args.spatial_permutations,
    )


if __name__ == "__main__":
    main()
