"""Locked variables used by the longitudinal extension."""

from __future__ import annotations

MARKER_SETS: dict[str, tuple[str, ...]] = {
    "memory_associated": ("CCR7", "IL7R", "LEF1", "LTB", "MAL", "SELL", "TCF7"),
    "effector_associated": ("CCL5", "GNLY", "GZMB", "IFNG", "NKG7", "PRF1"),
    "dysfunction_associated": ("ENTPD1", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TOX"),
}

EARLY_FAMILY_ID = "early_cxcr6_fraction_three_cohort"
PRIMARY_ENDPOINT = "cxcr6_detected_fraction"
DEFAULT_SEED = 20260902

SAMPLE_COLUMNS = (
    "dataset",
    "patient_id",
    "product",
    "stage",
    "biological_sample_id",
    "n_technical_records",
    "eligible_cells",
    "cxcr6_positive_cells",
    "cxcr6_umi",
    "total_umi",
    "cxcr6_fraction",
    "cxcr6_cpm",
)

CELL_COLUMNS = (
    "dataset",
    "patient_id",
    "product",
    "stage",
    "biological_sample_id",
    "cell_id",
    "total_umi",
    "cxcr6_umi",
)

MARKER_GENES = tuple(dict.fromkeys(gene for genes in MARKER_SETS.values() for gene in genes))
MARKER_COUNT_COLUMNS = tuple(f"gene_{gene}_umi" for gene in MARKER_GENES)
