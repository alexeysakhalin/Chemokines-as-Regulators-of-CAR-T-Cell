# Public single-cell and spatial reanalysis

This module reproduces the exploratory secondary analysis reported in Figure 6 of the revised review. It uses two independent public datasets and keeps their biological units and clinical contexts separate:

- **GSE125881**: 62,167 author-filtered CD8-positive CAR-T cells from four patients, measured in the infusion product and at early (days 12–21), contraction/late (days 28–38), and very-late (days 83–112) post-infusion phases.
- **GSE269379**: author-curated cerebrospinal-fluid single-cell metadata from five patients sampled during ICANS and four idiopathic-intracranial-hypertension comparators, plus one Visium HD section from a fatal ICANS case.

The analysis is deliberately narrow. For GSE125881, it summarizes transcript detection across defined genes and marker sets; for GSE269379, it uses the authors' CAR, CXCR6, and cell-type annotations. It does not re-cluster either dataset, infer a new cell atlas, or pool cells across patients as independent replicates. The GSE269379 spatial specimen is a single case. Its CXCL16 map is descriptive, and sparse CXCR6 transcript detection cannot establish transcript-level proximity, a functional protein gradient, direction of migration, or the membrane-bound versus soluble CXCL16 proteoform.

## Reproduce on Linux

From this directory:

```bash
mamba env create --file environment.yml
mamba activate cart-public-reanalysis
make fetch
make analysis
make verify
make test
```

The complete Snakemake route is:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 --rerun-incomplete
```

`make fetch` downloads the exact public files recorded in `config/sources.tsv`, verifies byte sizes and SHA-256 digests, and extracts only the spatial archive members required by the analysis. Third-party matrices and images remain below `data/raw/` and are excluded from Git. The committed `results/frozen/` directory contains only derived aggregate tables, the manuscript figure, and a run manifest so that a local run can be compared with the reviewed result.

`config/gse125881_sample_crosswalk.tsv` records the exact patient-specific collection days and GEO accessions used to interpret the four author-defined longitudinal phases. Its sample-day provenance is the checksum-pinned GSE125881 family SOFT record listed in `config/sources.tsv`.

## Statistical design

For GSE269379, the reported between-group single-cell endpoints are calculated within donor. The independent unit is the donor, not an individual cell. Acute-ICANS and comparator donor summaries are compared with exact two-sided Mann–Whitney tests; Benjamini–Hochberg values are reported for the two defined endpoints. The within-donor CAR-positive versus CAR-negative summary is paired and exploratory. The single longitudinal ICANS donor is descriptive.

For GSE125881, CXCR6-detected fractions and CXCR6 counts per million total UMIs are calculated for each patient and phase. With four selected patients, phase changes are reported as within-patient directions, medians, and observed ranges; no population-level efficacy claim is made.

The GSE125881 marker sets are memory-associated (CCR7, IL7R, LEF1, LTB, MAL, SELL, TCF7), effector-associated (CCL5, GNLY, GZMB, IFNG, NKG7, PRF1), and dysfunction-associated (ENTPD1, HAVCR2, LAG3, PDCD1, TIGIT, TOX). Within each patient and phase, marker detection is the number of nonzero marker–cell pairs divided by the number of genes in the set multiplied by the number of cells. Figure 6B shows the median of this value across the four patients. This is a descriptive detection summary, not a module score or cell-state classification.

For the single GSE269379 Visium HD section, the workflow reports transcript-positive 8-µm bins. A conditional depth-matched randomization checks whether CXCL16-positive bins occur within 24 µm of bins containing at least one of CD14, CD68, LST1, AIF1, TYROBP, or FCER1G. The Monte Carlo tail probability is conditional on bin-level exchangeability within library-size deciles; the null does not preserve spatial autocorrelation or tissue architecture. This is a within-section diagnostic, not patient-level inference, formal biological significance, or cell-type assignment. CXCR6 is displayed only as sparse transcript detection.

## Outputs

- `results/frozen/figure6_public_reanalysis.png`, `.pdf`, and `.svg`
- donor- and patient-level aggregate tables in `results/frozen/tables/`
- `results/frozen/analysis_manifest.json` with input digests and run parameters
- `results/frozen/SHA256SUMS` for the committed outputs

## Source studies

- Sheih A, Voillet V, Hanafi LA, et al. *Nat Commun.* 2020;11:219. doi:10.1038/s41467-019-13880-1. GEO: GSE125881.
- Lu IN, Müller-Miny L, Krekeler C, et al. *Genome Med.* 2025;17:71. doi:10.1186/s13073-025-01498-6. GEO: GSE269379.

The upstream files remain governed by their original repository terms. This repository does not redistribute the source matrices, cell-level metadata, or tissue images.
