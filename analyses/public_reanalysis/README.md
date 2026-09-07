# Public single-cell and spatial reanalysis

This module contains two related but analytically separate components:

1. the original descriptive GSE125881 longitudinal and GSE269379 ICANS/spatial reanalysis;
2. a patient-level longitudinal extension using independently audited public cohorts GSE197268, GSE235760, GSE162975, and GSE273170.

Cells are never treated as independent biological replicates. Patient-level values are constructed within each cohort, and cohorts are reported separately without pooled cell-level testing or cross-study batch correction.

## Main longitudinal result

The public extension does not support a universal early increase in the fraction of CAR-T cells with detectable CXCR6 transcript. Early cohort medians were positive in GSE197268 and GSE235760 but slightly negative in GSE162975, and none of the three exact patient-level sign tests was significant after Holm adjustment. The later GSE162975 and GSE273170 summaries were directionally lower but small and exploratory. These observations are consistent with temporal and study-context heterogeneity, not a causal trafficking mechanism.

The finalized paired results are:

| Cohort | Contrast | Eligible patients | Median change | Direction | Exact P | Holm P |
|---|---|---:|---:|---:|---:|---:|
| GSE197268 | D7-CART − infusion product | 9 | +23.63 percentage points | 7/9 increased | 0.1797 | 0.5391 |
| GSE235760 | expansion peak − infusion product | 5 | +1.36 percentage points | 4/5 increased | 0.3750 | 0.7500 |
| GSE162975 | T1 − T0 | 10 | −0.93 percentage points | 4/10 increased | 0.7539 | 0.7539 |
| GSE162975 | first T3/T4.5 − T1 | 7 | −15.70 percentage points | 6/7 decreased | 0.1250 | not in early family |
| GSE162975 | first T6-or-later − T1 | 5 | −1.60 percentage points | 3/5 decreased | 1.0000 | not in early family |
| GSE273170 | D14 − D7 | 4 | −14.70 percentage points | 3/4 decreased | 0.6250 | not in early family |

The first three rows form the three-test early comparison family in the finalized exploratory workflow. This wording is deliberate: the workflow was fixed after feasibility checks and is not a prospective preregistration.

## Why the paired n is small

The limiting quantity is not the number of deposited cells or the advertised study cohort size. A patient must have both required biological stages, recoverable CAR-T identity under a cohort-specific public definition, the locked cell stratum, and the minimum eligible-cell denominator at both stages. Missing stages, unsorted aliquots, retreatment courses, sparse CAR-transcript detection, and incomplete public annotations reduce the eligible patient count before inference.

For GSE197268 specifically, 21 patients had author-QC metadata at both infusion product and D7-CART. The outcome-blind denominator screen retained nine with at least 25 author-QC CAR-positive CD8 T cells at each stage and excluded 12 below that threshold. The repository records all 21 decisions; only the 18 raw bundles belonging to the nine eligible patients are fetched for CXCR6 reconstruction.

GSE290722 was specifically audited and excluded from the declared longitudinal CAR-T endpoint. Its repeated-timepoint public metadata do not identify CAR-T cells, while the relevant publisher source-data sheet contains only 128 week-4 cell records from 13 patients. It therefore contributes no infusion-product-to-post-infusion CAR-T pair under the locked definition. The workbook identity, checksum, sheet, and record counts are documented in [`docs/longitudinal_extension/GSE290722_SCREENING.md`](docs/longitudinal_extension/GSE290722_SCREENING.md). The search was targeted rather than an exhaustive meta-analysis; additional candidates require separate CAR-identity and timepoint harmonization.

## Reproduce the longitudinal extension on Linux

From `analyses/public_reanalysis`:

```bash
mamba env create --file environment.yml
mamba activate cart-public-reanalysis
make fetch-longitudinal
make longitudinal
make verify-longitudinal
make test
```

The workflow invokes four cohort-specific raw adapters:

- GSE197268: two checksum-pinned author annotation files and 18 checksum-pinned GEO 10x integer-count bundles for the nine patients passing the outcome-blind 25-cell rule; eligibility is author-QC `CAR == True` and author subtype `CD8 T`;
- GSE235760: the checksum-pinned CELLxGENE H5AD, with integer UMI counts from `raw.X` and `Transduction == "CAR+"`;
- GSE162975: the checksum-pinned dense UMI matrix and the committed 196-record GEO metadata crosswalk;
- GSE273170: the checksum-pinned GEO archive, with a cell called CAR-positive when the custom `CAR` feature has at least one UMI.

`make longitudinal` and the Snakemake workflow call all four adapters directly. The generic analysis layer validates raw integer counts, unique sample-cell identities, patient-stage aggregation, denominator rules, and one contribution per patient per contrast.

`results/longitudinal_extension/frozen/` contains only derived patient-level or cohort-level tables, an analysis report, source/code manifests, and exact audit tables. Third-party cell-level matrices are downloaded to `data/raw/` and are not redistributed.

`make verify-longitudinal` rejects a stale freeze: it checks the complete frozen-file inventory, current working-tree code, source-manifest and audit-table digests, content-addressed input provenance, and all nine declared canonical outputs byte for byte. Raw input keys contain SHA-256 digests rather than machine-specific paths, so an identical run from another Linux checkout or download cache produces the same analysis manifest.

## Original GSE125881/GSE269379 component

- **GSE125881:** 62,167 author-filtered CD8-positive CAR-T cells from four patients across infusion product, early, late, and very-late phases.
- **GSE269379:** author-curated CSF metadata from five ICANS patients and four idiopathic-intracranial-hypertension comparators, plus one Visium HD section from a fatal ICANS case.

For this component:

```bash
make fetch
make analysis
make verify
```

The repository-only spatial panel is a direct overlay on the deposited GSE269379 H&E image. Orange points denote 960 CXCL16-positive 8-µm bins; teal triangles denote 14 CXCR6-positive bins. It is not a synthetic image, and no smoothing, imputation, segmentation, or deconvolution was used. See `docs/PANEL_F_PROVENANCE.md` for exact source files, checksums, crop coordinates, and interpretation limits.

## Interpretation limits

Transcript detection does not establish protein abundance or a chemokine gradient. These data cannot distinguish soluble from membrane-bound CXCL16 and do not directly test chemotaxis, retention, egress, efficacy, toxicity, or causality. Product-stratified estimates are descriptive; the available cohorts do not support a formal time-by-product interaction.

## Source studies

- Sheih A, Voillet V, Hanafi LA, et al. *Nat Commun.* 2020;11:219. doi:10.1038/s41467-019-13880-1. GEO: GSE125881.
- Haradhvala NJ, Leick MB, Maurer K, et al. *Nat Med.* 2022;28:1848-1859. doi:10.1038/s41591-022-01959-0. GEO: GSE197268.
- Li Z, Zhao L, Zhang Y, et al. *Cell Rep.* 2023;42:113263. doi:10.1016/j.celrep.2023.113263. GEO: GSE162975.
- Maurer K, Grabski IN, Houot R, et al. *Blood.* 2024;144:2490-2502. doi:10.1182/blood.2024024381. GEO: GSE273170.
- Guerrero-Murillo M, Rill-Hinarejos A, Trincado JL, et al. *Cell Rep Med.* 2024;5:101803. doi:10.1016/j.xcrm.2024.101803. GEO: GSE235760.
- Lu IN, Müller-Miny L, Krekeler C, et al. *Genome Med.* 2025;17:71. doi:10.1186/s13073-025-01498-6. GEO: GSE269379.

Upstream files remain governed by their original repository terms. This repository does not redistribute complete source matrices, cell-level metadata, or source tissue images. Derived tables retain source-defined coded patient identifiers and the clinical categories needed to audit pairing. The spatial visualization includes a cropped, attributed histology background. See [`THIRD_PARTY_NOTICES.md`](../../THIRD_PARTY_NOTICES.md) for attribution and reuse conditions.
