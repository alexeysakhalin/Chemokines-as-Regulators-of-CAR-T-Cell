# Analysis specification: longitudinal public-cohort extension

Version: 0.1 (2026-09-02)

## Status and scope

This document fixes the analysis rules for the longitudinal extension before the repository workflow is implemented and frozen. It was prepared after public-data access and feasibility checks, some of which included provisional target-gene counts; it is therefore an analysis specification, not a prospective preregistration. Any later data-dependent change must be recorded in `DECISIONS.md` and labelled exploratory.

The extension tests whether CXCR6 transcription in CAR-T cells changes over time in independent public cohorts. It does not enlarge the GSE269379 ICANS-CSF donor comparison or its single-case spatial analysis, and it cannot test a CXCL16 protein gradient, chemotaxis, retention, egress, or causality.

## Cohort roles

- **GSE125881** is the four-patient discovery time course already present in Figure 6.
- **GSE197268** is the principal independent validation cohort for the infusion-product-to-day-7 contrast.
- **GSE162975** is an independent plate-based STRT-seq replication cohort and provides longer follow-up. It is analysed separately because its disease, CAR19/CAR22 cocktail, enrichment procedure, and platform differ from the 10x cohorts.
- **GSE273170** is an independent axi-cel sensitivity cohort for day 7 to day 14. Public matrices precede part of the authors' downstream QC and CAR-transcript-positive denominators are sparse, so this cohort is not promoted to a primary validation cohort.
- **GSE269379** remains an orthogonal ICANS-CSF and spatial dataset.

## Independent unit and sample construction

The patient is the independent biological unit. Cells, barcodes, technical libraries, plates, and repeated aliquots are never treated as independent replicates.

For each cohort, technical records belonging to the same patient and biological time point are collapsed before analysis. Each patient contributes at most once to a specified contrast. Matrices are not concatenated across studies, and no cross-study batch correction is used for inference.

### GSE197268

The first treatment course is used. The eligible sampling frame is the 21 patients with both an infusion-product matrix and a physically enriched D7-CART matrix. Matrices are joined to the authors' public cell metadata and CD4/CD8 classifications by the globalized 10x barcode. The principal cell stratum is author-QC-passing, product-matched CAR-transcript-positive, CD8-positive T cells at both time points. Patient 29 retreatment records, unsorted day-7 aliquots, and D7-CART aliquots from patients without an infusion product are excluded from the principal contrast and may be reported only as sensitivity analyses.

### GSE162975

The deposited gene-by-cell UMI matrix contains sequence-validated CD3-positive/CAR-positive cells. GEO sequencing records are collapsed to patient by canonical biological stage using the deposited `sampling stage` field; stage must not be inferred from the sample-title prefix because several titles disagree with the deposited stage. The principal stratum is all validated CAR-T cells because the public matrix does not provide a complete per-cell CD4/CD8 or CAR19/CAR22 assignment.

### GSE273170

Cells are called CAR-transcript-positive when the custom `CAR` feature has at least one UMI, matching the source study. The principal sensitivity contrast is paired day 7 versus day 14. Results are reported both without an additional denominator filter and with at least 10 CAR-transcript-positive cells at each time point. Baseline CAR signal is treated as background evidence, not as genuine pre-infusion CAR-T cells.

## Endpoints

For every patient and time point, the workflow records:

1. number of eligible cells;
2. number and fraction with at least one CXCR6 UMI;
3. total CXCR6 UMI;
4. total UMI in the eligible-cell pseudobulk;
5. CXCR6 counts per million total UMI;
6. detection summaries for the locked memory-, effector-, and dysfunction-associated marker sets already used in the original module.

The principal interpretable endpoint is the within-patient change in the CXCR6-detected cell fraction. Pseudobulk CXCR6 counts per million are a sensitivity endpoint. No missing time point is imputed, and a sample with no recoverable CAR-T cells is missing for within-CAR-T expression rather than assigned zero expression.

## Locked contrasts

- GSE197268: infusion product to D7-CART (principal early validation).
- GSE162975: T0 infusion product to T1 peak within the first month (independent early replication).
- GSE162975: T1 to the first available T3 or T4.5 sample (exploratory post-peak change).
- GSE162975: T1 to the first available T6-or-later sample (exploratory extended follow-up).
- GSE273170: day 7 to day 14 (sensitivity replication).

The advertised study cohort size is never substituted for the number of complete eligible patient pairs.

## Statistical summaries

Each cohort and contrast is analysed separately. The workflow reports the paired patient values, median paired change, interquartile range, number increasing/decreasing/unchanged, and a two-sided exact binomial sign-test P value after removing exact zero differences. A paired Wilcoxon result is reported only when its assumptions and exact handling of ties/zeros are explicit. Confidence intervals are patient-level bootstrap intervals with a fixed seed and are labelled descriptive when the number of pairs is small.

The GSE197268 and GSE162975 early contrasts form one two-test family; Holm-adjusted values are reported. Later contrasts and GSE273170 are supportive or exploratory and cannot rescue a failed early validation. Marker-set comparisons are adjusted by Benjamini-Hochberg within the prespecified marker family.

For GSE197268, product-stratified patient estimates for axi-cel and tisa-cel are always shown because product-specific CAR detection and biology differ. Response-stratified estimates are exploratory; the cohort is too small for a high-dimensional adjusted outcome model.

## Sensitivity analyses

- minimum eligible-cell thresholds of 25 and 100 for GSE197268;
- GSE197268 all CAR-positive T cells versus the locked CD8-positive stratum;
- GSE197268 physical D7-CART fraction versus product-matched transcript-positive cells;
- exclusion of any very small GSE162975 biological specimen, reported by an explicit threshold rather than by CXCR6 value;
- GSE273170 minimum CAR-positive denominators of 1 and 10;
- patient-level detection fraction versus pseudobulk counts per million;
- leave-one-patient-out median effects for cohorts with at least five eligible pairs.

## Cross-cohort interpretation

The primary synthesis is a cohort-level table, not a pooled cell analysis or a pooled P value. A common temporal programme is described only when cohort-specific directions are concordant. Discordance is reported as heterogeneity rather than hidden through pooling.

Even concordant transcript-level changes support only a temporal association. They do not establish enhanced tumour or CNS trafficking, a CXCL16 gradient, membrane-versus-soluble CXCL16 biology, treatment efficacy, neurotoxicity, or a causal role for CXCR6.

## Required reproducibility outputs

- checksum-pinned source manifests and sample crosswalks;
- per-patient/time-point QC and aggregate tables;
- one explicit inclusion table per contrast;
- cohort-specific paired results and sensitivity results;
- a deterministic Markdown report suitable for manuscript drafting;
- a machine-readable run manifest and SHA-256 manifest of frozen outputs;
- tests preventing duplicate patients, duplicate biological samples, and cell-level pseudoreplication.

