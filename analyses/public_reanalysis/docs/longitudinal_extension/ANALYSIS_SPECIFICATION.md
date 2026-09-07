# Analysis specification: longitudinal public-cohort extension

Version 0.2, finalized 2026-09-07

## Status and scope

This is an exploratory secondary analysis of public data, not a prospective preregistration or a confirmatory clinical study. The rules below were fixed after data-access and feasibility checks. Any data-dependent departure must be entered in `DECISIONS.md`, propagated to the run manifest, and labelled exploratory.

The primary question is whether the fraction of CAR-T cells with detectable CXCR6 transcript changes within patients between the manufactured product and an early post-infusion sample. The patient is the independent biological unit. The extension does not test a CXCL16 protein gradient, chemotaxis, retention, egress, treatment efficacy, neurotoxicity, or causality.

This is a documented, targeted public-data reanalysis rather than an exhaustive meta-analysis. A cohort enters the finalized workflow only after its public cell identity, CAR definition, biological time points, and patient pairing have been audited and harmonized. Additional public cohorts may merit separate future audit; their existence is not treated as negative evidence.

## Cohort roles

- **GSE125881** is the four-patient discovery time course in the original repository module. It remains descriptive and is not included in the three-test early family.
- **GSE197268** is an independent early comparison of infusion-product and physically enriched day-7 CAR-T cells in large B-cell lymphoma.
- **GSE235760** is an independent paired early comparison of manufactured varni-cel and the patient-specific peripheral-blood expansion peak in five adults with B-ALL.
- **GSE162975** is an independent plate-based STRT-seq replication cohort and provides later follow-up. It is kept separate because disease, CAR19/CAR22 cocktail, enrichment, and platform differ from the 10x cohorts.
- **GSE273170** is a supportive axi-cel day-7-to-day-14 sensitivity cohort. Deposited matrices precede part of the authors' downstream quality control, and CAR-transcript-positive denominators are sparse.
- **GSE269379** remains an orthogonal ICANS-CSF and single-section spatial dataset and is never pooled with blood time courses.

## Independent unit and sample construction

Cells, barcodes, plates, technical libraries, and repeated aliquots are not biological replicates. Technical records from the same patient and biological time point are collapsed before contrasts. Each patient contributes at most one paired value to a defined contrast. Studies are not concatenated, batch-corrected together, or combined into a pooled P value.

### GSE197268

The first treatment course is used. The sampling frame contains 21 patients with author-QC metadata for both an infusion product and a physically enriched D7-CART sample. Before inspecting CXCR6 counts, the two author annotation files are joined by globalized 10x barcode and the denominator is restricted to product-matched cells with `CAR == True` and author subtype `CD8 T`. Nine patients have at least 25 such cells at both stages; the other 12 fail this outcome-blind denominator rule. Only the corresponding 18 GEO processed matrices are downloaded for the primary raw-count reconstruction. The committed screening table reports all 21 patients and the exclusion reason. Patient 29 retreatment records, unsorted day-7 aliquots, and D7-CART aliquots without a matching product are excluded from the principal contrast.

### GSE235760

The checksum-pinned CELLxGENE H5AD contains 37,100 author-QC T cells and integer UMI counts in `raw.X`. `obs/donor_id` defines patient, `obs/Timepoint` defines `IP` and patient-specific `Peak`, `obs/Sample_id` identifies the biological sample, and `obs/Transduction == "CAR+"` defines the eligible stratum. CXCR6 is the unique `raw.var/feature_name == "CXCR6"` feature (ENSG00000172215). The final eligible set contains 18,978 CAR-positive cells: 10,516 at IP and 8,462 at Peak.

### GSE162975

The deposited gene-by-cell UMI matrix contains sequence-validated CD3-positive/CAR-positive cells. GEO sequencing records are collapsed by patient and deposited `sampling stage`; a stage is not inferred from title prefix because several titles disagree with the deposited stage. All validated CAR-T cells are used because complete per-cell CD4/CD8 and CAR19/CAR22 assignments are unavailable.

### GSE273170

The checksum-pinned GEO archive contains gzip-compressed, dense R-style gene-by-cell RNA tables: quoted identifiers are separated by whitespace, the header contains the cell barcodes, and each data row begins with a quoted gene symbol. The parser accepts this deposited layout directly and verifies the SHA-256 digest of every selected member against `config/gse273170_crosswalk.tsv`. Cells are CAR-transcript positive when the custom `CAR` feature has at least one UMI, matching the source study. Day 7 and day 14 are paired within patient. Results are reported for all recoverable pairs and for the locked minimum of 10 CAR-transcript-positive cells at both stages. Baseline CAR signal is not interpreted as genuine pre-infusion CAR-T biology.

### GSE290722 exclusion

GSE290722 was excluded from the declared longitudinal CAR-T endpoint because the public repeated-timepoint metadata do not identify CAR-T cells: the only public author-called CAR-positive barcode list covers 128 week-4 cells from 13 patients, leaving no infusion-product-to-post-infusion patient pairs under the locked CAR-T definition. The cohort may be useful for a separate native-immune or week-4 descriptive analysis, but it cannot increase the eligible n for the current endpoint.

## Endpoints

For each patient and stage, the workflow records the number of eligible cells, the number and fraction with at least one CXCR6 UMI, total CXCR6 UMI, total raw UMI where recoverable, and CXCR6 counts per million total UMI. The inferential endpoint is the within-patient change in CXCR6-detected cell fraction; CPM is retained as a descriptive normalization only.

Memory-, effector-, and dysfunction-associated marker sets are locked as follows:

- memory: CCR7, IL7R, LEF1, LTB, MAL, SELL, TCF7;
- effector: CCL5, GNLY, GZMB, IFNG, NKG7, PRF1;
- dysfunction: ENTPD1, HAVCR2, LAG3, PDCD1, TIGIT, TOX.

Marker-set detection is the number of nonzero gene-cell pairs divided by the number of possible gene-cell pairs. It is a descriptive summary, not a cell-state label or a module score. Any paired marker-set sign test is restricted to the same included patients, eligible cells, minimum-cell rule, and preselected stages as the corresponding CXCR6 contrast. No missing stage is imputed, and a stage with zero recoverable eligible CAR-T cells is missing rather than assigned zero CXCR6 expression. For GSE162975, the deposited-to-analysis stage mapping and the priority used to select the first available T3/T4.5 or T6-or-later sample are defined in `config/gse162975_crosswalk.tsv`; the same selected source stage is aliased in both the CXCR6 and marker aggregate tables.

## Locked contrasts

- GSE197268: infusion product to D7-CART, minimum 25 eligible cells at both stages;
- GSE235760: infusion product to patient-specific expansion peak, minimum 25 eligible cells at both stages;
- GSE162975: T0 infusion product to T1 first-month peak, minimum 25 eligible cells at both stages;
- GSE162975: T1 to first available T3/T4.5, exploratory post-peak comparison;
- GSE162975: T1 to first available T6-or-later stage, exploratory extended follow-up;
- GSE273170: day 7 to day 14, minimum 10 CAR-transcript-positive cells at both stages.

The advertised study cohort size is never substituted for the number of complete eligible patient pairs.

## Statistical summaries

Each cohort and contrast is analyzed separately. The workflow reports paired patient values, the median and interquartile range of paired changes, counts increasing/decreasing/unchanged, and a two-sided exact binomial sign-test P value after discarding exact zero changes. There are no cell-level P values.

The GSE197268, GSE235760, and GSE162975 early comparisons form the three-test early family in the finalized exploratory workflow; Holm-adjusted values are reported. Later GSE162975 contrasts and GSE273170 are supportive or exploratory. Marker-set tests are Benjamini-Hochberg adjusted within cohort and contrast.

GSE197268 product-stratified axi-cel and tisa-cel estimates are descriptive. The available patients do not support a formal time-by-product interaction or a high-dimensional adjusted outcome model.

## Sensitivity analyses

- product-stratified descriptive estimates for axi-cel and tisa-cel;
- thresholds of 25, 100, 500, and 1,000 cells for GSE235760;
- removal of explicitly small GSE162975 specimens by denominator rather than by CXCR6 value;
- GSE273170 minimum CAR-positive denominators of 1 and 10;
- leave-one-patient-out median effects for contrasts with at least five pairs.

## Cross-cohort interpretation

The synthesis is a cohort-level table, not a pooled cell analysis. A common temporal program is described only if independent cohort directions are concordant. Discordance is reported as heterogeneity. Even concordant transcript changes would support temporal association only and would not establish tumor or CNS trafficking, a functional gradient, soluble-versus-membrane CXCL16 biology, or causality.

## Required reproducibility outputs

- checksum-pinned source manifests and sample crosswalks;
- patient/stage aggregates and contrast-inclusion tables;
- cohort-specific paired summaries and sensitivity results;
- a deterministic Markdown report;
- a machine-readable run manifest and internal SHA-256 manifest;
- tests preventing duplicate patients, duplicate cells, and cell-level pseudoreplication.

The run manifest uses repository-relative paths for version-controlled code and configuration, and content-addressed keys for downloaded matrices. Verification compares the observed provenance with the frozen provenance, confirms that frozen code/source/audit digests still match the working tree, and checks every declared canonical output byte for byte.
