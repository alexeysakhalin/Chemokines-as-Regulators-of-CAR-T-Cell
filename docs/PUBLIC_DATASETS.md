# Public datasets relevant to the proposed analyses

## Scope and evidence boundary

This inventory identifies public data that may support parts of the prospective analysis plan. It was initially assembled from repository records and linked publications checked on 2026-09-01. The completed public-data reanalysis, including cohort-specific eligibility and results, is documented separately in [`analyses/public_reanalysis/README.md`](../analyses/public_reanalysis/README.md). This inventory retains the wider feasibility assessment; accessions listed here are not all included in the scientific analysis.

The five requested analyses are abbreviated as follows:

1. **A1 — single-cell state analysis:** scRNA-seq/CITE-seq of the CAR-T infusion product and serial blood, with patient-level memory-like, effector, proliferating, and dysfunction-associated states;
2. **A2 — spatial analysis:** spatial localization of chemokine sources relative to vessels, stroma, tumor nests, lymphatics, and CAR-T cells;
3. **A3 — protein validation:** orthogonal measurement of chemokine protein gradients, explicitly separating soluble and membrane-associated CXCL16;
4. **A4 — functional trafficking:** chemotaxis, retention, and egress assays;
5. **A5 — matched product comparison:** a chemokine-modified CAR-T product compared with an otherwise identical CAR-T product lacking only the chemokine module.

`Strong` means that the accession can address the named analysis directly, subject to the stated limitations. `Partial` means that it supports only a component or an exploratory case study. `No` means that the required measurement or causal comparison is absent.

## At-a-glance suitability

| Accession | Independent human or model units represented in the accession | Principal modality and compartment coverage | A1 | A2 | A3 | A4 | A5 |
|---|---:|---|---|---|---|---|---|
| [GSE125881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125881) | 4 patients | 5′ scRNA-seq + scTCR; product and three serial blood time points | Partial | No | No | No | No |
| [GSE197268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197268) | 32 patients | 5′ scRNA-seq; baseline blood, product, predominantly day-7 blood, sparse day-14 blood | Strong | No | No | No | No |
| [GSE273170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273170) | 10-patient single-cell subset | 5′ scRNA-seq, scTCR, and 32-marker ADT; partial baseline/product and day-7/day-14 blood | Strong | No | No | No | No |
| [GSE166352](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166352) | 3 patients | 3′ scRNA-seq; product and two serial blood time points | Partial | No | No | No | No |
| [GSE235760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235760) | 5 patients | 5′ scRNA-seq + scTCR; product and one patient-specific blood expansion peak, CAR-positive and CAR-negative fractions | Partial | No | No | No | No |
| [GSE197215](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197215) | 12 patients | scRNA-seq + 14-marker CITE-seq; infusion product under four ex vivo stimulation conditions | Partial | No | No | No | No |
| [GSE253352](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253352) | 18 patients represented by product scRNA-seq; 24 treated in the linked clinical study | scRNA-seq, CAR-PCR, and scTCR; product, selected PBMC, and selected tumor-biopsy samples | Partial | No | No | No | No |
| [GSE310353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310353), including [GSE282302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282302) and [GSE310352](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310352) | 39 patients in Visium; 9 tumors in the CosMx subset | Whole-transcriptome Visium with histology plus targeted single-cell CosMx and limited imaging proteins; untreated-status claim is on QC hold | No | Partial | No | No | No |
| [GSE332850](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332850) | 13 patients, 18 FFPE sections | Xenium targeted spatial transcriptomics of HCC biopsies | No | Strong | No | No | No |
| [GSE325706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325706) + [GSE325772](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325772) | 1 patient in these two accessions | Visium HD tumor sections plus 5′ scRNA-seq/scTCR of product and serial blood | Partial | Partial | No | No | No |
| [GSE335372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335372) | 2 spatial samples; independent mouse count not stated in the Series record | DBiT spatial transcriptomics of HT29 lung xenografts after CAR-T or CAR-T–drug-conjugate treatment | No | Partial | No | No | No |

No listed accession can execute A3, A4, or A5. Transcript abundance cannot distinguish soluble from membrane CXCL16, spatial RNA is not a protein gradient, changes in circulating cell abundance are not chemotaxis measurements, and comparisons between different clinical products or cohorts are not split-product comparisons that isolate a chemokine module.

## Dataset records

### GSE125881 — long serial follow-up of sorted CD8 CAR-T cells

- **Study and modality:** 10x Chromium 5′ gene expression with paired V(D)J/scTCR sequencing of sorted CD8-positive, EGFRt-positive CAR-T cells. The linked paper reports 62,167 quality-controlled cells.
- **Biological coverage:** four patients selected for durable CAR-T persistence: two with CLL and two with NHL. Each contributes an infusion product and three post-infusion blood samples. The 16 biological specimens appear as 32 GEO records because gene-expression and V(D)J libraries are separate. Exact series time points are CLL1: product/d21/d38/d112; CLL2: product/d12/d29/d83; NHL6: product/d12/d29/d102; and NHL7: product/d12/d28/d89.
- **Access:** processed matrices and clonotype outputs are public in GEO; raw sequence reads are public through SRA [SRP182902](https://www.ncbi.nlm.nih.gov/sra/?term=SRP182902) and BioProject [PRJNA517816](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA517816). No controlled tier is identified in the GEO record.
- **Local integrity check:** the retrieved processed expression matrix is 114,573,707 bytes with SHA-256 `243e6e2dd09158e08dcaecfbf35b23122ac3096c1039355b7590474e76d8c02c`; `GSE125881_RAW.tar` is 2,048,000 bytes with SHA-256 `769fa5b71b12b7efe433f51164d4370e0af00c20bebc647029cc4fbd9a943487`; and `GSE125881_family.soft.gz` is 5,747 bytes with SHA-256 `13a3085bc28fb27228ca858ed935d74b936dbd3f5efc376e485f73ef7b62eb1e`. The processed-matrix header contains 62,167 cell columns. Barcode suffixes run from `-1` through `-16`; their specimen mapping is defined by the GEO SOFT field `id in raw.expmatrix.csv` and must be used instead of assuming suffix order from patient names.
- **Suitability:** useful for within-patient CD8 CAR-T state and clonotype trajectories from product to late persistence. It is unusually valuable for long follow-up.
- **Exact limitations:** only four selected durable persisters are represented; only sorted CD8 CAR-T cells were profiled; there is no CITE-seq, host T-cell compartment, pretreatment PBMC, tumor, marrow, or spatial data. Serial blood abundance cannot distinguish proliferation, survival, dilution, and redistribution. The primary paper includes larger TCRB and integration-site cohorts that must not be attributed to this four-patient GEO single-cell cohort.
- **Paper:** Sheih et al., *Nature Communications* (2020), [DOI 10.1038/s41467-019-13880-1](https://doi.org/10.1038/s41467-019-13880-1), [PMID 31924795](https://pubmed.ncbi.nlm.nih.gov/31924795/). Published analysis code is linked at [ValentinVoillet/CAR-T](https://github.com/ValentinVoillet/CAR-T).

### GSE197268 — largest open processed product/early-blood cohort in this set

- **Study and modality:** 10x Chromium 5′ scRNA-seq with coupled scTCR in the study; no CITE-seq or ADT. The paper reports 602,577 high-quality cells from adults with refractory LBCL treated with axi-cel or tisa-cel.
- **Biological coverage:** 32 patients: 19 axi-cel and 13 tisa-cel. The paper's main analytic set contains 105 samples; the current GEO record has 109 sample records after additional retreatment records. Coverage includes 20 baseline PBMC samples, 31 initial infusion products, day-7 samples that are either sorted CAR-positive, corresponding CAR-negative, or unsorted, and four unsorted day-14 samples. Twenty-eight of 32 patients have an initial product and day-7 record, but the exact paired set depends on whether the analysis requires explicitly sorted CAR-positive cells.
- **Access:** per-sample processed gene-expression matrices are open in GEO. The GEO aggregate file contains processed 10x matrices; raw reads are not present in the GEO sample records. Raw scRNA/scTCR data are controlled through dbGaP [phs002922.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002922.v1.p1), with a broader version available as phs002922.v2.p1.
- **Suitability:** the strongest first discovery dataset for patient-level product-to-early-blood state analysis, response-stratified summaries, and product-specific sensitivity analyses. Public processed matrices are sufficient for replicate-aware state fractions and pseudobulk analyses if the required cell metadata and CAR assignment are available or can be reconstructed reproducibly.
- **Exact limitations:** substantial sample-pattern missingness; follow-up is predominantly day 7; only four day-14 samples; no tumor or spatial data; RNA-only phenotype in the public tier; raw TCR reprocessing requires controlled access. Axi-cel versus tisa-cel is a confounded product comparison, not an identical CAR with or without a module. Blood-frequency changes cannot be interpreted as trafficking.
- **Paper:** Haradhvala et al., *Nature Medicine* (2022), [DOI 10.1038/s41591-022-01959-0](https://doi.org/10.1038/s41591-022-01959-0), [PMID 36097221](https://pubmed.ncbi.nlm.nih.gov/36097221/). Published analysis code is linked at [getzlab/Haradhvala_et_al_2022](https://github.com/getzlab/Haradhvala_et_al_2022).

### GSE273170 — CITE-seq-supported early longitudinal LBCL subset

- **Study and modality:** 10x 5′ scRNA-seq, scTCR-seq, and CITE-seq/ADT. The GEO accession contains 81 assay-level records and the linked paper reports 87,432 quality-controlled cells for the new single-cell subset.
- **Biological coverage:** ten axi-cel-treated LBCL patients, four durable responders and six nonresponders. There are 30 biological specimens: five infusion products, five baseline/apheresis samples, ten day-7 PBMC samples, and ten day-14 PBMC samples. RNA, TCR, and ADT are not complete for every specimen. The paper's prospective clinical cohort comprises 50 patients, but the accession's new single-cell subset is ten patients; a combined 28-patient analysis in the paper also incorporates a prior cohort and must not be labeled as GSE273170 alone.
- **Access:** processed RNA, TCR, and ADT tables plus the antibody list are public in GEO. Raw and individual-level sequence data are controlled through dbGaP [phs002922.v2.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002922.v2.p1). Raw reads are not exposed by the GEO record.
- **Suitability:** strong for early state and clonotype analyses and for protein-supported high-level phenotype annotation. It can provide an independent early-time-point check of findings developed in GSE197268, using only those patients with the required paired specimens.
- **Exact limitations:** only five products and five baseline specimens; no tumor or spatial data; incomplete modality coverage. The 32-marker ADT panel does not include CXCL16, CXCR6, CXCR3, or CXCR4, so it cannot validate the key chemokine proteins or soluble/membrane CXCL16. It contains no functional migration measurements or module-matched control.
- **Paper:** Maurer et al., *Blood* (2024), [DOI 10.1182/blood.2024024381](https://doi.org/10.1182/blood.2024024381), [PMID 39241199](https://pubmed.ncbi.nlm.nih.gov/39241199/).

### GSE166352 — small serial cohort of non-viral PD1-targeted CAR-T cells

- **Study and modality:** 10x Chromium 3′ v3.1 scRNA-seq without CITE-seq/ADT or V(D)J sequencing.
- **Biological coverage:** three B-NHL patients treated with non-viral PD1-locus-integrated anti-CD19 PD1-19bbz CAR-T cells. There are nine biological samples: patient 1 product/d12/d29, patient 2 product/d12/d28, and patient 3 product/d7/d28. The linked clinical paper includes eight treated patients, but only three are represented in this accession.
- **Access:** processed MTX/TSV files are open in GEO; raw reads are open in SRA [SRP305312](https://www.ncbi.nlm.nih.gov/sra/?term=SRP305312) and BioProject [PRJNA700587](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA700587).
- **Suitability:** exploratory paired product-to-two-time-point blood analysis within a distinctive engineered product.
- **Exact limitations:** three patients, no protein or clonotype modality, no tumor or spatial data, and one late sample with too few CAR-positive cells for some reported state plots. All clinical samples received the PD1-19bbz product; there is no otherwise identical conventional-CAR control in the accession.
- **Paper:** Zhang et al., *Nature* (2022), [DOI 10.1038/s41586-022-05140-y](https://doi.org/10.1038/s41586-022-05140-y), [PMID 36045296](https://pubmed.ncbi.nlm.nih.gov/36045296/).

### GSE235760 — paired CAR-positive/CAR-negative fractions at product and expansion peak

- **Study and modality:** 10x 5′ scRNA-seq and paired scTCR-seq. The paper reports 37,100 cells.
- **Biological coverage:** five adults with B-ALL treated with varni-cel. For each patient, CAR-positive and CAR-negative T-cell fractions were sorted from the infusion product and from one peripheral-blood expansion-peak sample. This yields 20 biological fractions and 40 GEO records because RNA and TCR libraries are separate. The expansion peak is the patient-specific maximum within the first four weeks, not a common nominal day and not a serial time course.
- **Access:** processed files are open in GEO and raw reads are open through SRA/BioProject [PRJNA987581](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA987581). No controlled tier is stated for this accession.
- **Suitability:** useful for within-patient CAR-positive versus non-transduced T-cell state and clonotype contrasts from product to peak.
- **Exact limitations:** five patients and only one post-infusion time point; the day varies by patient; no CITE-seq, tumor, or spatial data. CAR-negative cells are non-transduced cells, not an otherwise identical CAR-T product lacking a chemokine module, and therefore cannot satisfy A5. Larger validation cohorts in the paper are not part of this GEO accession.
- **Paper:** Guerrero-Murillo et al., *Cell Reports Medicine* (2024), [DOI 10.1016/j.xcrm.2024.101803](https://doi.org/10.1016/j.xcrm.2024.101803), [PMID 39471818](https://pubmed.ncbi.nlm.nih.gov/39471818/). A published correction should be retained in the citation and reproduction record: [DOI 10.1016/j.xcrm.2025.102026](https://doi.org/10.1016/j.xcrm.2025.102026).

### GSE197215 — antigen-specific functional landscape of infusion products

- **Study and modality:** single-cell RNA profiling plus 14-marker CITE-seq of CAR-T infusion products under four ex vivo conditions: CD19-specific CAR stimulation, TCR stimulation, non-target APC stimulation, and unstimulated control. The paper reports 101,326 transcriptomes and 97,981 cells after filtering.
- **Biological coverage:** infusion products from 12 pediatric or young-adult ALL patients: five durable complete responders, five patients with CD19-positive relapse, and two nonresponders. GEO contains 120 assay-level records reflecting patients, conditions, and modalities.
- **Access:** an integrated processed RDS object is open in GEO and the GEO records link raw data through SRA/BioProject [PRJNA809371](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA809371). No controlled component is identified in the Series record.
- **Suitability:** strong for product-only memory/effector and antigen-response heterogeneity, including patient-level comparisons across controlled stimulation conditions. It is a useful product-state reference, not a longitudinal dataset.
- **Exact limitations:** no post-infusion blood, tumor, or spatial samples; no serial persistence analysis; no scTCR in the stated accession design. The 14 ADTs cover T-cell identity, differentiation, activation, and inhibitory markers but not soluble CXCL16 or a tissue gradient. Ex vivo stimulation responses must not be presented as in vivo trafficking.
- **Paper:** Bai et al., *Science Advances* (2022), [DOI 10.1126/sciadv.abj2820](https://doi.org/10.1126/sciadv.abj2820), [PMID 35675405](https://pubmed.ncbi.nlm.nih.gov/35675405/).

### GSE253352 — solid-tumor CAR-T product, blood, and biopsy transcriptomes

- **Study and modality:** scRNA-seq with CAR-PCR and scTCR libraries for selected compartments. The linked first-in-human study evaluated conventional GPC3 CAR-T and IL-15-coexpressing GPC3 CAR-T cells.
- **Biological coverage:** the clinical study treated 24 patients in two cohorts. The GEO Series contains 83 assay-level records. Product mRNA records represent 12 IL-15-CAR patients and six conventional-CAR patients. PBMC mRNA is available for selected patients, and tumor-biopsy mRNA is present for eight IL-15-CAR patients; TCR and CAR-PCR coverage is also incomplete by patient and compartment. The record does not provide a balanced product–blood–tumor series for all 24 treated patients, and sample titles do not establish a common serial time-point schedule.
- **Access:** raw sequence data are linked through SRA/BioProject [PRJNA1056257](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1056257); GEO provides open raw-count, scaled-count, metadata, and UMAP-coordinate tables for several integrated objects. No controlled component is identified in the GEO record, but clinical variables outside the deposited metadata remain governed by the study and publication terms.
- **Suitability:** valuable as a solid-tumor, multi-compartment single-cell dataset and for exploratory study of product versus circulating versus tumor-infiltrating states where the same patient can be linked. It is the best listed transcriptomic accession for testing whether a state signature developed in hematologic malignancies transfers to tumor-infiltrating CAR-T cells.
- **Exact limitations:** coverage is highly unbalanced across compartments and modalities; collection-day metadata must be verified before calling any comparison longitudinal. It is not spatial transcriptomics. Conventional and IL-15 products belong to different clinical cohorts and patients, and IL-15 is a cytokine armoring module rather than the chemokine module requested in A5. The comparison is neither paired nor isolated from trial-cohort, dose, and patient differences. There is no soluble-versus-membrane CXCL16 assay or functional migration assay.
- **Paper:** Steffin et al., *Nature* (2025), [DOI 10.1038/s41586-024-08261-8](https://doi.org/10.1038/s41586-024-08261-8), [PMID 39604730](https://pubmed.ncbi.nlm.nih.gov/39604730/).

### GSE310353 — multi-modal PDAC spatial reference with a metadata QC hold

- **Study and modality:** this is an observational spatial atlas of primary resected PDAC, not a CAR-T study. The paper reports 341,949 Visium spots across 108 H&E sections from 39 tumors and 531,718 CosMx cells across 147 tumor-enriched fields of view from 9 tumors. Visium provides transcriptome-wide 55-µm spots with an estimated mean of 18 cells per spot. CosMx used a fixed 969-gene human panel with 10 added PDAC markers, DAPI, and CD298, PANCK, CD45, and CD3 imaging proteins. Sequential sections supplied H&E and, where available, trichrome-derived collagen measurements.
- **Access:** [GSE310353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310353) is the superseries. [GSE282302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282302) contains 108 Visium records and a 4.9-GiB series archive. Each record provides filtered and raw feature-barcode H5 matrices, a high-resolution tissue PNG, tissue-position CSV, and scale-factor JSON. [GSE310352](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310352) contains eight CosMx slide records and a 5.9-GiB series archive. Each record provides cell-by-gene counts, cell metadata, a raw-core archive with molecule detections, segmentation polygons and field positions, and a Seurat RDS object; the GEO record states that morphology images were omitted, while the coordinates needed to reconstruct cell and molecule positions were retained. [GSE310388](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310388) is a separate hypoxia-versus-normoxia cell-line RNA-seq series and is not required for the primary spatial source map.
- **Panel audit:** the header of the exact public file `GSM9294406_slide5b3_cell_by_gene_counts.csv.gz` was inspected. It contains `CXCL16`, `CXCR6`, `CXCL9`, `CXCL10`, `CXCL12`, `CXCL14`, `CXCR4`, `CCL2`, `CCL19`, `CCL21`, and `CCR7`; `CXCL11` is absent. The panel is fixed across the experiment, but the workflow must still assert feature presence in every imported file before analysis.
- **Published code:** the paper links analysis repositories under [anvaly](https://github.com/anvaly), including [SpaRTa_workfow](https://github.com/anvaly/SpaRTa_workfow), [SpaRTa_AirFlow](https://github.com/anvaly/SpaRTa_AirFlow), and [SpatialPortalV2](https://github.com/anvaly/SpatialPortalV2). Archived releases are identified by Zenodo DOIs [10.5281/zenodo.17817599](https://doi.org/10.5281/zenodo.17817599), [10.5281/zenodo.17817567](https://doi.org/10.5281/zenodo.17817567), [10.5281/zenodo.17817587](https://doi.org/10.5281/zenodo.17817587), and [10.5281/zenodo.17818928](https://doi.org/10.5281/zenodo.17818928).
- **Suitability:** this is the preferred external PDAC reference for the non-CAR-T portion of A2: mapping candidate ligand-producing tumor, fibroblast, endothelial, myeloid, and lymphoid compartments; measuring their relation to tumor nests, stroma, vessels, and collagen; and validating Visium source assignments at single-cell resolution with CosMx. CXCL16 expression can be analyzed as an exploratory transcript-source endpoint. The multi-resolution Visium-to-CosMx design, concentric tumor-distance bins, Delaunay neighborhoods, and registered histology are methods worth adapting.
- **Inference adaptation:** the published spot- and cell-level workflows are useful for annotation and feature discovery, but the publication analysis must summarize or model results at the patient level. Sections, regions, fields, spots, and cells remain nested observations. Use patient and section random effects where justified, patient-blocked tissue-constrained permutations, patient-level or hierarchical bootstrap intervals, and multiplicity control across the prespecified ligand-receptor family. A single-cell Wilcoxon test or a spot-level regression alone must not be presented as the cohort-level inferential result.
- **Metadata QC hold:** the article repeatedly describes the human cohort as untreated or treatment-naive and states that patient Visium used the human transcriptome v1 probe set. In contrast, the current GSE282302 Series and all 108 sample records label the tissue as post-neoadjuvant, state treatment with at least three cycles of FOLFIRINOX, and identify the Visium human probe set v2. GSE310352 describes treatment-naive resections. This unresolved conflict blocks any treatment-status-stratified or confirmatory biological claim. Before scientific use, obtain a correction or written clarification from the authors/GEO and record the resolved sample-level treatment and reference version in the manifest. Until then, the data may be used only as a treatment-agnostic spatial method/reference analysis with the discrepancy disclosed.
- **Exact limitations:** no cells in this atlas are CAR-T cells, no infusion product or serial blood is present, and the four CosMx imaging proteins do not measure CXCL16. RNA cannot distinguish membrane-bound from shed CXCL16, quantify a soluble protein gradient, or establish chemotaxis, retention, or egress. The dataset therefore cannot address A1, A3, A4, or A5 and cannot measure CAR-T localization within A2.
- **Paper:** Lyubetskaya et al., *Cell Reports* (2026), [DOI 10.1016/j.celrep.2025.116827](https://doi.org/10.1016/j.celrep.2025.116827), [PMID 41533516](https://pubmed.ncbi.nlm.nih.gov/41533516/).

### GSE332850 — patient-level Xenium spatial HCC biopsies

- **Study and modality:** Xenium in situ targeted spatial transcriptomics using the Xenium Prime 5K Human Pan Tissue and Pathways panel, with single-cell coordinates, transcript locations, cell and nucleus boundaries, morphology images, and processed feature matrices.
- **Biological coverage:** 18 FFPE needle-biopsy sections from 13 patients with HCC treated with dnTGFβRII-armored GPC3 CAR-T cells. GEO labels sections as pretreatment, partial response, or non-partial-response. Only a minority of patients have more than one response-status section, so most contrasts are cross-sectional rather than paired.
- **Access:** the Xenium data and images are open in GEO. The linked paper states that deidentified individual patient data are under controlled access in Science Data Bank [10.57760/sciencedb.39391](https://doi.org/10.57760/sciencedb.39391), and that a separate Digital Spatial Profiler dataset is deposited in GSA-Human as [HRA010444](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA010444). Source data accompany the paper.
- **Suitability:** strongest listed dataset for A2. It can support cell-type-resolved ligand-source maps, patient-stratified distances, and tissue-constrained neighborhood tests if the required ligand, receptor, endothelial, stromal, tumor, and lymphatic markers are present in the targeted panel and sufficiently detected.
- **Exact limitations:** targeted rather than whole-transcriptome panel; candidate-gene presence must be checked before preregistration. Eighteen sections do not imply 18 independent patients. Response-status groups are not a randomized product comparison, and repeated sections must be nested within patient. The mapping from 18 GEO samples to numbered exported file bundles must be validated. RNA signal cannot establish secreted protein, proteolytic processing, glycosaminoglycan presentation, or a directional protein gradient. No product, serial blood, functional migration, or parent-CAR control is present in this accession.
- **Paper:** Zhang et al., *Nature* (2026), [DOI 10.1038/s41586-026-10786-z](https://doi.org/10.1038/s41586-026-10786-z), [PMID 42457964](https://pubmed.ncbi.nlm.nih.gov/42457964/).

### GSE325706 and GSE325772 — linked single-patient GCAR1 spatial and serial data

- **Study and modality:** GSE325706 contains two Visium HD spatial gene-expression samples from the ASPS-02 participant. GSE325772 contains eight 5′ scRNA-seq libraries and eight matched V(D)J/scTCR libraries from the same participant: apheresis, infusion product, and blood on days 8, 15, 17, 22, 31, and 46.
- **Biological coverage:** one participant with metastatic alveolar soft-part sarcoma treated with GPNMB-directed GCAR1. The publication identifies GSE325706 as two ASPS-02 spatial samples. A separate xenograft spatial dataset is GSE282057 and must not be counted as part of GSE325706. GSE325772 is explicitly a one-participant serial blood/product series even though the linked clinical report discusses a broader study context.
- **Access:** raw and processed files are open in GEO/SRA for both accessions. The paper additionally identifies controlled deep-sequencing data at EGA [EGAS50000001696](https://ega-archive.org/studies/EGAS50000001696) and states that deidentified participant-level clinical data require reviewed request and a data-use agreement.
- **Suitability:** together they provide an unusually coherent illustration of product, serial circulation, and tumor spatial context for one participant. They can be used for workflow feasibility, qualitative concordance, and a clearly labeled case study.
- **Exact limitations:** one patient cannot support population inference or patient-level confidence intervals. The scRNA record lists only apheresis, product, and serial blood sample names; any claim that tumor cells are contained in GSE325772 must be verified rather than inferred from the Series narrative. Visium HD spots are not automatically single cells and require validated segmentation/deconvolution. The spatial pre/post interpretation requires exact section timing and treatment metadata. GCAR1 is not a chemokine-modified product, and no identical product control is present.
- **Paper:** Zemp et al., *Nature Cancer* (2026), [DOI 10.1038/s43018-026-01194-3](https://doi.org/10.1038/s43018-026-01194-3). The paper's data-availability statement explicitly links both accessions and the companion [GSE282057](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282057). Published analysis code is linked at [MorrissyLab/GCAR1](https://github.com/MorrissyLab/GCAR1).

### GSE335372 — two-sample xenograft spatial comparison

- **Study and modality:** DBiT spatial transcriptomics on lung sections from an HT29-GL xenograft experiment in NSG mice, with serial H&E sections. The GEO record is organized under human transcriptomic mapping.
- **Biological coverage:** two GEO samples: one CAR-T-treated sample and one CAR-T–drug-conjugate-treated sample, collected 24 hours after treatment. The number of independent mice contributing to each sample is not stated in the Series metadata, so the two GEO records must not be assumed to be replicated biological groups.
- **Access:** merged processed RDS and per-sample processed archives are open in GEO; raw reads are open in SRA/BioProject [PRJNA1478025](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1478025). This is a SubSeries of [GSE335592](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335592).
- **Suitability:** useful only as a preclinical spatial workflow test or descriptive contrast of two conditions.
- **Exact limitations:** likely one spatial specimen per condition unless independent replication is demonstrated elsewhere; treatment and sample are therefore confounded. It is an NSG xenograft rather than a human tumor cohort, uses a single 24-hour endpoint, and has no product or serial blood samples. Human-only mapping may not represent the murine host microenvironment. A CAR-T–drug conjugate is not a chemokine module, and this record cannot establish A3, A4, or A5.
- **Paper:** the GEO record lists no citation. No publication DOI was verified for this accession as of the inventory date; one should not be assigned until the Series or a peer-reviewed data-availability statement links it explicitly.

## Ranked minimal defensible first analysis

### Rank 1 — patient-level product-to-early-blood state analysis

Use GSE197268 as the discovery cohort. Restrict the primary analysis to patients with a valid infusion-product sample and an explicitly attributable post-infusion CAR-T sample. Freeze one high-level state annotation and one primary contrast before examining outcomes, for example the within-patient change in memory-like versus effector/dysfunction-associated fractions from product to day 7. Summarize counts and fractions at patient level; use patient-level pseudobulk or paired estimates rather than treating cells as replicates. Stratify or adjust for axi-cel versus tisa-cel only when the sample size supports it.

This is the smallest analysis that is both clinically relevant and statistically defensible with currently open processed data. It addresses A1 only. It must be described as state dynamics in product and blood, not as chemotaxis, tissue recruitment, retention, or egress.

### Rank 2 — orthogonal early-time-point validation

Use GSE273170 to test the direction and robustness of the prespecified state result in a separate axi-cel subset. Restrict a product-to-blood validation to the patients with both required compartments. Use the available ADTs to support broad T-cell phenotype annotation, while explicitly noting that the panel lacks CXCL16, CXCR6, CXCR3, and CXCR4. Do not pool this cohort with GSE197268 until cohort identity, overlap, batch, and endpoint harmonization have been audited.

### Rank 3 — disease and time-scale sensitivity analyses

Use GSE235760 for a five-patient CAR-positive product-to-expansion-peak sensitivity analysis, and GSE125881 for long-term within-patient CD8 CAR-T trajectories. These analyses test portability across disease, product, sampling window, and sorting strategy; they are not confirmatory replacements for the primary cohort. GSE197215 is a complementary product-only stimulation reference.

### Rank 4 — separate spatial workstream

Use GSE310353 as the primary external PDAC reference for treatment-agnostic ligand-source, tumor-nest, stromal, vascular, and collagen analyses only after the treatment/probe-set metadata conflict is resolved or explicitly retained as a QC limitation. Use GSE282302 Visium for discovery and GSE310352 CosMx for orthogonal cell-level validation; do not treat the two resolutions as independent cohorts. This atlas contains no CAR-T cells.

Use GSE332850 as the CAR-T-treated spatial cohort only after confirming the targeted panel contains the prespecified chemokines, receptors, and anatomical markers and after validating sample-to-file mapping. The primary spatial endpoint should be a patient-level summary, such as the distance from CAR-T cells to ligand-positive sources or tumor nests under a tissue-constrained null. Report both studies as spatial RNA association, not a protein gradient. GSE325706/GSE325772 may provide a single-patient end-to-end case study; GSE335372 is suitable only for pipeline testing or descriptive preclinical illustration.

### Not recoverable from these public accessions

A3 requires non-permeabilized surface measurement of membrane CXCL16 and a separately calibrated soluble-CXCL16 assay. A4 requires raw migration counts or tracks with gradient, uniform-ligand, reverse-gradient, receptor-blockade, viability, and recovery controls. A5 requires donor- or patient-matched products that differ only by the chemokine module. None of these requirements can be reconstructed from the listed expression datasets.

## Repository and redistribution policy

Third-party biological data are inputs, not repository source code.

- Do not commit third-party FASTQ, BAM, CRAM, MTX, H5, H5AD, RDS, parquet, image, source-count, scaled-count, or source cell-metadata files to GitHub.
- Keep `data/raw/` and `data/processed/` excluded by `.gitignore`. Download and processing scripts may create these directories locally, but CI and release jobs must not upload them as caches or artifacts.
- Version only accession metadata, retrieval scripts, transformation code, schemas, configurations, and checksums. For every retrieved file, the local manifest must record the repository accession, exact source URL, UTC retrieval date, byte size, SHA-256 digest, access class, license or terms URL, required citation/DOI, reference version, and the script/version used to retrieve it.
- Do not invent a checksum for a file that has not been retrieved. A public-dataset manifest row becomes complete only after local retrieval has produced the size and SHA-256 value. Repository checksums supplied by the source should be stored in a separate field and labeled with their algorithm.
- Public availability does not automatically grant redistribution rights. Record and follow GEO/SRA, publication, and study-specific terms before sharing any source or derivative file.
- Small, non-sensitive aggregate or figure-source tables may be versioned only after an explicit license, consent, data-use-agreement, and privacy review confirms that redistribution is permitted. Otherwise, deposit permitted outputs in an appropriate restricted or DOI-bearing archive and version only their accession and provenance. A publication release should ideally archive the code and permitted source-data tables with a versioned Zenodo DOI or an equivalent long-term repository record.
- Controlled data, including dbGaP, EGA, GSA-Human, Science Data Bank controlled records, and any derivative covered by a data-use agreement, must never be uploaded to GitHub, a public release, CI, an external artifact store, or a shared container image. Store only the accession, access-status metadata, approved local path, and compliant analysis scripts. Do not commit credentials, approval documents, participant metadata, or access tokens.
- A reproducible public release should fail clearly when an external input is absent and print the accession and expected checksum. It must not silently substitute synthetic data. The deterministic demonstration dataset remains explicitly labeled `synthetic` and is never cited or used as scientific evidence.

## Minimum pre-analysis audit

Before any accession is added to an active analysis manifest:

1. resolve one row per actual sample/file unit with `patient_id`, `sample_id`, `section_id` where applicable, `timepoint_id`, compartment, product/arm, modality, batch, and source accession;
2. verify that patient identifiers are accession-local and do not collide across studies;
3. distinguish biological specimens from modality-specific GEO records;
4. verify raw versus processed status from the sample record, not from a filename such as `RAW.tar`;
5. map every downloaded file to an exact source URL and checksum;
6. identify missing and unbalanced compartments before defining paired contrasts;
7. confirm CAR-positive assignment, feature/reference version, and state-annotation inputs;
8. for targeted spatial panels, verify every prespecified gene before running the workflow;
9. preserve the patient or independent model as the unit of inference; and
10. record corrections, retractions, superseded repository versions, and controlled-access obligations.
