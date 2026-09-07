# Analysis plan

## Purpose and inferential boundary

This document prespecifies analyses for evaluating chemokine-dependent CAR-T-cell trafficking and state. It is an analysis plan, not a report of completed experiments. No biological conclusion can be drawn until the required raw data, sample metadata, quality-control results, and provenance records are available.

The five linked workstreams are:

1. scRNA-seq/CITE-seq of the final infusion product and longitudinal blood samples;
2. spatial mapping of ligand sources, vessels, stromal boundaries, tumor nests, and CAR-T cells;
3. orthogonal protein measurement, including soluble and membrane-associated CXCL16;
4. chemotaxis, retention, and egress assays;
5. paired evaluation of a chemokine-modified product against an otherwise identical CAR-T control.

## Study design

### Biological units

The independent biological unit is the patient for clinical analyses and the independent donor or independently established biological model for preclinical analyses. A cell, sequencing library, technical well, field of view, image tile, region of interest, organoid, or repeated measurement is not an independent biological unit when nested within the same patient or donor.

### Paired product comparison

The preferred preclinical design is a split-product experiment. The same leukapheresis or donor starting material is divided before engineering into:

- `modified`: the CAR-T product containing the prespecified chemokine receptor, chemokine ligand, or regulatory module;
- `control`: an otherwise identical CAR-T product lacking only that module.

The following must be matched or recorded: CAR antigen specificity; scFv; hinge, transmembrane, costimulatory, and CD3ζ domains; promoter; vector backbone; manufacturing protocol; activation and cytokine conditions; culture duration; target vector-copy-number range; cell dose; viability threshold; tumor target; and assay batch. A mock-transduced or untransduced T-cell arm may answer additional questions but cannot replace the identical-CAR control.

If modified and control products are administered to different patients, the comparison is not paired. Clinical arm effects then require randomization or a justified adjustment set and remain vulnerable to unmeasured confounding. Paired pre/post samples from the same patient estimate within-patient temporal change, not the causal effect of a product modification.

### Time points and compartments

Time points are protocol-defined before analysis. A typical design includes the final product, pre-infusion baseline, early expansion, contraction, and later persistence windows, but nominal days are not silently pooled. Each sample records actual collection time relative to infusion. Blood, bone marrow, tumor, draining lymph node, and other tissues are analyzed as distinct compartments.

### Endpoint hierarchy

One primary endpoint must be selected before unblinding based on the intended mechanism of the module. Recommended primary endpoints are:

- trafficking module: patient- or donor-level difference in CAR-T density within antigen-positive tumor nests;
- product-state module: patient-level difference in the prespecified memory-like CAR-T fraction in the infusion product;
- retention module: donor-level difference in residence-time area under the curve in the validated retention assay;
- egress module: donor-level difference in the fraction of viable CAR-T cells exiting the tumor-side compartment per unit time.

All other endpoints are secondary or exploratory. Changing the primary endpoint after inspecting the data must be recorded as exploratory.

## Workstream 1: scRNA-seq/CITE-seq

### Input and preprocessing

The complete study requires raw FASTQ files or repository-generated unfiltered count matrices, feature-barcode matrices for antibody-derived tags, sample-level metadata, and—when available—paired V(D)J data and a validated CAR-transgene capture feature. The root package begins from a validated annotated long table and does not perform FASTQ alignment, ambient RNA correction, cell calling, doublet detection, CITE-seq normalization, or V(D)J assembly. Their external commands, thresholds, software versions, reference bundles, and output checksums must be supplied as upstream provenance. Samples are not excluded solely because their biology differs from the cohort.

RNA and antibody modalities are quality-controlled separately. Antibody-derived tags are normalized with a method appropriate to the experimental controls, such as centered log-ratio normalization or a background-aware method when empty droplets and isotype controls are available. Integration is used for visualization and annotation, not as a substitute for replicate-aware inference.

The current implementation expects `cell_state` and CAR-versus-endogenous `t_cell_origin` to be assigned upstream. It summarizes those labels and does not derive them from gene expression. The annotation file, CAR gate, marker/reference set, and version must therefore be supplied as inputs.

### Cell identity and CAR-T assignment

CAR-T identity is assigned using a validated CAR-transgene feature, construct-specific oligonucleotide tag, or an experimentally validated surrogate. CD3 expression alone does not distinguish infused CAR-T cells from endogenous T cells. When CAR capture is unavailable, analyses are labelled as total T-cell analyses unless clonotype or another validated feature permits attribution.

Cell-state annotation combines canonical markers, protein features, reference mapping, and manual review. The prespecified high-level states are:

- naïve/stem-like or central-memory-like;
- effector-memory/effector;
- proliferating;
- interferon-responsive;
- dysfunctional/exhaustion-associated;
- senescence-associated.

Exhaustion and senescence are not treated as synonyms. A state label requires a coherent multigene or multimodal program and cannot be assigned from a single checkpoint or cell-cycle marker. Marker sets, reference version, scoring direction, and cutoffs are frozen in the run configuration.

### Primary single-cell endpoints

- fraction of CAR-positive CD8 and CD4 cells in each prespecified state;
- patient-level pseudobulk expression of the prespecified state program;
- expansion and contraction of CAR-associated clonotypes across serial time points;
- surface abundance of relevant chemokine receptors measured by CITE-seq or flow cytometry;
- within-patient transition in state composition relative to the infusion product.

### Statistical analysis

Differential expression uses patient-level pseudobulk counts within a prespecified cell type or state. A typical paired model is:

```text
expression ~ patient_id + experimental_arm + timepoint_id + experimental_arm:timepoint_id + batch_id
```

Terms are included only when supported by the design matrix. For repeated serial observations, an equivalent generalized linear mixed model may include a patient random intercept. Cell counts are offsets where appropriate.

Cell-state proportions are analyzed with a compositional model or a binomial/beta-binomial mixed model. The denominator is stated explicitly. A typical model is:

```text
state_count / total_CAR_T ~ experimental_arm * timepoint_id + covariates + (1 | patient_id)
```

Clonal expansion is summarized per patient using preregistered diversity and dominance measures. Clone-level descriptive plots do not replace patient-level inference. Differential-expression and marker-testing families use Benjamini-Hochberg false-discovery-rate control. Effect sizes and confidence intervals are reported with adjusted P values.

Single-cell differential-expression methods can otherwise produce false precision when cells are treated as replicates; the replicate-aware strategy follows the evidence summarized by Squair et al., *Nature Communications* (2021), DOI: [10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2).

## Workstream 2: spatial mapping

### Required annotations

Each specimen requires registered expression or imaging data, tissue mask, segmentation, spatial coordinates, scale, and blinded annotations. The root package begins from a segmented object/spot table and does not segment or register raw images. Those upstream transformations must therefore be versioned and checksum-verified. Required annotations include:

- perfused and non-perfused vessels when perfusion information is available;
- endothelial and lymphatic endothelial compartments;
- fibroblast-rich stroma and extracellular-matrix boundaries;
- antigen-positive and antigen-negative tumor nests;
- necrotic or non-evaluable regions;
- endogenous T cells and CAR-T cells, when distinguishable;
- candidate chemokine-producing cells.

RNA abundance is not called a protein gradient. A functional gradient requires orthogonal spatial protein or functional evidence.

### External PDAC reference atlas

GSE310353 is designated as a treatment-agnostic external reference for the non-CAR-T portion of this workstream. Its GSE282302 Visium series is used for whole-transcriptome discovery across tumor nests, stromal regions, vascular regions, and collagen-registered tissue; its GSE310352 CosMx series is used as orthogonal single-cell spatial validation. The Visium and CosMx layers are linked parts of one study and are not counted as independent cohorts. The atlas contains no CAR-T cells, infusion product, or serial blood.

The CosMx feature audit confirmed `CXCL16`, `CXCR6`, `CXCL9`, `CXCL10`, `CXCL12`, `CXCL14`, `CXCR4`, `CCL2`, `CCL19`, `CCL21`, and `CCR7`; `CXCL11` is absent. Feature presence and detection rate are checked per imported file before any endpoint is frozen. CXCL16 is interpreted only as transcript expression from a candidate source compartment. The dataset does not separate membrane and soluble CXCL16 and does not measure a protein gradient.

The paper describes the human tumors as untreated/treatment-naive and patient Visium as human transcriptome probe set v1, whereas current GSE282302 Series and sample metadata describe post-neoadjuvant FOLFIRINOX tissue and probe set v2. This conflict is a QC hold. No treatment-status claim, treatment comparison, or confirmatory biological conclusion is made until an author or GEO correction resolves treatment and reference version at sample level. A treatment-agnostic exploratory run may proceed only with the discrepancy carried into the provenance record and report.

The published concepts to adapt are multi-resolution Visium discovery followed by CosMx validation, tissue registration, concentric tumor-distance bins, and cell-neighborhood graphs. The implementation does not copy the original inferential unit. Ligand-source abundance, distance, and neighborhood endpoints are summarized per section and patient; cohort inference uses patient-level summaries or mixed models with sections and regions nested within patient. The current runner permutes neighborhood labels within patient and section and uses a section-then-patient equal-weight statistic over sections containing at least one graph edge. Its p/q values describe a conditional random-label test in the observed tissues and are not a population-level patient-effect test. It reports total and edge-bearing section/patient counts but does not yet enforce a prespecified minimum number of cells or graph edges per section. Before real-data inference, the canonical spatial table and permutation strata must be extended with validated tissue-compartment or mask labels so that shuffling also preserves compartment boundaries and local density, minimum graph-support rules must be frozen, and a prespecified patient-level model must be fitted; unrestricted cross-compartment shuffling is not considered confirmatory.

### Spatial endpoints

- CAR-T density per mm² in viable antigen-positive tumor nests;
- fraction of CAR-T cells located in tumor nests, stroma, perivascular zones, and invasive margins;
- nearest-neighbor distance from each CAR-T cell to a chemokine-positive source, perfused vessel, stromal boundary, and antigen-positive tumor cell;
- crossing ratio across the stromal-tumor boundary;
- enrichment of CAR-T cells around ligand-positive sources relative to a tissue-constrained null distribution;
- spatial association between membrane CXCL16, soluble-CXCL16 proxy measurements, CXCR6-positive CAR-T cells, and compartment identity;
- distance of CXCR4-positive CAR-T cells to lymphatic vessels and observed exit structures when egress is being tested.

Distances are reported in micrometres, not pixels. Regions with failed segmentation, folds, necrosis, saturation, or insufficient marker quality are excluded using blinded rules and retained in the exclusion log.

### Spatial statistics

Cell-level distances are summarized per region and then per independent specimen. Cohort inference is performed at the patient or donor level. Hierarchical models may retain region-level information with random intercepts for patient and specimen. For confirmatory real-data inference, permutation tests must preserve the tissue mask, compartment areas, and observed cell density; unrestricted shuffling across non-tissue space is invalid. The current runner does not yet implement those mask- or compartment-constrained permutations.

A typical paired model is:

```text
endpoint ~ experimental_arm * timepoint_id + compartment + batch_id + (1 | patient_id/specimen_id)
```

Spatial autocorrelation is measured and reported. Confidence intervals are obtained by patient-level bootstrap or a justified hierarchical bootstrap, not by resampling individual cells as if independent. When multiple distances or ligand-receptor pairs are tested, false-discovery-rate control is applied within the prespecified spatial family.

Spatial interpretation is motivated by the context dependence of CXCL16-CXCR6 signaling reported by Chia et al., *Frontiers in Immunology* (2023), DOI: [10.3389/fimmu.2023.1331287](https://doi.org/10.3389/fimmu.2023.1331287), and by the spatial confinement of IFNγ signaling described by Centofanti et al., *PNAS* (2023), DOI: [10.1073/pnas.2304190120](https://doi.org/10.1073/pnas.2304190120).

## Workstream 3: protein validation and CXCL16 proteoforms

Membrane and soluble CXCL16 are separate analytes:

- membrane CXCL16: surface staining on viable, non-permeabilized cells by flow cytometry or quantitative tissue imaging;
- soluble CXCL16: validated immunoassay or targeted mass spectrometry in clarified cell-free supernatant, plasma, or tissue fluid;
- total CXCL16: lysate or permeabilized-cell measurement, reported separately and not used as a surrogate for either surface or soluble CXCL16.

Primary measurements are surface molecules or fluorescence intensity normalized to a validated reference and soluble concentration normalized to sample volume and viable producer-cell count where appropriate. Soluble concentration and membrane fluorescence must not be tested against each other as if they shared a scale. Direct form-to-form testing is allowed only within one validated assay and a common calibrated unit, with that comparability explicitly declared. Otherwise, each form is summarized and compared between products separately. The soluble-to-surface ratio is exploratory because it is affected by expression, cell number, death, internalization, and cleavage. ADAM10 perturbation can support a shedding mechanism only when viability, total CXCL16, and assay recovery are controlled. The mechanistic basis for ADAM10-dependent shedding is supported by Abel et al., *Journal of Immunology* (2004), DOI: [10.4049/jimmunol.172.10.6362](https://doi.org/10.4049/jimmunol.172.10.6362).

Protein endpoints are analyzed at the patient/donor level using paired models for split products and repeated-measures models for serial samples. Plate, operator, and acquisition batch are recorded. Values below quantification are handled using a prespecified censored-data rule; they are not replaced by zero.

## Workstream 4: functional trafficking assays

### Chemotaxis

The minimum design includes no-ligand, forward-gradient, uniform-ligand, reverse-gradient, receptor-blockade or receptor-deficient, and viability controls. A forward gradient is compared with uniform ligand to distinguish chemotaxis from chemokinesis. Primary endpoints are migrated fraction and, for live imaging, chemotactic index and velocity. Total recovery is reported to distinguish migration from loss or death.

### Retention

Retention is evaluated after an explicit loading phase followed by controlled washout, perfusion, or ligand withdrawal. Primary endpoints are retained viable fraction at a prespecified time and residence-time area under the curve. Adhesion under flow and persistence within a three-dimensional tumor/stromal model answer different questions and are reported separately.

### Egress

Egress requires an anatomically directed model, such as migration across a lymphatic endothelial barrier or exit from a three-dimensional tumor compartment into a defined collecting compartment. Primary endpoints are viable egress fraction per unit time and time to exit. Low recovery cannot be interpreted as retention unless death, nonspecific adhesion, sampling loss, and barrier failure are excluded.

The CXCL12-CXCR4 axis can contribute to entry, positioning, retention, or lymphatic egress depending on context. The egress analysis therefore does not assign a universal direction to this axis. Steele et al. demonstrated antigen- and CXCR4-dependent regulation of T-cell egress in a preclinical melanoma setting (*Nature Immunology*, 2023; DOI: [10.1038/s41590-023-01443-y](https://doi.org/10.1038/s41590-023-01443-y)); extrapolation to a CAR-T product requires direct testing.

### Functional statistical model

Technical replicates are summarized within each independent donor and run after prespecified quality control. A typical paired model is:

```text
endpoint ~ experimental_arm * ligand_condition + run_order + (1 | patient_id)
```

Concentration-response experiments use a prespecified nonlinear model when supported by the observed range. The number of donors is chosen from the minimal biologically relevant effect and pilot-derived between-donor variance. Adding technical wells does not increase biological sample size.

## Workstream 5: integrated modified-versus-control assessment

The integration is mechanistic, not a single unvalidated composite score. Evidence is considered concordant when the modified product shows the prespecified directional effect in orthogonal measurements, for example receptor surface expression, gradient-dependent migration, and intratumoral localization, without unacceptable loss of viability, expansion, target-specific killing, or phenotypic stability.

The following safety and product-quality endpoints are retained alongside trafficking endpoints:

- viability, expansion, transduction, and vector-copy-number distribution;
- CAR surface expression and antigen-specific cytotoxicity;
- cytokine release under antigen-positive and antigen-negative conditions;
- tonic activation and activation-induced cell death;
- off-tumor migration toward ligand-positive normal-cell or tissue models when biologically relevant;
- memory/effector/dysfunctional-state composition.

An increase in tissue entry is not considered beneficial if it is accompanied by loss of antigen specificity, excessive antigen-independent activation, impaired persistence, or redistribution to an unsafe compartment.

## Covariates, missingness, and multiplicity

Prespecified covariates may include disease type, tumor site, baseline burden, prior therapy, lymphodepletion, manufacturing batch, sequencing batch, tissue compartment, and actual sampling day. Covariates are limited by the number of independent biological units; high-dimensional adjustment is not used in a small cohort.

Reasons for missing samples and failed assays are categorized before outcome analysis. No last-observation-carried-forward imputation is used. Primary analyses use observed data under the stated model; sensitivity analyses address informative attrition when the data support them.

Endpoint families are defined as single-cell state, pseudobulk expression, spatial localization, protein, chemotaxis, retention, egress, and safety. False-discovery-rate correction is performed within each family. Primary endpoints retain a prespecified two-sided significance threshold, while secondary and exploratory findings are interpreted through effect sizes, uncertainty, and multiplicity-adjusted evidence.

## Power and reporting

Power is based on independent patients or donors and the variance of the patient-level endpoint. Pilot data are used to simulate the planned paired or mixed-effects model. Cell number, image tiles, or technical replicates cannot compensate for too few independent patients.

Every report includes:

- biological and technical sample counts at each processing step;
- exclusions and reasons;
- effect size with 95% confidence interval;
- raw and adjusted P values where applicable;
- model formula and contrast;
- software, reference, manifest, and configuration versions;
- patient-level or donor-level points in quantitative plots;
- explicit labels for confirmatory and exploratory analyses.

## Limitations before data acquisition

The prospective protein, functional, and matched-product workstreams require additional experimental inputs. The separate [`public reanalysis`](../analyses/public_reanalysis/README.md) reports patient-level transcript-detection summaries from deposited counts and author annotations, but does not establish a clinical protein gradient, CXCL16 proteoforms, retention, egress, or the effect of adding a chemokine module. Raw-sequence processing, cell annotation, image registration/segmentation, and microscopy-track extraction remain upstream responsibilities. Public expression datasets cannot replace a matched experiment if they lack the same CAR construct, tissue context, time points, and paired control.

## Key references

- Deng Q et al. CAR-T infusion-product states associated with efficacy and toxicity. *Nature Medicine* (2020). DOI: [10.1038/s41591-020-1061-7](https://doi.org/10.1038/s41591-020-1061-7)
- Sheih A et al. Longitudinal clonal kinetics and single-cell profiling of CAR-T cells. *Nature Communications* (2020). DOI: [10.1038/s41467-019-13880-1](https://doi.org/10.1038/s41467-019-13880-1)
- Sarén T et al. Patient-level CAR-T single-cell states associated with response. *Clinical Cancer Research* (2023). DOI: [10.1158/1078-0432.CCR-23-0178](https://doi.org/10.1158/1078-0432.CCR-23-0178)
- Steffin D et al. Paired blood/tumor profiling of IL-15-armored GPC3 CAR-T cells. *Nature* (2025). DOI: [10.1038/s41586-024-08261-8](https://doi.org/10.1038/s41586-024-08261-8)
- Lesch S et al. CXCR6 engineering in adoptive T-cell therapy. *Nature Biomedical Engineering* (2021). DOI: [10.1038/s41551-021-00737-6](https://doi.org/10.1038/s41551-021-00737-6)
- Foeng J et al. Chemokine-system engineering for CAR-T homing. *Cell Reports Medicine* (2022). DOI: [10.1016/j.xcrm.2022.100543](https://doi.org/10.1016/j.xcrm.2022.100543)
