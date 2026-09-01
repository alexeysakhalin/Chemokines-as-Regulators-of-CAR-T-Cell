# Experimental validation plan

## Objective

The purpose of validation is to determine whether a defined chemokine module changes CAR-T trafficking, positioning, retention, or egress while preserving antigen-dependent function and product quality. Transcript abundance and computational ligand-receptor scores are hypothesis-generating; they do not establish a functional protein gradient or causal migration mechanism.

No experiments described here have been completed as part of the repository. Exact concentrations, sampling times, and sample size must be optimized and frozen before confirmatory testing.

## Matched-product design

### Required comparator

For every independent donor or patient starting material, generate:

- chemokine-modified CAR-T product;
- identical CAR-T control differing only by absence of the chemokine module.

Use the same source material, activation, vector backbone, CAR cassette, multiplicity target, culture medium, cytokines, vessels, culture duration, harvest criterion, cryopreservation, and assay day. Randomize arm position across plates, devices, and imaging chambers. Analysts should remain blinded to product arm until primary quality control and endpoint extraction are locked.

An untransduced T-cell or irrelevant-receptor control can test background behavior but does not replace the identical-CAR control. If the chemokine module changes vector size or expression, quantify CAR abundance and vector copy number and interpret trafficking effects conditional on these differences.

### Release and comparability panel

Before functional comparison, document:

- viable-cell recovery and expansion;
- CD3, CD4, CD8, and CAR surface abundance;
- transduction efficiency and vector copy number where applicable;
- memory/effector phenotype and activation markers;
- sterility, endotoxin, and mycoplasma status for translational experiments;
- antigen-specific killing against target-positive cells;
- background activity against target-negative cells;
- cytokine release with and without cognate antigen.

Large differences in viability, CAR expression, cell dose, or subset composition are mechanistic confounders and must not be attributed to chemokine signaling without adjustment and follow-up experiments.

## CXCL16: soluble and membrane-associated forms

CXCL16 is synthesized as a transmembrane molecule and can be released by proteolytic shedding. Membrane-associated and soluble CXCL16 must therefore be measured independently. ADAM10-dependent constitutive shedding and cytokine-inducible CXCL16 were demonstrated by Abel et al., *Journal of Immunology* (2004), DOI: [10.4049/jimmunol.172.10.6362](https://doi.org/10.4049/jimmunol.172.10.6362).

### Surface CXCL16

Measure surface CXCL16 on viable, non-permeabilized producer cells by flow cytometry. Include fluorescence-minus-one controls, an appropriate negative biological control, viability dye, compensation or spectral-unmixing controls, and a validated positive control. Record antibody clone and lot, staining temperature and duration, acquisition settings, and time from harvest to fixation or acquisition.

Report:

- fraction of CXCL16-positive viable cells;
- median fluorescence intensity with the chosen reference normalization;
- absolute molecules per cell when calibrated beads and a validated quantitative workflow are used;
- producer-cell identity and, for tissue, anatomic compartment.

Internalized or total cellular CXCL16 requires a separate permeabilized assay and cannot be combined with the surface measurement.

### Soluble CXCL16

Measure soluble CXCL16 in clarified conditioned medium, plasma, or tissue fluid using a validated immunoassay or targeted mass-spectrometry assay. Include a full standard curve, blank, matrix spike, dilution linearity, replicate acceptance rule, and lower and upper quantification limits. Conditioned-medium results are normalized to collection volume, conditioning time, and viable producer-cell number.

Do not directly test soluble CXCL16 concentration against membrane CXCL16 fluorescence. They have different units and measurement models. Compare modified and control products within each molecular form. A form-to-form comparison is permitted only after both measurements have been placed on one validated calibrated scale within the same assay.

At collection, record cell viability and a prespecified cell-lysis marker. Centrifugation, filtration, storage, and freeze-thaw conditions must be identical across arms. Values below the quantification limit are retained as censored measurements rather than set to zero.

### Shedding mechanism

To test whether a change reflects shedding rather than altered transcription or cell loss, combine:

- surface CXCL16;
- soluble CXCL16;
- total CXCL16 RNA and/or cellular protein;
- ADAM10 abundance or activity;
- a pharmacological or genetic ADAM10 perturbation;
- matched viability and lysis measurements.

An ADAM10 inhibitor requires a concentration-response window that suppresses target activity without materially reducing viability. A rescue or independent genetic perturbation strengthens causal interpretation. The soluble-to-surface ratio is supportive but not definitive because it is also influenced by synthesis, internalization, cell number, and degradation.

### Spatial protein validation

In tissue, combine CXCL16 with markers for tumor cells, endothelial cells, lymphatic endothelium, myeloid cells, fibroblasts, CXCR6, CAR or transgene, and tumor antigen. Quantify membrane-associated signal at segmented cell boundaries separately from extracellular or diffuse signal. Serial sections and orthogonal assays are used when the imaging platform cannot distinguish extracellular soluble protein from membrane staining.

The context dependence of CXCL16-CXCR6 localization is supported by Chia et al., *Frontiers in Immunology* (2023), DOI: [10.3389/fimmu.2023.1331287](https://doi.org/10.3389/fimmu.2023.1331287). CXCR6 engineering improved adoptive T-cell accumulation and activity in pancreatic tumor models in Lesch et al., *Nature Biomedical Engineering* (2021), DOI: [10.1038/s41551-021-00737-6](https://doi.org/10.1038/s41551-021-00737-6); this does not establish the same effect in another tumor or CAR construct.

## Chemotaxis assays

### Minimum condition matrix

For each product arm and independent donor, include:

1. no ligand in either compartment;
2. forward ligand gradient;
3. the same ligand concentration in both compartments;
4. reverse gradient;
5. forward gradient with receptor blockade or receptor-deficient control;
6. vehicle and isotype controls appropriate to the perturbation;
7. target-positive and target-negative tumor-conditioned medium when testing complex cues.

The uniform-ligand condition distinguishes directional chemotaxis from increased nondirectional motility. A reverse gradient tests directional interpretation. A receptor-blockade condition tests pathway dependence. Ligand concentrations should span a prespecified physiologic and assay-responsive range; a single supraphysiologic concentration is insufficient.

### Transwell or transendothelial assay

Use a validated pore size and membrane coating suitable for activated human T cells. When vascular entry is the question, add a confluent endothelial monolayer and document barrier integrity before and after the experiment. Load equal numbers of viable CAR-T cells. Quantify migrated cells with an absolute-counting method and measure input, migrated, non-migrated, adherent, and dead fractions where feasible.

Primary endpoint:

```text
migrated viable cells / viable input cells
```

Secondary endpoints include total recovery, receptor dependence, and migration relative to the no-ligand and uniform-ligand controls. Report donor-level values; replicate wells are technical replicates.

### Microfluidic live-cell assay

Confirm gradient formation with a fluorescent tracer of similar diffusion behavior or a validated device model. Track cells using a locked segmentation and tracking pipeline. Prespecify minimum track duration and rules for boundary contacts, divisions, and lost tracks.

Report:

- chemotactic index along the gradient axis;
- migration speed;
- persistence and displacement;
- fraction of directionally responding cells;
- viable-cell recovery at the end of imaging.

Do not infer chemotaxis from speed alone.

## Retention assays

Retention is operationally defined as persistence within a specified compartment after entry and after removal or reversal of the recruitment cue.

### Adhesion and washout

Load CAR-T cells onto ligand-presenting tumor, stromal, endothelial, or matrix surfaces under a standardized attachment phase. Apply a defined washout or shear profile. Record the retained viable fraction over time and the distribution of residence times. Include ligand-negative, receptor-blocked, and matrix-only controls.

### Three-dimensional organotypic retention

Use matched three-dimensional tumor/stromal cultures with defined tumor-antigen expression and a documented chemokine source. Image entry and subsequent position after external ligand withdrawal. Measure:

- viable CAR-T density in the model;
- residence-time area under the curve;
- distance to tumor nests, stromal boundaries, and ligand-producing cells;
- target-cell engagement and killing;
- redistribution after receptor blockade or ligand withdrawal.

Persistent signal from nonviable or immobilized cells is not retention. Viability and track continuity are required.

## Egress assays

Egress is distinct from failure to enter. It requires CAR-T cells to begin inside a defined tumor-side or tissue compartment and then cross into a collecting compartment.

### Lymphatic endothelial model

Establish a polarized lymphatic endothelial barrier with validated junctional integrity. Place CAR-T cells in the tissue-side compartment and quantify movement to the lymphatic-side compartment over time. Orient the chemokine condition according to the biological question and include uniform and reverse-gradient controls.

Measure:

- viable egressed fraction per unit time;
- time to barrier crossing by live imaging when feasible;
- barrier integrity and permeability;
- non-specific retention on the membrane or chamber;
- cell death in each compartment;
- CXCR4, CXCR6, and other prespecified receptor abundance before and after the assay.

### Tissue or organoid exit

Preload labelled CAR-T cells into a three-dimensional tumor model, remove unattached cells, and collect cells leaving the model at defined intervals. Quantify viable cells remaining, exiting, dying, and adhering to the external surface. An outward chemokine cue can be varied, but a lower recovered count is not evidence of retention until technical loss and death are excluded.

The role of CXCL12-CXCR4 in lymphatic egress is context-dependent. Steele et al. reported that tumor-associated lymphatic CXCL12 positioned CXCR4-positive CD8 T cells and promoted egress in a preclinical melanoma model (*Nature Immunology*, 2023; DOI: [10.1038/s41590-023-01443-y](https://doi.org/10.1038/s41590-023-01443-y)). The direction and magnitude must be retested in the selected CAR-T and tumor system.

## Spatial mapping in tissue

### Sampling strategy

Obtain paired pretreatment and on-treatment tissue when ethically and clinically feasible. Multisite sampling is preferred because one core cannot represent a spatially heterogeneous lesion. Record anatomic site, biopsy guidance, viable tumor content, necrosis, tissue orientation, fixation, ischemia time, and section order.

### Required compartments and markers

The panel should resolve:

- CAR-T cells versus endogenous T cells;
- antigen-positive tumor cells;
- vascular and lymphatic endothelium;
- fibroblasts and extracellular matrix;
- myeloid cells;
- CXCL16 and CXCR6;
- other module-specific ligand-receptor pairs;
- perfusion or vessel-function marker where the protocol permits.

Spatial transcriptomics can nominate ligand-producing cells but cannot distinguish secreted, glycosaminoglycan-bound, cleaved, and membrane-associated protein. Multiplex protein imaging, in situ hybridization, and soluble-protein assays provide complementary evidence.

### Prespecified quantitative readouts

- CAR-T cells per mm² of viable antigen-positive tumor;
- percentage of CAR-T cells in tumor nests versus stroma and perivascular zones;
- distance from CAR-T cells to ligand-positive cells and perfused vessels;
- stromal-boundary crossing ratio;
- ligand-source/CAR-T colocalization relative to a compartment-preserving spatial null;
- within-patient change from baseline to on-treatment biopsy.

Analysis is blinded to product arm during segmentation and primary endpoint extraction. Regions are excluded using prespecified image-quality rules, and all exclusions remain in the audit table.

## Functional integration with cytotoxicity

Trafficking improvement is interpreted together with antigen-dependent activity. In the same donor-matched products, quantify:

- short-term and serial target-cell killing;
- killing of antigen-negative controls;
- CAR-T expansion and viable recovery;
- degranulation and cytokine release;
- activation-induced cell death;
- phenotype after repeated antigen exposure.

Where possible, test spatially structured models in which tumor nests are separated by stroma or endothelium. This determines whether increased migration produces wider tumor coverage rather than only higher accumulation at the tissue boundary.

## Replication, randomization, and statistics

Biological replication is defined by independent donors, patients, or independently established models. Technical replicate wells or fields estimate assay variability but do not increase the inferential sample size.

Before confirmatory work:

1. define the minimal biologically relevant effect;
2. estimate between-donor variance from pilot data;
3. simulate the paired model and attrition;
4. lock the primary endpoint, exclusion rules, and contrasts;
5. randomize plate, device, and acquisition order;
6. blind endpoint extraction where feasible.

Use paired analyses for split-product designs and hierarchical models for repeated technical measurements. Display all donor-level or patient-level points, effect sizes, and confidence intervals. Apply multiplicity correction within prespecified endpoint families.

## Minimum evidence for a mechanistic claim

A claim that the chemokine module improves trafficking requires all of the following:

- verified module and receptor surface expression;
- a gradient-dependent functional effect distinct from chemokinesis;
- loss or attenuation of the effect with pathway blockade or genetic control;
- an orthogonal spatial or three-dimensional localization readout;
- preservation of viability and antigen-specific CAR function;
- replication across independent biological donors or models.

A claim concerning CXCL16 shedding additionally requires separate surface and soluble measurements with a perturbation consistent with the proposed protease mechanism. A clinical-benefit claim requires an appropriately designed clinical study; ex vivo migration or xenograft localization alone is insufficient.
