# Chemokines as Regulators of CAR-T Cell Therapy

Reproducible analysis framework for testing how chemokine signals influence CAR-T-cell state, trafficking, spatial positioning, retention, and egress. The repository is designed for five linked workstreams:

1. single-cell RNA-seq/CITE-seq of the infusion product and serial blood samples;
2. spatial mapping of chemokine-producing cells relative to vessels, stroma, tumor nests, and CAR-T cells;
3. protein-level validation of chemokine gradients, with explicit separation of soluble and membrane-associated CXCL16;
4. functional chemotaxis, retention, and egress assays;
5. paired comparison of a chemokine-modified CAR-T product with an otherwise identical CAR-T control.

## Project status

This is a development scaffold on `draft/reproducibility-v0.1`. It contains analysis code, data contracts, configuration, a synthetic demonstration dataset, and reporting templates. The root demonstration remains non-evidentiary and is used only for software testing. A separate checksum-pinned module now reports an exploratory secondary reanalysis of the public GSE125881 and GSE269379 datasets.

No new primary data or controlled-access patient files are committed. Third-party cell-level matrices, metadata, and tissue images remain in their source repositories and are retrieved by verified download commands.

The framework will support a public release only after the input datasets, accession identifiers, checksums, analysis configuration, quality-control decisions, and frozen software environment are available and independently verified.

### Public-data reanalysis

The reviewed analysis is isolated in [`analyses/public_reanalysis/`](analyses/public_reanalysis/). It contains exact source URLs and SHA-256 checksums, a patient/day crosswalk, a pinned environment, a deterministic Python workflow, tests, aggregate result tables, and Figure 6 in PNG, PDF, and SVG formats. The analysis keeps donor or patient as the biological unit and treats the single Visium HD specimen as a descriptive case.

From the module directory, reproduce and verify the frozen outputs on Linux with:

```bash
mamba env create --file environment.yml
mamba activate cart-public-reanalysis
make fetch
make analysis
make verify
make test
```

The module README documents every endpoint and the limits of interpretation, including the IIH comparator, the single spatial case, sparse CXCR6 transcript detection, and the inability of RNA data to distinguish soluble from membrane-bound CXCL16.

## Scientific scope

Chemokine abundance alone does not establish a directionally active gradient. Interpretation requires ligand source, spatial presentation, receptor surface abundance, proteolytic processing, tissue architecture, and time. This is particularly important for CXCL16, which can function as a transmembrane adhesion/scavenger molecule or as a soluble CXCR6 ligand after proteolytic shedding. The workflow therefore keeps RNA, surface protein, soluble protein, spatial localization, and migration readouts as distinct measurements.

The paired comparison is defined at the biological-unit level. Whenever possible, starting material from the same donor or patient is split into two manufacturing arms: chemokine-modified CAR-T and an otherwise identical control. Technical wells, fields of view, regions of interest, and single cells are nested measurements; they are never treated as independent patients.

## Repository layout

```text
.
├── config/                  # analysis configuration and input manifest
├── data/
│   └── demo/               # synthetic inputs for software testing
├── docs/                    # analysis, data, validation, and reproducibility plans
├── resources/               # machine-readable TSV input contracts
├── scripts/                 # reproducibility audit utilities
├── src/chemokine_cart/      # analysis modules and command-line interface
├── workflow/                # Snakemake rules
├── tests/                   # unit and integration tests
├── Makefile                 # reproducible project commands
├── Snakefile                # workflow DAG and declared outputs
├── environment.yml          # Conda bootstrap using the Python lock
├── requirements.lock.txt    # hash-locked Python dependency graph
├── pyproject.toml           # Python package definition
└── CITATION.cff
```

## Quick start on Linux

The commands below run the synthetic demonstration and do not download controlled patient data.

```bash
git clone --branch draft/reproducibility-v0.1 \
  https://github.com/alexeysakhalin/Chemokines-as-Regulators-of-CAR-T-Cell.git
cd Chemokines-as-Regulators-of-CAR-T-Cell

mamba env create --file environment.yml
mamba activate chemokine-cart
python -m pip install --no-deps --no-build-isolation --editable .

chemokine-cart make-demo --output-dir data/demo --seed 20260901 --force
chemokine-cart validate \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --check-files
chemokine-cart run \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901
```

For a hash-locked pip installation on Linux with Python 3.11.11:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip==25.0.1
python -m pip install --require-hashes --requirement requirements.lock.txt
python -m pip install --no-deps --no-build-isolation --editable .
```

Equivalent stage-oriented execution is available through the Python module:

```bash
python -m chemokine_cart.pipeline \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --stage all \
  --random-seed 20260901
```

Available stages are `validate`, `single-cell`, `spatial`, `functional`, `figures`, `report`, and `all`. Before running real data, use a dry run:

```bash
chemokine-cart run \
  --config config/config.yaml \
  --manifest config/manifest.tsv \
  --output-dir results \
  --random-seed 20260901 \
  --dry-run
```

`--random-seed` overrides `project.random_seed` for direct CLI execution. Make targets pass
`SEED`, and Snakemake receives the same value through `--config random_seed=...`. Snakemake
accepts exactly one complete analysis configuration through `--configfile`; partial
configuration overlays are not supported.

### Current analytical contract

- The single-cell module summarizes upstream-assigned `cell_state`, T-cell subtype, and CAR-origin labels; keeps sample, arm, compartment, and batch provenance; reports low-cell and missing-timepoint QC; computes descriptive patient-level mean and paired-change bootstrap intervals; and plots only observed patient trajectories that pass the configured cell threshold, in an explicit timepoint order. Strata without an eligible sample remain in the bootstrap table with `no_eligible_patients` status and blank estimates. RNA and ADT raw counts are never combined: state-specific RNA and ADT pseudobulk tables are written separately from either canonical long fields (`feature_id`, `feature_type`, `count`) or explicitly typed wide count columns. The module does not infer cell states or CAR identity from expression.
- The spatial module validates coordinates, calculates distances to vascular, stromal, tumor, and ligand-source landmarks, summarizes edge-bearing sections with equal section weight per patient, and calculates patient-equal neighborhood enrichment with false-discovery-rate correction under a patient/section-blocked label-permutation null. Its p/q values test conditional spatial association in the observed tissues, not a population-level patient effect. Real-data confirmatory inference additionally requires validated compartment or tissue-mask strata, prespecified minimum graph support per section, and a separate patient-level model; those steps are not yet implemented by the runner.
- The functional module aggregates technical replicates within independent biological replicates, compares paired products, and generates dose-response, migration, and paired-comparison figures. Soluble CXCL16 concentration and membrane CXCL16 fluorescence are summarized separately unless a common calibrated unit and assay have been explicitly validated.

The pipeline writes single-cell state fractions, mean and paired-change confidence intervals, separate RNA/ADT pseudobulk tables, sample QC, and missingness tables; spatial distance and patient-equal neighborhood tables; and functional/protein biological summaries and QC files. Exact filenames are declared in `Snakefile`. Figure 6 panels are exported as PDF, SVG, and 300-dpi PNG together with `figure_manifest.tsv`. The synthetic-only demonstration is allowed but explicitly marked non-evidentiary. With the default evidence policy, any other analysis stops before scientific stages when even one manifest record has `evidence_eligible=false`.

The draft does not yet process raw FASTQ files, call cells, normalize CITE-seq, segment/register raw tissue images, or extract tracks from microscope files. Those upstream operations must be supplied as versioned, checksum-verified inputs until dedicated workflow stages are implemented and tested.

## Required inputs

The manifest is the single source of truth for sample identity and pairing. At minimum it records a stable sample identifier, patient or donor identifier, experimental arm, time point, compartment, assay type, batch, input path, file checksum, and reference-build or assay-panel version. Spatial experiments additionally require image, coordinate, segmentation, region-of-interest, and tissue-annotation metadata. Functional experiments require plate maps or device/run identifiers and raw event- or track-level measurements.

Raw data are not committed to Git. Controlled-access files remain in their authorized repository or institutional storage, and the manifest records their accession and SHA-256 digest. See [Data requirements](docs/DATA_REQUIREMENTS.md).

### Public and third-party datasets

The repository does not redistribute third-party raw counts, FASTQ files, clinical tables, or tissue images. Instead, [`config/public_datasets.tsv`](config/public_datasets.tsv) records an exact repository URL, retrieval date, expected byte size, SHA-256 digest, access class, source terms, primary citation, and ignored local destination for each file that has been independently verified. [`docs/PUBLIC_DATASETS.md`](docs/PUBLIC_DATASETS.md) separates what each accession can support from what it cannot establish.

List, retrieve, and verify the currently checksum-pinned public inputs with:

```bash
make data-list
make data-fetch DATASET=GSE125881
make data-verify DATASET=GSE125881
```

Downloads are written below `data/raw/`, which is ignored by Git. The downloader refuses path traversal, does not replace a mismatching local file without explicit review, and verifies both size and SHA-256 before accepting a download. Controlled-access records are never downloaded by this command and must be obtained under the originating repository's authorization and data-use agreement.

An external dataset is cited by both its primary publication and persistent accession. The Methods and Data Availability sections must identify the exact files or accession release, retrieval date, preprocessing decisions, reference genome or panel, exclusion rules, and software versions. The independent patient or donor remains the unit of inference; cells, spots, fields, and technical replicates are nested observations rather than additional biological replicates. Only small, non-sensitive derived tables may be released when the source license, consent, privacy conditions, and repository terms permit redistribution.

## Analysis principles

- The patient or independent donor is the unit of inference.
- Modified and control products are paired only when they originate from the same biological starting material and have matched manufacturing conditions.
- The current runner writes patient-level pseudobulk count tables for downstream differential-expression modelling; it does not yet fit that model. Future differential-expression inference must use patient-level pseudobulk or an explicitly hierarchical model rather than treating cells as independent replicates.
- Bootstrap intervals emitted by the current draft are descriptive patient-level summaries. Confirmatory cell-state comparisons require the prespecified compositional or binomial/beta-binomial model described in the analysis plan; that model is not yet implemented by this runner.
- Spatial endpoints are summarized per patient or independent specimen before cohort inference. Confirmatory real-data analyses must use tissue-compartment- or mask-constrained spatial null models; the current runner implements only patient/section blocking.
- Multiplicity is controlled within prespecified endpoint families using the Benjamini-Hochberg false-discovery rate.
- Missingness, sample attrition, excluded regions, and failed assays are reported explicitly; no result is generated from absent raw data.

The complete endpoint hierarchy and model formulas are described in [Analysis plan](docs/ANALYSIS_PLAN.md).

## Experimental validation

The computational workflow cannot establish chemokine function by itself. The validation plan includes:

- non-permeabilized flow cytometry or quantitative imaging for membrane CXCL16;
- immunoassay or targeted mass spectrometry of cell-free supernatant for soluble CXCL16;
- gradient, uniform-ligand, reverse-gradient, receptor-blockade, and viability controls in migration assays;
- transendothelial and three-dimensional models for retention and egress;
- matched cytotoxicity, phenotype, expansion, and cytokine-release measurements for the modified and control products.

See [Experimental validation](docs/EXPERIMENTAL_VALIDATION.md).

## Reproducibility and provenance

Every complete run through the `report` stage writes a machine-readable provenance record containing the Git commit and dirty-worktree status, configuration and manifest file records, source-input checksums, the process arguments and a canonical reproduction command, the effective random seed, installed package versions, SHA-256 records for `requirements.lock.txt` and `environment.yml`, and generated-output checksums. Logs, workflow markers, and previously generated provenance files are excluded from scientific-output checksums. A result is considered reproducible only when these records resolve to immutable inputs and the workflow passes from a clean Linux environment.

See [Reproducibility](docs/REPRODUCIBILITY.md).

## Evidence base

The design is informed by evidence that chemokine effects are source-, form-, and context-dependent, and by patient-level CAR-T single-cell studies:

- Ozga AJ, Chow MT, Luster AD. *Immunity* (2021). DOI: [10.1016/j.immuni.2021.01.012](https://doi.org/10.1016/j.immuni.2021.01.012)
- Foeng J, Comerford I, McColl SR. *Cell Reports Medicine* (2022). DOI: [10.1016/j.xcrm.2022.100543](https://doi.org/10.1016/j.xcrm.2022.100543)
- Abel S et al. *Journal of Immunology* (2004). DOI: [10.4049/jimmunol.172.10.6362](https://doi.org/10.4049/jimmunol.172.10.6362)
- Lesch S et al. *Nature Biomedical Engineering* (2021). DOI: [10.1038/s41551-021-00737-6](https://doi.org/10.1038/s41551-021-00737-6)
- Steele MM et al. *Nature Immunology* (2023). DOI: [10.1038/s41590-023-01443-y](https://doi.org/10.1038/s41590-023-01443-y)
- Deng Q et al. *Nature Medicine* (2020). DOI: [10.1038/s41591-020-1061-7](https://doi.org/10.1038/s41591-020-1061-7)
- Sheih A et al. *Nature Communications* (2020). DOI: [10.1038/s41467-019-13880-1](https://doi.org/10.1038/s41467-019-13880-1)
- Sarén T et al. *Clinical Cancer Research* (2023). DOI: [10.1158/1078-0432.CCR-23-0178](https://doi.org/10.1158/1078-0432.CCR-23-0178)
- Steffin D et al. *Nature* (2025). DOI: [10.1038/s41586-024-08261-8](https://doi.org/10.1038/s41586-024-08261-8)

These references justify the measurement strategy; they do not substitute for validation in the present experimental system.

## Documentation

- [Analysis plan](docs/ANALYSIS_PLAN.md)
- [Data requirements](docs/DATA_REQUIREMENTS.md)
- [Experimental validation](docs/EXPERIMENTAL_VALIDATION.md)
- [Reproducibility](docs/REPRODUCIBILITY.md)
- [BioRender figure prompts](docs/BIORENDER_PROMPTS.md)

## Citation and license

Citation metadata are provided in [CITATION.cff](CITATION.cff). Code is released under the [MIT License](LICENSE). Data remain subject to the terms of their originating studies, repositories, consent documents, and data-use agreements.
