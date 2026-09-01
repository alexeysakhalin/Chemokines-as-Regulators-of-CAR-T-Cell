# Data requirements and provenance contract

## Scope

This document defines the minimum information needed to reproduce the analyses. Raw clinical data are not stored in Git. The repository contains only code, schemas, configuration, documentation, and explicitly labelled synthetic demonstration data.

The machine-readable TSV contracts in `resources/` are authoritative for validation. The human-readable requirements below explain the scientific meaning of each field and the additional assay-specific metadata needed for interpretation.

## Core manifest

Each row represents one immutable input object or one assay-level sample, as defined by the schema. Identifiers must be stable, unique within their scope, and free of direct personal identifiers.

### Required schema columns

| Field | Scientific meaning |
|---|---|
| `record_id` | Unique manifest-row identifier; never recycled. |
| `sample_id` | Stable assay-sample identifier. |
| `patient_id` | Pseudonymous independent biological unit; for preclinical work this may identify the independent donor or model. |
| `section_id` | Histological section or spatial unit; use the schema-defined missing value for non-spatial assays. |
| `timepoint_id` | Protocol-defined product, baseline, or post-infusion time point. |
| `modality` | RNA, ADT, V(D)J, spatial RNA, image, protein, or functional modality. |
| `assay` | Specific assay within the modality, such as scRNA-seq, CITE-seq, chemotaxis, retention, or egress. |
| `condition` | Additional experimental condition, such as ligand, stimulation, or perturbation. |
| `comparator` | Prespecified comparator or contrast label. |
| `experimental_arm` | Canonical arm label, including modified and identical-control assignments. |
| `compartment` | Product, blood, marrow, tumor, lymph node, or another prespecified compartment. |
| `batch_id` | Manufacturing, library, imaging, plate, or run batch as appropriate. |
| `reference_version` | Versioned genome, annotation, panel, construct, or assay reference required for interpretation. |
| `source_path` | Project-relative path for an authorized local input; leave empty only when an accession resolves the input. |
| `accession` | Stable repository accession at the most specific available level. |
| `sha256` | Lowercase SHA-256 digest; mandatory for a local source. |
| `data_origin` | Controlled provenance term, including explicit identification of synthetic data. |
| `evidence_eligible` | Whether the record may contribute to scientific evidence; must be `false` for synthetic data. |
| `notes` | Non-identifying structured clarification; not a substitute for required metadata. |

`record_id` must be unique. At least one of `source_path` or `accession` is required. Local paths are project-relative and are validated against the configured project root. A synthetic record is never evidence-eligible.

Each real-data row represents one sample/file unit. The aggregate `MULTI` value is reserved for the synthetic demonstration and is not permitted as a shortcut for real samples. `protocol_version` and `assay_panel_version` are optional schema fields but become scientifically required whenever protocol or panel differences could change the measurement.

### Pairing and assay metadata

The core manifest is accompanied by versioned assay metadata keyed by `sample_id` or `record_id`. Fields needed for a particular analysis include:

| Field | Requirement | Scientific meaning |
|---|---:|---|
| `specimen_id` | recommended | Physical specimen from which libraries, images, or assays were derived. |
| `product_id` | required for product studies | Manufactured CAR-T lot or split-product identifier. |
| `pair_id` | required for paired product analyses | Link between modified and identical-control products from the same starting material. |
| `chemokine_module` | required for product studies | Exact receptor, ligand, or regulatory module; use `none` for the identical-CAR control. |
| `car_construct` | required for product studies | Versioned construct identifier linked to a construct specification. |
| `collection_time` | recommended | Actual time relative to infusion, with unit. |
| `file_size_bytes` | required for provenance | Exact byte size used to detect truncation or replacement. |
| `file_format` | required | FASTQ, MTX, H5, CSV, TSV, OME-TIFF, Parquet, or another validated format. |
| `assay_version` | required | Chemistry, panel, kit, instrument method, or acquisition protocol version. |
| `reference_id` | required for sequence data | Versioned reference bundle identifier. |
| `source_version` | required | Repository release, data-freeze date, or immutable source version. |
| `license_or_dua` | required | Data license or controlled-access/data-use agreement identifier. |

`pair_id` is valid only when both arms originate from the same biological starting material and the CAR backbone and manufacturing conditions are otherwise matched. Similar donors or sequential manufacturing runs are not automatically paired. Do not use a mutable web download URL as the only provenance. Record the accession, repository, object name, version, byte size, and checksum.

## scRNA-seq/CITE-seq and V(D)J

### Accepted primary inputs

Preferred inputs are raw FASTQ files together with sample sheets and the exact reference bundle. Repository-generated unfiltered feature-barcode matrices may be used when FASTQ files are unavailable, but the producing software, version, parameters, and reference must be recorded. Filtered matrices alone are insufficient to reproduce cell calling and ambient-background correction.

Required or conditionally required metadata include:

- library chemistry and kit version;
- sequencing instrument and read structure;
- genome assembly and gene-annotation release;
- feature-barcode reference with antibody clone, fluorophore or oligonucleotide tag, supplier, lot, and panel version;
- CAR-transgene sequence and capture-feature definition;
- sample multiplexing or hashing design and demultiplexing rules;
- V(D)J reference version and contig-level output when clonotypes are analyzed;
- expected-cell recovery and loaded-cell count;
- viability before loading;
- tissue-processing interval and cryopreservation status;
- raw RNA and ADT count layers.

Upstream annotation must provide the `state` and CAR-versus-endogenous `origin` labels used by the current single-cell summarization module. The module does not infer those labels from expression. Each label set must therefore be accompanied by a versioned annotation file, marker/reference definition, CAR-gating method, reviewer, and date.

### Minimum longitudinal design

For serial analyses, the product and at least two post-infusion windows are recommended, together with actual collection times. Missing time points are retained as missing; nearby samples are not relabelled into the same nominal day without a prespecified window. Clinical response, toxicity, lymphodepletion, cell dose, tumor burden, and relevant prior therapies should be supplied in a separate controlled-access clinical table keyed by pseudonymous patient identifier.

## Spatial transcriptomics and multiplex imaging

Each spatial sample requires:

- raw expression counts or raw image channels;
- full-resolution image with physical pixel size;
- spot, cell, or molecule coordinates in the image coordinate system;
- transformation matrices for registered modalities;
- segmentation masks and segmentation-software version;
- tissue and quality masks;
- region-of-interest definitions with inclusion/exclusion status;
- marker panel, clone, lot, cycle, exposure, and channel assignment;
- antigen, CAR/transgene, vascular, lymphatic, stromal, tumor, viability, and candidate-chemokine markers where applicable;
- annotation provenance and blinded-review status;
- section thickness, orientation, anatomic site, and distance between serial sections;
- perfusion marker and acquisition timing if claims concern perfused vessels.

Spatial RNA and spatial protein measurements must not share the same analyte label without a `modality` field. Pixel distances must be convertible to micrometres. Cropped images are insufficient unless the uncropped source image and transformation are available.

## CXCL16 and other protein measurements

Membrane, soluble, and total CXCL16 require separate records:

| Analyte record | Required source | Essential metadata |
|---|---|---|
| membrane CXCL16 | viable non-permeabilized cells or quantitative tissue imaging | antibody clone/lot, staining temperature, acquisition settings, live-cell gate, compensation/unmixing, reference control |
| soluble CXCL16 | clarified supernatant, plasma, or tissue fluid | assay kit/lot, standard curve, lower/upper quantification limits, dilution, recovery, collection volume, processing time |
| total CXCL16 | lysate or permeabilized cells | extraction method, total protein or cell normalization, antibody/peptide definition |

Soluble concentration and membrane fluorescence are not commensurate variables and must not be directly compared or ratioed without a validated common scale. The default analysis stratifies by molecular form, response unit, and assay. A direct form comparison requires one assay, one common calibrated unit, and an explicit comparability declaration in the analysis configuration.

For conditioned medium, record producer-cell count, viable-cell fraction, medium volume, conditioning duration, centrifugation/filtration, freeze-thaw cycles, and storage temperature. Cell lysis can artifactually increase measured soluble protein; lactate dehydrogenase release or another prespecified lysis control is recommended.

If ADAM10 perturbation is used, record reagent or genetic construct, concentration, exposure, on-target validation, viability, and total CXCL16. A change in soluble-to-surface CXCL16 alone is not sufficient to establish cleavage.

## Chemotaxis, retention, and egress data

### Chemotaxis

Required raw data are event counts or cell tracks, not only endpoint means. Metadata include device or insert type, pore size, membrane coating, cell input, ligand identity and concentration, compartment volumes, incubation time, temperature, flow or static condition, receptor expression, viability, and recovery-counting method. Conditions must distinguish:

- no ligand;
- forward gradient;
- uniform ligand;
- reverse gradient;
- receptor blockade, receptor-negative, or genetic control;
- positive migration control when available.

### Retention

Record the loading phase, ligand and matrix presentation, washout or perfusion profile, shear stress when applicable, imaging interval, segmentation/tracking version, input and recovered viable cells, and the prespecified definition of a retained cell. Raw time series are required for residence-time analyses.

### Egress

Record the source and receiving compartments, lymphatic or endothelial barrier identity, barrier-integrity measurement, gradient orientation, collection schedule, viable input and output counts, imaging tracks when available, and losses to nonspecific surfaces. A reduced recovered count is not labelled retention unless death and technical loss are quantified.

### Replicate identifiers

Every record distinguishes:

- `biological_replicate`: independent patient, donor, or independently established biological model;
- `technical_replicate`: repeat well, field, device channel, staining, or library from the same biological material;
- `run_id`: acquisition or experimental run;
- `operator_id`: pseudonymous operator code when operator effects are relevant.

In the canonical functional table, `product` is `identical_control` or `receptor_matched`. A different module class requires a versioned schema extension rather than recoding it as receptor matching.

## Construct and manufacturing metadata

The modified and identical-control products require versioned construct maps and a structured difference record. At minimum document:

- antigen-binding domain and target;
- hinge, transmembrane, costimulatory, and activation domains;
- promoter and vector backbone;
- chemokine-module sequence and linkage strategy;
- selectable or reporter elements;
- release criteria and assay methods;
- transduction efficiency, CAR surface abundance, vector copy number, viability, identity, and sterility status where applicable;
- activation reagent, cytokines, basal medium, supplements, culture vessel, culture duration, harvest criteria, and cryopreservation;
- deviations from the matched manufacturing procedure.

## Accession and version rules

1. Use stable accessions at the most specific resolvable level: study, sample, experiment, and file/object.
2. Record the repository and retrieval date, but do not use retrieval date as a substitute for a version.
3. When a repository silently replaces an object, treat the new checksum as a new input version and preserve the earlier provenance record.
4. Reference genomes use assembly plus patch release and annotation provider plus release, for example `GRCh38.p14` with a specific GENCODE release.
5. Antibody and imaging panels receive explicit semantic versions; any clone, barcode, or channel change increments the version.
6. Construct, annotation, and analysis configurations are immutable once used for a tagged release.
7. Controlled data accessions may be listed publicly only when permitted by the relevant data-use agreement.

## Checksum procedure

Generate checksums without modifying the raw files:

```bash
sha256sum path/to/input.fastq.gz > path/to/input.fastq.gz.sha256
sha256sum --check path/to/input.fastq.gz.sha256
stat --format='%s' path/to/input.fastq.gz
```

For a collection of files:

```bash
find data/raw -type f -print0 \
  | sort --zero-terminated \
  | xargs -0 sha256sum > checksums.raw.sha256
sha256sum --check checksums.raw.sha256
```

The checksum manifest itself is version-controlled. Raw files are mounted read-only during analysis.

## Privacy and controlled access

No names, medical-record numbers, dates of birth, contact information, free-text clinical notes, facial images, or unshifted identifying dates are accepted. Patient identifiers must be study-specific pseudonyms. The re-identification key remains outside this repository and outside the analysis environment.

Genomic and clinical data are handled under the original consent, institutional approval, repository terms, and data-use agreement. Publication of code does not authorize redistribution of source data.

## Validation before analysis

Run schema and checksum validation before any scientific stage:

```bash
chemokine-cart validate \
  --config config/config.yaml \
  --manifest config/manifest.tsv \
  --check-files
```

Validation failure is a stop condition. Missing raw data, unresolved accessions, checksum mismatches, invalid pairs, or unknown assay versions are not bypassed by creating placeholder values.
