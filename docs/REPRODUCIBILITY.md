# Reproducibility guide

## Reproducibility standard

A figure or table is reproducible only when it can be regenerated from immutable inputs in a clean Linux environment using a recorded command, configuration, code commit, software environment, and random seed. A rendered image without its source table and provenance is not considered reproducible.

The public-data scientific analysis is documented in [`analyses/public_reanalysis/README.md`](../analyses/public_reanalysis/README.md). It has its own environment, input manifests, deterministic commands, and frozen results. The guide below describes the separate root package and its synthetic demonstration. The two environments have different dependency versions and must be created separately.

## Supported platform

- Linux x86_64
- Python 3.11.11
- Mamba or Conda-compatible environment manager
- Snakemake 8.28.0 and the resolved Python dependency graph as pinned in `requirements.lock.txt`
- `environment.yml` as the Python/Conda bootstrap for that lock
- Git with long-path and symbolic-link support

Run all commands from the repository root. Do not execute the workflow from a directory containing another checkout with the same package name.

## Clean installation

Clone the publication version:

```bash
git clone --branch v1.0.0 \
  https://github.com/alexeysakhalin/Chemokines-as-Regulators-of-CAR-T-Cell.git
cd Chemokines-as-Regulators-of-CAR-T-Cell
git status --short
git rev-parse HEAD
```

Create the pinned environment and install the local package without resolving unpinned dependencies:

```bash
mamba env create --file environment.yml
mamba activate chemokine-cart
python --version
python -m pip install --no-deps --no-build-isolation --editable .
chemokine-cart --help
```

Alternatively, install the hash-locked Python runtime and development dependencies:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip==25.0.1
python -m pip install --require-hashes --requirement requirements.lock.txt
python -m pip install --no-deps --no-build-isolation --editable .
chemokine-cart --help
```

The equivalent project wrappers are:

```bash
make environment
mamba activate chemokine-cart
make install
```

If `environment.yml` changes, create a new environment rather than updating an old one in place. Record the explicit package export:

```bash
mamba list --explicit > environment.explicit.txt
sha256sum environment.yml requirements.lock.txt environment.explicit.txt
```

The Conda environment stores the deterministic variables below. For a plain virtual
environment or another direct CLI setup, export the same values before execution:

```bash
export PYTHONHASHSEED=20260901
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export TZ=UTC
export SOURCE_DATE_EPOCH=1788220800
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
```

## Synthetic end-to-end check

Generate deterministic synthetic inputs:

```bash
chemokine-cart make-demo --output-dir data/demo --seed 20260901 --force
```

Validate the configuration, manifest, and file checksums:

```bash
chemokine-cart validate \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --check-files
```

Inspect the planned execution before running:

```bash
chemokine-cart run \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901 \
  --dry-run
```

Run the complete workflow:

```bash
snakemake \
  --snakefile Snakefile \
  --cores 1 \
  --configfile config/config.yaml \
  --config manifest_path=data/demo/manifest.tsv \
           output_dir=results/demo \
           random_seed=20260901 \
  --rerun-incomplete \
  --printshellcmds
```

The equivalent wrapper is:

```bash
make demo
make validate
make dry-run
make workflow
make check
```

The synthetic dataset tests execution and expected file structure. It is not evidence for a biological effect.

Verify that the scientific tables and figure files are identical across two clean executions:

```bash
python scripts/verify_determinism.py --random-seed 20260901
```

This check intentionally excludes runtime provenance and workflow markers, which contain execution-specific timestamps and paths. It compares every file under `tables/` and `figures/`, including the figure manifest, by SHA-256 digest.

## Stage-by-stage execution

Available stages are `validate`, `single-cell`, `spatial`, `functional`, `figures`, `report`, and `all`.
Run the validation stage first; the individual scientific-stage commands below do not replace
schema and checksum validation.

```bash
python -m chemokine_cart.pipeline \
  --stage validate \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901
```

```bash
python -m chemokine_cart.pipeline \
  --stage single-cell \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901

python -m chemokine_cart.pipeline \
  --stage spatial \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901

python -m chemokine_cart.pipeline \
  --stage functional \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901

python -m chemokine_cart.pipeline \
  --stage figures \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901

python -m chemokine_cart.pipeline \
  --stage report \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901
```

For a complete run through the package interface:

```bash
chemokine-cart run \
  --config config/config.yaml \
  --manifest data/demo/manifest.tsv \
  --output-dir results/demo \
  --random-seed 20260901
```

## Running authorized real data

### Prepare an immutable input area

Raw data are never committed. Because the validator rejects absolute and path-escaping inputs, authorized storage must be staged or bind-mounted under the ignored project-relative directory `data/raw/`. A Linux read-only bind mount avoids copying large files:

```bash
export CHEMOKINE_CART_RAW=/absolute/path/to/authorized/read_only_inputs
test -d "${CHEMOKINE_CART_RAW}"
mkdir -p data/raw
sudo mount --bind "${CHEMOKINE_CART_RAW}" data/raw
sudo mount --options remount,bind,ro data/raw

find data/raw -type f -print0 \
  | sort --zero-terminated \
  | xargs -0 sha256sum > checksums.raw.sha256
```

If bind mounts are unavailable, copy only authorized inputs into `data/raw/` and make that directory read-only for the run. Symlinks that resolve outside the project root are intentionally rejected. Populate `config/manifest.tsv` with project-relative paths such as `data/raw/file.tsv`, pseudonymous identifiers, immutable accessions, file sizes, SHA-256 digests, assay versions, and pairing metadata. Confirm that use of every input is permitted by its consent and data-use agreement.

### Validate before computation

```bash
chemokine-cart validate \
  --config config/config.yaml \
  --manifest config/manifest.tsv \
  --check-files
```

Checksum mismatch, unresolved accession, invalid product pairing, absent biological-unit identifier, or unknown reference version is a stop condition. Do not bypass validation by editing the generated outputs.

The default evidence policy also blocks every non-synthetic scientific stage unless all manifest records have `evidence_eligible=true`. Validation itself remains available so that records on QC hold can be inspected. The synthetic-only demonstration is the sole exception and its figures remain explicitly labelled as non-evidentiary software tests.

### Dry run and execute

```bash
chemokine-cart run \
  --config config/config.yaml \
  --manifest config/manifest.tsv \
  --output-dir results \
  --random-seed 20260901 \
  --dry-run

snakemake \
  --snakefile Snakefile \
  --cores 8 \
  --configfile config/config.yaml \
  --config manifest_path=config/manifest.tsv \
           output_dir=results \
           random_seed=20260901 \
  --rerun-incomplete \
  --printshellcmds
```

The seed is explicit for stochastic operations. Changing it creates a distinct run and must be recorded. Do not use `--touch`, manual timestamp changes, or edited intermediate files to force completion.

## Reproducing figures locally

Figures are generated from versioned analysis tables rather than manually transcribed values. Run the upstream analysis stages before the figure stage:

```bash
python -m chemokine_cart.pipeline \
  --stage all \
  --config config/config.yaml \
  --manifest config/manifest.tsv \
  --output-dir results \
  --random-seed 20260901
```

To regenerate figures after upstream tables have been verified:

```bash
python -m chemokine_cart.pipeline \
  --stage figures \
  --config config/config.yaml \
  --manifest config/manifest.tsv \
  --output-dir results \
  --random-seed 20260901
```

Each figure must have:

- a source table in a non-proprietary format;
- a script or workflow rule;
- a configuration record for labels, contrasts, and thresholds;
- vector output (`.svg` or `.pdf`) when supported;
- a raster preview with recorded resolution;
- output SHA-256 checksum;
- a caption stating the biological unit and statistical model.

Patient-level or donor-level quantitative plots must show the independent biological-unit points. Single-cell points, tracks, or regions may be displayed descriptively but may not replace biological replicates.

## Provenance record

A complete run through the `report` stage automatically records:

- exact Git commit and dirty-worktree status;
- configuration and manifest paths, byte sizes, and SHA-256 digests;
- local source-input paths, byte sizes, and SHA-256 digests;
- process arguments, a canonical reproduction command, and the effective random seed;
- operating system, architecture, Python version, and selected package versions;
- `requirements.lock.txt` and `environment.yml` byte sizes and SHA-256 digests;
- generated-output byte sizes and SHA-256 digests.

For an independently reproduced analysis, retain:

- repository URL, exact commit, and release tag;
- run start and finish timestamps and the exact outer Snakemake or Make command;
- explicit Conda export or immutable container digest;
- rule-level logs;
- reference genome, annotation, panel, and construct identifiers;
- excluded samples, regions, cells, tracks, and exclusion reasons.

Recommended capture before a real run:

```bash
mkdir -p results/provenance
git rev-parse HEAD > results/provenance/git_commit.txt
git status --porcelain=v1 > results/provenance/git_status.txt
sha256sum config/config.yaml config/manifest.tsv \
  > results/provenance/config_manifest.sha256
mamba list --explicit > results/provenance/environment.explicit.txt
uname -a > results/provenance/system.txt
```

A non-empty `git_status.txt` indicates that results were produced from uncommitted code. Such a run may be used during development but must not support a release figure.

## Determinism and numerical tolerance

Set random seeds in the top-level configuration and propagate them to Python, NumPy, statistical libraries, bootstraps, spatial permutations, and any parallel worker. Preserve stable sorting by explicit identifiers before seeded operations.

Byte-identical output is expected for deterministic text tables in the same locked environment. Floating-point outputs generated on different hardware may be compared using a documented tolerance. Tolerance-based validation must compare scientific values and row identities, not only file size or plotting appearance.

## Tests

Run all checks before sharing a commit:

```bash
make check
```

At minimum the test suite should cover:

- manifest and configuration validation;
- checksum failure and missing-file failure;
- invalid modified/control pair detection;
- patient-level rather than cell-level aggregation;
- soluble and membrane CXCL16 separation;
- chemotaxis control orientation;
- retention and egress denominator definitions;
- deterministic bootstrap or permutation output under a fixed seed;
- successful synthetic end-to-end workflow;
- figure generation from source tables.

## Branch and release policy

`main` is the publication branch. Cite the immutable `v1.0.0` release or its exact commit when reproducing this version. The retained development branch records the analysis history; it is not the preferred citation target.

The release includes fixed endpoints, audited sample crosswalks, source-file checksums, pinned environments, tests, and frozen output tables. Source data remain in their originating repositories. The source studies and any applicable reuse conditions must be cited separately from this software.

Independent reproduction should use a fresh checkout and the environment for the relevant module, then compare all generated outputs with the frozen tables using the documented verifier. An archival DOI may be added once that exact release has been deposited; no DOI is asserted here.

## Independent reproduction checklist

An independent analyst should be able to answer yes to each item:

- Can the exact source files be obtained under the stated access conditions?
- Do all byte sizes and SHA-256 digests match?
- Is the biological-unit and modified/control pairing unambiguous?
- Are the CAR construct and chemokine module versioned?
- Are state and CAR-origin annotations versioned and auditable?
- Are soluble and membrane CXCL16 represented as distinct analytes?
- Are chemotaxis, retention, and egress controls identifiable in the manifest?
- Does a clean environment pass validation and tests?
- Can every figure be regenerated from a source table by a recorded command?
- Do captions state the biological unit, model, and multiplicity procedure?

If any answer is no, the affected output remains developmental and must not be presented as fully reproducible.
