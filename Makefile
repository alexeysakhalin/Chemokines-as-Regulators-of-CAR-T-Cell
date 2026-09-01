SHELL := /bin/bash

PYTHON ?= python3.11
SNAKEMAKE ?= $(PYTHON) -m snakemake
CONFIG ?= config/config.yaml
MANIFEST ?= data/demo/manifest.tsv
OUTPUT_DIR ?= results
CORES ?= 1
SEED ?= 20260901
PUBLIC_MANIFEST ?= config/public_datasets.tsv
DATASET ?=

export PYTHONHASHSEED := $(SEED)
export LANG := C.UTF-8
export LC_ALL := C.UTF-8
export TZ := UTC
export SOURCE_DATE_EPOCH := 1788220800
export OMP_NUM_THREADS := 1
export OPENBLAS_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export NUMEXPR_NUM_THREADS := 1

.PHONY: help environment install demo data-list data-fetch data-verify validate dry-run workflow determinism lint format-check test check

help:
	@echo "environment   Create the pinned Conda environment"
	@echo "install       Install the package and development dependencies"
	@echo "demo          Create a deterministic demonstration dataset"
	@echo "data-list     List checksum-pinned public input files"
	@echo "data-fetch    Download checksum-pinned public inputs into data/raw"
	@echo "data-verify   Verify locally downloaded public inputs"
	@echo "validate      Validate configuration, manifest, and input files"
	@echo "dry-run       Show the planned Snakemake jobs"
	@echo "workflow      Run the complete analysis workflow"
	@echo "determinism   Verify byte-identical scientific outputs"
	@echo "lint          Run static checks"
	@echo "format-check  Verify source formatting"
	@echo "test          Run the test suite"
	@echo "check         Run lint, formatting, tests, and determinism verification"

environment:
	conda env create --file environment.yml

install:
	$(PYTHON) -m pip install --require-hashes --requirement requirements.lock.txt
	$(PYTHON) -m pip install --no-deps --no-build-isolation --editable .

demo:
	$(PYTHON) -m chemokine_cart.cli make-demo --output-dir data/demo --seed $(SEED) --force

data-list:
	$(PYTHON) scripts/fetch_public_data.py --manifest "$(PUBLIC_MANIFEST)" $(if $(DATASET),--dataset "$(DATASET)",) --list

data-fetch:
	$(PYTHON) scripts/fetch_public_data.py --manifest "$(PUBLIC_MANIFEST)" $(if $(DATASET),--dataset "$(DATASET)",)

data-verify:
	$(PYTHON) scripts/fetch_public_data.py --manifest "$(PUBLIC_MANIFEST)" $(if $(DATASET),--dataset "$(DATASET)",) --verify-only

validate:
	$(PYTHON) -m chemokine_cart.cli validate --config "$(CONFIG)" --manifest "$(MANIFEST)" --check-files

dry-run:
	$(SNAKEMAKE) --snakefile Snakefile --cores $(CORES) --dry-run \
		--configfile "$(CONFIG)" --config manifest_path="$(MANIFEST)" output_dir="$(OUTPUT_DIR)" random_seed=$(SEED)

workflow:
	$(SNAKEMAKE) --snakefile Snakefile --cores $(CORES) --rerun-incomplete --printshellcmds \
		--configfile "$(CONFIG)" --config manifest_path="$(MANIFEST)" output_dir="$(OUTPUT_DIR)" random_seed=$(SEED)

determinism: demo
	$(PYTHON) scripts/verify_determinism.py --random-seed $(SEED)

lint:
	$(PYTHON) -m ruff check .

format-check:
	$(PYTHON) -m ruff format --check .

test:
	$(PYTHON) -m pytest

check: lint format-check test determinism
