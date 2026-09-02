PYTHON ?= python
SEED ?= 20260902
PERMUTATIONS ?= 10000

.PHONY: list fetch verify-inputs analysis verify test clean-results

list:
	$(PYTHON) scripts/fetch_public_inputs.py --list

fetch:
	$(PYTHON) scripts/fetch_public_inputs.py

verify-inputs:
	$(PYTHON) scripts/fetch_public_inputs.py --verify-only

analysis:
	PYTHONHASHSEED=$(SEED) MPLBACKEND=Agg $(PYTHON) scripts/run_public_reanalysis.py \
		--gse125881-matrix data/raw/GSE125881/GSE125881_raw.expMatrix.csv.gz \
		--gse125881-metadata data/raw/GSE125881/meta-data.csv \
		--gse125881-crosswalk config/gse125881_sample_crosswalk.tsv \
		--gse125881-soft data/raw/GSE125881/GSE125881_family.soft.gz \
		--metadata data/raw/GSE269379/GSE269379_metadata.csv.gz \
		--raw-dir data/raw/GSE269379/extracted \
		--output-dir results/local \
		--seed $(SEED) \
		--spatial-permutations $(PERMUTATIONS)

verify:
	$(PYTHON) scripts/verify_results.py --observed results/local --expected results/frozen

test:
	$(PYTHON) -m pytest -q

clean-results:
	$(PYTHON) scripts/clean_results.py --path results/local
