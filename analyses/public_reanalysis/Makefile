PYTHON ?= python
SEED ?= 20260902
PERMUTATIONS ?= 10000

.PHONY: list fetch verify-inputs analysis verify fetch-longitudinal longitudinal verify-longitudinal test clean-results clean-longitudinal

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

fetch-longitudinal:
	$(PYTHON) scripts/fetch_longitudinal_inputs.py

longitudinal:
	PYTHONHASHSEED=$(SEED) $(PYTHON) scripts/run_longitudinal_extension.py \
		--gse197268-dir data/raw/GSE197268 \
		--gse197268-obs data/raw/GSE197268/Fig1_CART_global_obs.csv \
		--gse197268-subtype data/raw/GSE197268/Fig2_ALLT_CD_classifications.txt \
		--gse162975-matrix data/raw/GSE162975/GSE162975_CART_BALL_UMI.csv.gz \
		--gse235760-h5ad data/raw/GSE235760/2aad61a0-31b8-41f5-8b8a-dbd971269f7b.h5ad \
		--gse273170-archive data/raw/GSE273170/GSE273170_RAW.tar \
		--output-dir results/longitudinal_extension/local \
		--seed $(SEED)

verify-longitudinal:
	$(PYTHON) scripts/verify_longitudinal_results.py \
		--observed results/longitudinal_extension/local \
		--expected results/longitudinal_extension/frozen

test:
	$(PYTHON) -m pytest -q

clean-results:
	$(PYTHON) scripts/clean_results.py --path results/local

clean-longitudinal:
	$(PYTHON) scripts/clean_longitudinal_results.py \
		--path results/longitudinal_extension/local
