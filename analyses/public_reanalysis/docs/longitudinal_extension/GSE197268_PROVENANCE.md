# GSE197268 verified source and reconstruction specification

## Provenance

- GEO series: `GSE197268`.
- Primary report: Haradhvala NJ, Leick MB, Maurer K, et al. *Nature Medicine* (2022) 28:1848-1859. DOI: `10.1038/s41591-022-01959-0`.
- Author code repository: `getzlab/Haradhvala_et_al_2022`; the README blob inspected for the data links had Git SHA `2f1be445a24cb5e19c9cd07366aa86b7fbce1b6b`.
- Author shared folder: `1vw7J8HqUX22ICZmJ0UjAYEBpVjRJ9U9-`; cell-metadata subfolder: `1Vnmq8L-I9PF6-7ESftbbjLQFmIAad--F`.
- Author metadata file `Fig1_CART_global_obs.csv`: Google Drive file ID `1nmWoNKaJF8Rcj1NL_n8XPjL57hylFjEY`, 59,018,996 bytes, SHA-256 `9671ec215a2897a4bbc5283500f334a90f192d720aa007108452ca6bee249c83`.
- Author subtype file `Fig2_ALLT_CD_classifications.txt`: Google Drive file ID `1qcG_CQ86jiVg_QA6bpeHY_Dk4c5MlnUy`, 6,560,656 bytes, SHA-256 `9230a1c99ab1bfa13a1f77d598ba6ab77e6ac400ce7c03c1450d8ccaf4a95f88`.
- HCA/Azul project UUID: `c0d82ef2-1504-4ef0-9e5e-d8a13e45fdec`; source UUID: `4d97368a-d5ed-436e-b7a1-f2243cccb08b`. The locally inspected project response was 88,397 bytes with SHA-256 `7bd768f4646e643895103fec2e2af846bc57fd0d67041c4c6264feedbd02873a` and identified 109 contributor tar archives sourced from GEO. The 18 selected archives matched its recorded names, sizes, and SHA-256 digests exactly.

## Locked primary reconstruction

The sampling frame contains 21 patients with author-QC metadata at both the first-course infusion product (`IP`) and the physically enriched day-7 CAR-T fraction (`D7-CART`). Before reading CXCR6 counts, the author annotations were joined and the number of product-matched cells with `CAR == True` and subtype `CD8 T` was counted at each stage. Nine patients passed the locked minimum of 25 cells at both stages; 12 failed. The exact all-patient screen is committed as `results/longitudinal_extension/frozen/audit/gse197268_screening_exclusions.tsv`. Only the 18 GEO bundles belonging to the nine eligible patients are downloaded, which prevents the raw reconstruction from becoming a convenience sample selected on the CXCR6 result.

For those nine patients, a raw cell is retained only when its deposited 10x barcode is present in the author `Fig1_CART_global_obs.csv` index, its patient, stage, and product agree with the sample crosswalk, `CAR == True`, and the joined author classification in `Fig2_ALLT_CD_classifications.txt` is exactly `CD8 T`. The endpoint is the fraction of these eligible cells with an integer `CXCR6` UMI count greater than zero. Each patient contributes one paired difference, calculated as `D7-CART fraction - IP fraction`.

Directional consistency is evaluated with an exact two-sided sign test after removing zero differences. Cells are not treated as independent observations. The primary result is `n = 9`, seven positive and two negative differences, median change `0.23626373626373628` (23.6264 percentage points), exact sign-test `P = 0.1796875`, and Holm-adjusted `P = 0.5390625` in the three-cohort early family.

The product strata are descriptive only: axi-cel `n = 5`, five increases, median `0.26120981387478853`, `P = 0.0625`; tisa-cel `n = 4`, two increases and two decreases, median `-0.007292857971591002`, `P = 1`. These strata do not support a time-by-product interaction test.

## Audit files and exclusions

- `config/longitudinal_sources.tsv` contains the two author annotation sources and the 18 exact GEO sample URLs, sizes, and SHA-256 values.
- `config/gse197268_sample_crosswalk.tsv` records one row per included GEO sample, including patient, product, stage, raw matrix cells, author-QC matches, eligible cells, CXCR6-positive cells, and source hash.
- `results/longitudinal_extension/frozen/audit/gse197268_screening_exclusions.tsv` contains all 21 patients in the outcome-blind denominator screen.
- `results/longitudinal_extension/frozen/audit/gse197268_verified_patient_aggregates.tsv` contains the exact nine patient rows used in the contrast.
- `results/longitudinal_extension/frozen/audit/gse197268_product_stratified_results.tsv` contains the two descriptive product strata.

Six local `.part` files are incomplete download prefixes. They fail gzip integrity testing and must never appear in the source manifest, analysis inputs, frozen outputs, or checksum file. The source data themselves are not redistributed; only exact source identifiers, URLs, sizes, hashes, crosswalks, and derived patient-level summaries should be committed.
