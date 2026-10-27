# Chemokines as Regulators of CAR-T Cell Therapy

Reproducible analysis framework for testing how chemokine signals influence CAR-T-cell state, trafficking, spatial positioning, retention, and egress. The repository is designed for five linked workstreams:

1. single-cell RNA-seq/CITE-seq of the infusion product and serial blood samples;
2. spatial mapping of chemokine-producing cells relative to vessels, stroma, tumor nests, and CAR-T cells;
3. protein-level validation of chemokine gradients, with explicit separation of soluble and membrane-associated CXCL16;
4. functional chemotaxis, retention, and egress assays;
5. paired comparison of a chemokine-modified CAR-T product with an otherwise identical CAR-T control.

## Publication analysis

Version 1.0.0 contains an exploratory public-data reanalysis of GSE125881, GSE269379, GSE197268, GSE235760, GSE162975, and GSE273170.

## Repository layout

```text
.
├── analyses/public_reanalysis/ # public-cohort scientific analysis and frozen results
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

## Citation and license

Citation metadata are provided in [CITATION.cff](CITATION.cff). Code is released under the [MIT License](LICENSE). Data remain subject to the terms of their originating studies, repositories, consent documents, and data-use agreements.
