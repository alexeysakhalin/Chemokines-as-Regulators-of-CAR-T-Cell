# Third-party data and attribution

The [MIT license](LICENSE) covers the original repository software and documentation. It does not relicense third-party data, source annotations, or the histology background included in a derived visualization. Public availability is not a blanket grant of unrestricted rights; NCBI states this explicitly in its [website and data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).

## Source studies

| Dataset | Original study | Material used in this analysis |
|---|---|---|
| [GSE125881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125881) | Sheih et al., *Nature Communications* (2020), [10.1038/s41467-019-13880-1](https://doi.org/10.1038/s41467-019-13880-1) | Deposited processed expression counts and author annotations |
| [GSE197268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197268) | Haradhvala et al., *Nature Medicine* (2022), [10.1038/s41591-022-01959-0](https://doi.org/10.1038/s41591-022-01959-0) | Public processed 10x counts and annotations linked by the [author code repository](https://github.com/getzlab/Haradhvala_et_al_2022) |
| [GSE235760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235760) | Guerrero-Murillo et al., *Cell Reports Medicine* (2024), [10.1016/j.xcrm.2024.101803](https://doi.org/10.1016/j.xcrm.2024.101803); correction [10.1016/j.xcrm.2025.102026](https://doi.org/10.1016/j.xcrm.2025.102026) | Public CELLxGENE H5AD with author annotations and integer counts |
| [GSE162975](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162975) | Li et al., *Cell Reports* (2023), [10.1016/j.celrep.2023.113263](https://doi.org/10.1016/j.celrep.2023.113263) | Deposited processed UMI matrix and GEO sample metadata |
| [GSE273170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273170) | Maurer et al., *Blood* (2024), [10.1182/blood.2024024381](https://doi.org/10.1182/blood.2024024381) | Public processed RNA matrices and sample metadata |
| [GSE269379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269379) | Lu et al., *Genome Medicine* (2025), [10.1186/s13073-025-01498-6](https://doi.org/10.1186/s13073-025-01498-6) | Author-curated CSF metadata, Visium HD counts and coordinates, and the deposited registered H&E image |

The GSE235760 source-study correction concerns Figure 2L: the grouping of patients above or below median event-free survival and CAR-T persistence was corrected, and a TIM3-positive CD8 T-cell comparison changed from P=0.0001 to P=0.12. The present CXCR6 IP-to-peak analysis uses deposited integer counts, donor/sample identity, timepoint, and CAR-transduction labels; it does not use that median-based outcome grouping or the corrected Figure 2L statistic. See the [published correction](https://pmc.ncbi.nlm.nih.gov/articles/PMC11970379/).

The accession identifies the source study; file URLs, sizes, and digests in the analysis manifests identify the exact inputs used. The software citation does not replace these source citations. No controlled-access raw reads were used or redistributed.

## Derived outputs

Committed scientific tables contain patient- or cohort-level numerical summaries and source-defined coded identifiers needed to reconstruct inclusion and pairing. They are derived analyses of the above studies, not newly acquired clinical data. The complete third-party count matrices, cell-level metadata, and source tissue images are retrieved locally by the download scripts and are excluded from Git.

GSE235760 was retrieved through CELLxGENE Discover, whose [public-data contribution terms](https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data) specify [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The original contributors and linked publication are credited above. The remaining source records retain their originating repository and contributor terms; this repository does not assert that their input files are MIT-licensed.

## Spatial-image attribution

The histology background in `analyses/public_reanalysis/results/frozen/figure6_public_reanalysis.*`, panel F, is adapted from the deposited GSE269379 sample GSM8968967 associated with Lu et al. (2025). That study is published under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as stated in its [rights and permissions notice](https://link.springer.com/article/10.1186/s13073-025-01498-6#rightslink). The source image was cropped and overlaid with transcript-positive bin coordinates; it was not rotated, smoothed, imputed, or segmented. Construction and source-file attribution are documented in [PANEL_F_PROVENANCE.md](analyses/public_reanalysis/docs/PANEL_F_PROVENANCE.md). The figure filename is retained for reproducibility; it does not designate a Figure 6 in the manuscript.

## Synthetic test inputs

Files under `data/demo/` were constructed for software tests and are explicitly identified by `SYNTHETIC_DATA_NOTICE.json`. They are not patient data and must not be represented as biological evidence.
