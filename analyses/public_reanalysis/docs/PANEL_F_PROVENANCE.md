# Provenance and construction of repository-only spatial panel F

The panel is a direct overlay on a public histology image, not a synthetic or generative image.

## Source

- Study: Lu IN, Müller-Miny L, Krekeler C, et al. *Genome Medicine*. 2025;17:71. doi:10.1186/s13073-025-01498-6.
- GEO series: GSE269379.
- Sample: GSM8968967, `ICANS6, VisiumHD`, a fatal grade-4 ICANS case.
- Histology member: `GSM8968967_tissue_hires_image.png.gz`; decompressed PNG dimensions 3200 × 3000 RGB; SHA-256 `b9238bf89e13b05cff6085e61326b8c3b33c12543b5f5d7fd9170363a1029b58`.
- Coordinates: `GSM8968967_tissue_positions.parquet.gz`; deposited `tissue_hires_scalef = 1.0`.
- Counts: 10x filtered feature-barcode matrix with 18,085 features, 302,297 in-tissue 8-µm bins, and 15,905,676 nonzero entries.

## Deterministic construction

The deposited FFPE hematoxylin-and-eosin image is displayed without rotation, re-registration, segmentation, deconvolution, smoothing, or imputation. A bin is called positive when its integer feature count is greater than zero. The plotted crop is image x=400–1900 and y=760–2280; all matrix bins whose deposited image coordinates fall inside this crop are eligible.

- CXCL16: 960 positive bins, 976 total UMI, shown as orange points.
- CXCR6: 14 positive bins, 14 total UMI, shown as teal triangles.

Symbol sizes are enlarged for visibility and do not represent the physical 8-µm footprint. The background is the public registered H&E image. The workflow does not generate tissue texture or add unobserved expression.

## Interpretation boundary

The markers denote transcript-positive spatial bins, not individually segmented cells and not proven CAR-T identities. A single section cannot establish a protein gradient, distinguish soluble from membrane-bound CXCL16, or demonstrate migration, retention, egress, efficacy, neurotoxicity, or causality. The panel remains repository-only and is not required for the manuscript's one-paragraph bioinformatic summary.

## Attribution and reuse

The histology background is adapted from Lu et al. (2025), GSE269379/GSM8968967. The associated publication is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); its [rights notice](https://link.springer.com/article/10.1186/s13073-025-01498-6#rightslink) requires attribution and identification of modifications. The crop and coordinate overlay are the modifications described above. The original source image is downloaded locally and is not redistributed as a standalone file. See [third-party notices](../../../THIRD_PARTY_NOTICES.md).
