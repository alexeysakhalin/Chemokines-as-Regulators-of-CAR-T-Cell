# GSE290722 screening record

## Decision

GSE290722 was screened but was not included in the declared patient-level
longitudinal CXCR6 endpoint. The exclusion was made before examination of a
within-patient CXCR6 contrast: the public materials do not provide the two
CAR-T stages required by the analysis definition.

## Sources examined

- GEO series: GSE290722.
- Source study: Zhang et al., *Nature Communications* (2025),
  DOI `10.1038/s41467-025-59904-x`.
- Publisher source-data workbook:
  `41467_2025_59904_MOESM5_ESM.xlsx`.
- Workbook size: 16,864,047 bytes.
- Workbook SHA-256:
  `c5a4c2e11624dc0541870875e6df5f959ad7eb4f23ea91bbf5f630046e10eef9`.
- Relevant worksheet: `Supplementary Figure1b`.

The workbook is not redistributed here. It is available from the article's
publisher page and remains subject to the publisher's terms.

## Audit result

The relevant worksheet contains 128 non-empty cell records from 13 patients.
Every record in that sheet is labelled `4 weeks after CAR T`; no infusion-product
or alternative post-infusion CAR-T stage is represented. The patient record
counts in the sheet are:

| Patient | Records |
|---|---:|
| P2 | 4 |
| P5 | 1 |
| P6 | 79 |
| P11 | 5 |
| P14 | 5 |
| P15 | 1 |
| P19 | 1 |
| P20 | 1 |
| P22 | 1 |
| P24 | 4 |
| P26 | 12 |
| P27 | 11 |
| P29 | 3 |

The repeated-timepoint metadata in the broader public deposit do not provide a
CAR-T identity field that can be applied consistently across the required
stages. Therefore, this cohort contributes zero eligible patient-matched pairs
under the locked endpoint, which requires recoverable CAR-T identity and both
specified biological stages for the same patient.

## Interpretation

This is a feasibility exclusion, not a negative biological result. It does not
imply that CXCR6 is absent in GSE290722, and no value from this cohort enters a
sign test, Holm adjustment, effect-size summary, or manuscript conclusion.
