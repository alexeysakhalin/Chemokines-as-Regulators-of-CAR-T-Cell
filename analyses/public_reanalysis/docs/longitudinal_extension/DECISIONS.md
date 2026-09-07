# Longitudinal extension decision log

## 2026-09-02: independent unit

The patient, not the cell, barcode, library, or plate, is the independent unit. Technical records from one patient-stage are collapsed before any contrast.

## 2026-09-02: endpoint and tests

The primary endpoint is the within-patient change in the fraction of eligible CAR-T cells with at least one CXCR6 UMI. Each cohort receives a two-sided exact sign test. There is no cross-study pooled P value and no cell-level inferential test.

## 2026-09-07: expansion with GSE235760

GSE235760 was added after direct audit of the checksum-identified CELLxGENE H5AD. All five varni-cel patients have matched manufactured-product and patient-specific expansion-peak samples and pass the locked 25-cell denominator. Because this choice followed feasibility review, it is part of a finalized exploratory family and is not described as prospectively preregistered.

The early Holm family consequently contains three cohort-specific tests: GSE197268, GSE235760, and GSE162975.

## 2026-09-07: GSE273170 role

GSE273170 remains supportive. The `CAR >= 1 UMI` definition follows the source study, but only four of ten public day-7/day-14 pairs pass the locked minimum of 10 CAR-transcript-positive cells at both stages. Public matrices precede part of downstream author QC.

## 2026-09-07: GSE290722 exclusion

GSE290722 is not used for the declared longitudinal CAR-T endpoint. Public repeated-timepoint metadata do not identify CAR-T cells, and the relevant publisher source-data sheet contains only 128 week-4 cell records from 13 patients. It therefore contributes zero eligible patient-matched CAR-T pairs under the locked definition. The source workbook, checksum, worksheet, and patient-level record counts are documented in [`GSE290722_SCREENING.md`](GSE290722_SCREENING.md).

## 2026-09-07: marker-pair sampling frame

Inferential marker-set contrasts reuse the complete CXCR6 contrast-inclusion table. Consequently, a patient excluded for a missing stage or an insufficient eligible-cell denominator cannot enter a marker-set sign test. GSE162975 stage aliases are generated from the committed priority table for both endpoint families before pairing.

## Search boundary

The workflow is a documented public-data secondary analysis, not a claim that every potentially relevant public cohort has been exhaustively harmonized. Further candidates require separate audit of CAR identity, time points, treatment course, and patient pairing before inclusion.
