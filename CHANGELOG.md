# Release history

## 1.0.0 — 2026-09-07

Publication version of the exploratory CAR-T chemokine reanalysis.

- Patient-level longitudinal CXCR6 transcript-detection analysis of GSE197268, GSE235760, GSE162975, and GSE273170, with cohort-specific CAR definitions, denominator thresholds, exclusion records, and paired comparisons.
- Exact two-sided sign tests, Holm adjustment across the three early comparisons, and explicitly exploratory later and product-stratified summaries.
- Descriptive GSE125881 longitudinal analysis and GSE269379 ICANS-CSF and single-section spatial analysis retained with source provenance.
- Download manifests, pinned Python environments, deterministic workflows, frozen derived tables, and output-verification commands for Linux reproduction.
- Reusable root-package functions and a separately labeled synthetic demonstration for software testing.
- Updated publication documentation, software citation, and third-party attribution.

The analysis does not support a universal early increase in CXCR6-detected CAR-T-cell fraction across the three early cohorts. The observed temporal and between-cohort variation does not establish a chemokine protein gradient, trafficking, retention, egress, treatment efficacy, or causality. The spatial overlay remains a descriptive repository visualization, not an additional manuscript figure.

Third-party source matrices are downloaded from their originating repositories rather than redistributed. See [third-party notices](THIRD_PARTY_NOTICES.md) and the [public-analysis README](analyses/public_reanalysis/README.md).
