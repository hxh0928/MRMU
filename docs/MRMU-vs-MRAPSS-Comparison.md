---
title: "Urate to CAD_UKB: MRMU vs MRAPSS Comparison"
author: "Generated with MRMU project workflow"
date: "2026-04-30"
---

# Purpose

This report adds a single-exposure MRAPSS result as a comparison to the MRMU
analysis.

Both analyses use:

- Exposure: `Biomarker_Urate`
- Outcome: `CAD_UKB`
- IV threshold: `5e-5`
- EUR LDSC scores
- local 1000 Genomes EUR PLINK clumping

The difference is:

| method | confounders included |
|---|---|
| MRMU | `Metabolic_SBP`, `Metabolic_DBP` |
| MRAPSS | none |

# Results

| method | beta | SE | p-value | IVs | effective IVs | pi0 |
|---|---:|---:|---:|---:|---:|---:|
| MRMU(SBP+DBP confounders) | 0.01232753 | 0.01526020 | 0.41919 | 696 | 322.2309 | 0.4629754 |
| MRAPSS(no confounders) | 0.04490893 | 0.01613808 | 5.3893e-03 | 696 | 324.7741 | 0.4666295 |

# Interpretation

The single-exposure MRAPSS result, which does not adjust for SBP and DBP,
estimates a larger urate effect on CAD_UKB and is statistically significant.

After adding SBP and DBP as confounders in MRMU, the estimated urate effect is
smaller and no longer statistically significant.

This comparison suggests that part of the apparent urate-CAD_UKB association in
the single-exposure analysis may be explained by shared genetic components with
blood pressure traits.

# Reproducibility

The MRAPSS comparison was run with:

```bash
MRMU_EXAMPLE_DIR=/path/to/MRMU_R_github Rscript scripts/run_urate_cad_mrapss_comparison.R
```

Key output files:

| file | description |
|---|---|
| `results/urate_cad_mrapss_comparison/urate_cad_mrapss_result.csv` | MRAPSS-only result |
| `results/urate_cad_mrapss_comparison/urate_cad_mrmu_vs_mrapss_comparison.csv` | side-by-side comparison |
| `results/urate_cad_mrapss_comparison/urate_cad_mrapss_fit.rds` | fitted MRAPSS object |
| `results/urate_cad_mrapss_comparison/urate_cad_mrapss_input.rds` | clumped MRAPSS input data |
| `results/urate_cad_mrapss_comparison/run.log` | LDSC, clumping, and model fitting log |
