---
title: "MRMU Example: Urate to CAD_UKB With LDSC Parameters and Local PLINK Clumping"
author: "Generated with MRMU R package"
date: "2026-04-30"
---

# Summary

This report reruns the urate-CAD example using both newly supplied production
resources:

- 1KG EUR PLINK bfile for local LD clumping
- EUR LDSC scores for estimating `Omega`, `C`, and using real SNP-level `L2`

Analysis setup:

- Exposure: `Biomarker_Urate`
- Outcome: `CAD_UKB`
- Confounders: `Metabolic_SBP`, `Metabolic_DBP`
- IV selection threshold: `5e-5`
- Clumping: local PLINK 1.9, `clump_kb = 1000`, `clump_r2 = 0.001`
- LDSC: pairwise estimation across `Biomarker_Urate`, `Metabolic_SBP`,
  `Metabolic_DBP`, and `CAD_UKB`

| exposure | outcome | confounders | beta | SE | p-value | clumped IVs | effective IVs | pi0 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Biomarker_Urate | CAD_UKB | Metabolic_SBP, Metabolic_DBP | 0.012328 | 0.015260 | 0.41919 | 696 | 322.231 | 0.462975 |

# Data

The analysis used the organised files under:

```text
path/to/MRMU_R_github
```

| folder | role |
|---|---|
| `raw_gwas/` | GWAS summary statistics for urate, SBP, DBP, and CAD_UKB |
| `plink_reference/` | `all_1000G_EUR_Phase3.bed/.bim/.fam` |
| `ldscore/eur_w_ld_chr/` | LDSC score files `1-22.l2.ldscore.gz` and `1-22.l2.M_5_50` |
| `tools/plink_mac/` | macOS PLINK 1.9 executable |

After harmonising alleles across the four traits with real LDSC `L2`, the
merged data contained 1,012,146 SNPs.

# Local Clumping

PLINK log summary:

| item | value |
|---|---:|
| candidate urate IVs at p <= 5e-5 | 14,090 |
| LD clumps formed | 696 |
| removed due to LD/reference absence | 13,394 |

The difference from the previous clumped fallback run is expected: this run uses
the real LDSC-merged SNP set, while the earlier fallback used `L2 = 1`.

# LDSC Parameters

Pairwise LDSC results:

| trait1 | trait2 | rg | omega | C |
|---|---|---:|---:|---:|
| Biomarker_Urate | Biomarker_Urate | 1.0000 | 1.163e-07 | 1.1123 |
| Biomarker_Urate | Metabolic_SBP | 0.1553 | 1.634e-08 | 0.0579 |
| Biomarker_Urate | Metabolic_DBP | 0.1316 | 1.331e-08 | 0.0666 |
| Biomarker_Urate | CAD_UKB | 0.2735 | 2.238e-08 | 0.0529 |
| Metabolic_SBP | Metabolic_SBP | 1.0000 | 9.822e-08 | 1.2026 |
| Metabolic_SBP | Metabolic_DBP | 0.7748 | 7.293e-08 | 0.8589 |
| Metabolic_SBP | CAD_UKB | 0.3625 | 2.747e-08 | 0.0460 |
| Metabolic_DBP | Metabolic_DBP | 1.0000 | 9.110e-08 | 1.1833 |
| Metabolic_DBP | CAD_UKB | 0.2798 | 2.044e-08 | 0.0163 |
| CAD_UKB | CAD_UKB | 1.0000 | 5.970e-08 | 1.0296 |

# MR-MU Result

| quantity | value |
|---|---:|
| beta | 0.01232753 |
| standard error | 0.0152602 |
| p-value | 0.41919 |
| clumped IVs | 696 |
| effective IVs | 322.2309 |
| pi0 | 0.4629754 |
| tau.sq | 2.368078e-07 |
| threshold | 5e-5 |
| method label | MR-MU(LDSC Omega/C, local PLINK clump) |

# Interpretation

After adding real LDSC parameters (`Omega`, `C`) and real SNP-level `L2`, the
estimated urate effect on CAD_UKB is small and not statistically significant in
this MR-MU run.

This LDSC-based result supersedes the earlier fallback reports because it uses
the LD score data and the supplied local PLINK reference.

# Reproducibility

The example was run from the package root:

```bash
MRMU_EXAMPLE_DIR=/path/to/MRMU_R_github Rscript scripts/run_urate_cad_mrmu_ldsc_example.R
```

Key output files:

| file | description |
|---|---|
| `results/urate_cad_mrmu_ldsc_clumped/urate_cad_mrmu_ldsc_result.csv` | LDSC-based MR-MU result table |
| `results/urate_cad_mrmu_ldsc_clumped/urate_cad_mrmu_ldsc_fit.rds` | fitted `MRMU()` object |
| `results/urate_cad_mrmu_ldsc_clumped/urate_cad_mrmu_ldsc_input.rds` | clumped MR-MU input object |
| `results/urate_cad_mrmu_ldsc_clumped/MRMU_ldsc_paras_urate_sbp_dbp_cadukb.rds` | estimated `Omega`, `C`, and `Rg` |
| `results/urate_cad_mrmu_ldsc_clumped/ldsc_pairwise_results.csv` | pairwise LDSC table |
| `results/urate_cad_mrmu_ldsc_clumped/run.log` | full run log |
