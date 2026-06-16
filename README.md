# MRMU

MRMU is an R package for Mendelian randomization analyses with multiple
exposures/confounders using GWAS summary statistics.

## Installation

```r
# install.packages("remotes")
remotes::install_github("hxh0928/MRMU")
```

The package depends on `MRAPSS`, `Matrix`, `expm`, and `data.table`.
`MRAPSS` is available from the MRC IEU r-universe repository.

## Demo: Urate to CAD With SBP and DBP Confounders

```r
library(MRMU)

# Load packaged example data:
# - urate_cad_mrmu_iv: clumped IV data
# - urate_cad_ldsc_parameters: LDSC-estimated C, Omega, and Rg matrices
data("urate_cad_mrmu_iv")
data("urate_cad_ldsc_parameters")

fit <- MRMU(
  MRdat = urate_cad_mrmu_iv,
  exposure = "Biomarker_Urate",
  confounders = c("Metabolic_SBP", "Metabolic_DBP"),
  outcome = "CAD_UKB",
  C = urate_cad_ldsc_parameters$C,
  Omega = urate_cad_ldsc_parameters$Omega
)

data.frame(
  exposure = fit$exposure,
  outcome = fit$outcome,
  beta = fit$beta1,
  se = fit$beta1.se,
  p_value = fit$beta1.pvalue,
  n_iv = fit$nIV,
  effective_iv = fit$nvalid,
  pi0 = fit$pi0
)
```

Expected result:

| exposure | outcome | beta | se | p-value | IVs | effective IVs | pi0 |
|---|---|---:|---:|---:|---:|---:|---:|
| Biomarker_Urate | CAD_UKB | 0.01232753 | 0.01526020 | 4.1919e-01 | 696 | 322.2309 | 0.4629754 |

The packaged demo data include the final clumped IV input and LDSC parameters,
so this example can run without downloading the full GWAS, LD score,
or PLINK reference files.

If you want to rebuild the example from raw summary statistics, first harmonise
the four traits with `MV_harmonise()`, estimate `C` and `Omega`, and then pass
`bfile` and `plink_bin` to `MRMU_Input()` for
PLINK-based clumping.

## Tutorial and example reports

The urate-CAD tutorial with SBP and DBP as confounders is available in:

- `docs/MRMU-Rpackage-Tutorial.md`
- `tutorials/MRMU_Rpackage_Tutorial.md`
- `tutorials/MRMU_Rpackage_Tutorial.pdf`

The comparison between MRMU and MRAPSS without confounders is available in:

- `docs/MRMU-vs-MRAPSS-Comparison.md`
- `reports/urate_cad_mrmu_vs_mrapss_comparison.md`
- `reports/urate_cad_mrmu_vs_mrapss_comparison.pdf`

The full raw GWAS, LD score, and PLINK reference data are not stored in this
repository. The tutorial includes a placeholder table for shared download links.

## Reproducing the paper analyses

Code, processed inputs, and selected result tables for reproducing the paper
analyses are maintained in the companion repository:

- <https://github.com/hxh0928/MRMU-reproduce>
