---
title: "A Real Example for Performing Summary-Data MR Analysis with the MRMU Package"
author: "Xianghong Hu"
date: "2026-04-30"
---

# Introduction

MRMU is an R package for Mendelian randomization analyses with one primary
exposure, one outcome, and multiple potential confounder traits using
genome-wide summary statistics.

This tutorial walks through a complete real-data example:

- Exposure: `Biomarker_Urate`
- Outcome: `CAD_UKB`
- Confounders: `Metabolic_SBP`, `Metabolic_DBP`
- Instrument threshold: `p <= 5e-5`
- LD clumping: local PLINK with the 1000 Genomes EUR reference panel
- LDSC parameters: estimated from EUR LD scores

The full workflow has four practical steps:

1. Install and load packages.
2. Prepare and harmonise summary statistics.
3. Estimate LDSC nuisance parameters and select clumped instruments.
4. Fit the MRMU model and read the output.

# Step 0: Installation and Loading Packages

Install MRMU from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hxh0928/MRMU", dependencies = TRUE)
```

Load the package and a fast table-reading package:

```r
library(MRMU)
library(data.table)
```

## Data Download Links

This tutorial assumes the example data have been downloaded and organised in the
folder structure shown below. The shared download links will be added after they
are available.

| data type | link |
|---|---|
| GWAS summary statistics | To be added |
| 1000 Genomes EUR PLINK reference | To be added |
| EUR LD score files | To be added |
| macOS/Linux PLINK executable | To be added |

For this tutorial, all input files are organised under one project folder:

```r
bundle_dir <- "path/to/MRMU_R_github"

data_dir <- file.path(bundle_dir, "raw_gwas")
ldscore_dir <- file.path(bundle_dir, "ldscore", "eur_w_ld_chr")
bfile_prefix <- file.path(bundle_dir, "plink_reference", "all_1000G_EUR_Phase3")
plink_bin <- file.path(bundle_dir, "tools", "plink_mac", "plink")
```

Replace `path/to/MRMU_R_github` with the folder where you downloaded and
organised the example data. The four paths mean:

| object | meaning |
|---|---|
| `data_dir` | formatted GWAS summary-statistic files |
| `ldscore_dir` | LDSC files `1-22.l2.ldscore.gz` and `1-22.l2.M_5_50` |
| `bfile_prefix` | PLINK reference prefix without `.bed/.bim/.fam` suffix |
| `plink_bin` | local PLINK executable |

# Step 1: Prepare and Harmonise Data

## 1.1. Required Input Columns

Each GWAS file should already be formatted with the following columns:

| column | meaning |
|---|---|
| `SNP` | rsID |
| `A1` | effect allele |
| `A2` | non-effect allele |
| `Z` | signed z-score |
| `N` | sample size |
| `chi2` | chi-square statistic |
| `P` | p-value |

For example, the urate file starts like this:

```text
SNP A1 A2 Z N chi2 P
rs1000000 A G -0.1127863 362754 0.01272074 0.9102
rs10000010 C T -0.5980598 362754 0.35767558 0.5498
rs10000023 T G 2.4577769 362754 6.04066718 0.01398
```

Create a file map for the four traits:

```r
traits <- c("Biomarker_Urate", "Metabolic_SBP", "Metabolic_DBP", "CAD_UKB")

data.files <- data.frame(
  trait.name = traits,
  file.dir = file.path(data_dir, traits),
  stringsAsFactors = FALSE
)
```

## 1.2. Read LD Scores

MRMU uses SNP-level LD scores in the LDSC parameter model. Read and combine all 22
chromosomes:

```r
ld <- rbindlist(lapply(1:22, function(chr) {
  fread(file.path(ldscore_dir, paste0(chr, ".l2.ldscore.gz")))
}), use.names = TRUE)
```

The combined `ld` table must contain at least:

```text
SNP, L2
```

## 1.3. Harmonise All Traits

Use `MV_harmonise()` to align all traits to the same effect allele. This
function:

- merges the GWAS files by `SNP`;
- checks whether alleles are aligned or flipped;
- flips `Z` when needed;
- converts `Z` and `N` to SNP effect estimates and standard errors;
- merges in the LD score column `L2`.

```r
dat <- MV_harmonise(data.files, ld)
```

In this example, the harmonised data contained 1,012,146 SNPs.

The output columns are:

```text
SNP, A1, A2,
b_Biomarker_Urate, se_Biomarker_Urate, pval_Biomarker_Urate,
b_Metabolic_SBP, se_Metabolic_SBP, pval_Metabolic_SBP,
b_Metabolic_DBP, se_Metabolic_DBP, pval_Metabolic_DBP,
b_CAD_UKB, se_CAD_UKB, pval_CAD_UKB,
L2
```

For each trait, the effect estimate and standard error are computed as:

```text
b  = Z / sqrt(N)
se = 1 / sqrt(N)
```

# Step 2: Estimate LDSC Parameters

MRMU needs two nuisance matrices for the exposure, confounders, and outcome:

| matrix | dimension | interpretation |
|---|---:|---|
| `Omega` | 4 by 4 | covariance matrix of polygenic effects |
| `C` | 4 by 4 | covariance/intercept matrix capturing sample structure and correlated errors |

Here the trait order is:

```r
traits_for_model <- c("Biomarker_Urate", "Metabolic_SBP", "Metabolic_DBP", "CAD_UKB")
```

Estimate the LDSC-based parameters:

```r
MRMU_ldsc_paras <- est_MRMU_ldsc_paras(
  data.files = data.files,
  ldscore.files = ldscore_dir
)

rownames(MRMU_ldsc_paras$Omega) <- colnames(MRMU_ldsc_paras$Omega) <- traits
rownames(MRMU_ldsc_paras$C) <- colnames(MRMU_ldsc_paras$C) <- traits
rownames(MRMU_ldsc_paras$gc) <- colnames(MRMU_ldsc_paras$gc) <- traits

diag(MRMU_ldsc_paras$C) <- pmax(diag(MRMU_ldsc_paras$C), 1)
```

The last line prevents diagonal elements of `C` from being smaller than 1,
matching the production workflow.

For this example, the estimated genetic correlations were:

| trait1 | trait2 | genetic correlation |
|---|---|---:|
| Biomarker_Urate | Metabolic_SBP | 0.1553 |
| Biomarker_Urate | Metabolic_DBP | 0.1316 |
| Biomarker_Urate | CAD_UKB | 0.2735 |
| Metabolic_SBP | Metabolic_DBP | 0.7748 |
| Metabolic_SBP | CAD_UKB | 0.3625 |
| Metabolic_DBP | CAD_UKB | 0.2798 |

The estimated `Omega` matrix was:

|  | Biomarker_Urate | Metabolic_SBP | Metabolic_DBP | CAD_UKB |
|---|---:|---:|---:|---:|
| Biomarker_Urate | 1.163e-07 | 1.634e-08 | 1.331e-08 | 2.238e-08 |
| Metabolic_SBP | 1.634e-08 | 9.822e-08 | 7.293e-08 | 2.747e-08 |
| Metabolic_DBP | 1.331e-08 | 7.293e-08 | 9.110e-08 | 2.044e-08 |
| CAD_UKB | 2.238e-08 | 2.747e-08 | 2.044e-08 | 5.970e-08 |

The estimated `C` matrix was:

|  | Biomarker_Urate | Metabolic_SBP | Metabolic_DBP | CAD_UKB |
|---|---:|---:|---:|---:|
| Biomarker_Urate | 1.1120 | 0.0579 | 0.0666 | 0.0529 |
| Metabolic_SBP | 0.0579 | 1.2026 | 0.8589 | 0.0460 |
| Metabolic_DBP | 0.0666 | 0.8589 | 1.1833 | 0.0163 |
| CAD_UKB | 0.0529 | 0.0460 | 0.0163 | 1.0296 |

# Step 3: Select Instruments by Local LD Clumping

Use urate p-values to select candidate instruments and then LD-clump them with
the local 1000 Genomes EUR reference panel.

```r
p_thresh <- 5e-5

MRMUInput <- MRMU_Input(
  dat = dat,
  exposure = "Biomarker_Urate",
  confounders = c("Metabolic_SBP", "Metabolic_DBP"),
  outcome = "CAD_UKB",
  p_thresh = p_thresh,
  bfile = bfile_prefix,
  plink_bin = plink_bin,
  clump_kb = 1000,
  clump_r2 = 0.001
)
```

The clumping parameters mean:

| argument | value | meaning |
|---|---:|---|
| `p_thresh` | `5e-5` | initial IV selection p-value threshold |
| `clump_kb` | `1000` | remove nearby correlated SNPs within 1000 kb |
| `clump_r2` | `0.001` | LD threshold for pruning correlated SNPs |
| `bfile` | local 1KG EUR reference | reference genotypes for LD calculation |

This run used:

```text
PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)
```

The clumping summary was:

| item | value |
|---|---:|
| candidate urate IVs at p <= 5e-5 | 14,090 |
| clumped IVs | 696 |
| removed due to LD/reference absence | 13,394 |

`MRMUInput` is a list. The important fields are:

| field | dimension | meaning |
|---|---:|---|
| `b.exp` | `m by 3` | SNP effects on urate, SBP, and DBP |
| `se.exp` | `m by 3` | standard errors for `b.exp` |
| `pval.exps` | `m by 3` | p-values for urate, SBP, and DBP |
| `b.out` | `m by 1` | SNP effects on CAD_UKB |
| `se.out` | `m by 1` | standard errors for CAD_UKB effects |
| `L2` | `m by 1` | SNP-level LD scores |
| `Threshold` | scalar | adjusted IV threshold used for selection-bias correction |

Here `m = 696` after clumping.

# Step 4: Fit MRMU

Now fit the MRMU model. The order of rows and columns in `Omega` and `C` must
match `c(exposure, confounders, outcome)`.

```r
traits_for_model <- c(
  "Biomarker_Urate",
  "Metabolic_SBP",
  "Metabolic_DBP",
  "CAD_UKB"
)

fit.MRMU <- MRMU(
  MRdat = MRMUInput,
  exposure = "Biomarker_Urate",
  confounders = c("Metabolic_SBP", "Metabolic_DBP"),
  outcome = "CAD_UKB",
  Omega = MRMU_ldsc_paras$Omega[traits_for_model, traits_for_model],
  C = MRMU_ldsc_paras$C[traits_for_model, traits_for_model],
  Cor.SelectionBias = TRUE,
  cut.confounders = FALSE,
  tol = 1e-8,
  a = 1,
  b = 1
)
```

Important arguments:

| argument | meaning |
|---|---|
| `Cor.SelectionBias = TRUE` | use the adjusted `Threshold` to correct IV selection bias |
| `cut.confounders = FALSE` | keep both SBP and DBP in this focused example |
| `tol = 1e-8` | convergence tolerance for the variational EM algorithm |
| `a = 1, b = 1` | prior hyperparameters for confounder effects |

The fitted object includes:

| field | meaning |
|---|---|
| `beta1` | estimated causal effect of urate on CAD_UKB |
| `beta1.se` | standard error |
| `beta1.pvalue` | likelihood-ratio-test p-value |
| `nIV` | number of clumped IVs |
| `nvalid` | effective number of foreground IVs |
| `pi0` | estimated foreground probability |
| `SigmaX` | fitted covariance matrix for exposure/confounder foreground effects |
| `tau.sq` | fitted variance of direct outcome effects |
| `post` | posterior quantities, including per-SNP posterior probabilities |

## Step 4 Result: MRMU

The LDSC-based MRMU result for this example was:

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
| method | MR-MU with LDSC `Omega/C` and local PLINK clumping |

After accounting for SBP and DBP and using LDSC-estimated
parameters, the estimated effect of urate on CAD_UKB was small and not
statistically significant.

# Step 5: Run the No-Confounder MRAPSS Comparison

To understand how much the SBP and DBP confounders change the result, we also
run a no-confounder analysis. This comparison uses the same exposure, outcome,
LD scores, local PLINK reference, and IV threshold, but does not include
`Metabolic_SBP` or `Metabolic_DBP`.

This comparison is useful because it answers a concrete question:

```text
What would the urate-CAD_UKB result look like if we ignored SBP and DBP?
```

## 5.1. Reuse the LDSC Parameters From Step 2

There is no need to repeat LDSC estimation for the no-confounder comparison.
Step 2 already estimated `Omega` and `C` for all four traits. The
no-confounder model only needs the 2 by 2 submatrices for:

```r
exposure <- "Biomarker_Urate"
outcome <- "CAD_UKB"
```

Extract the relevant rows and columns from the Step 2 output:

```r
C.no.conf <- MRMU_ldsc_paras$C[c(exposure, outcome), c(exposure, outcome)]
Omega.no.conf <- MRMU_ldsc_paras$Omega[c(exposure, outcome), c(exposure, outcome)]
```

For this example, the extracted matrices are:

```r
Omega.no.conf
```

```text
                 Biomarker_Urate   CAD_UKB
Biomarker_Urate        1.163e-07 2.238e-08
CAD_UKB                2.238e-08 5.970e-08
```

```r
C.no.conf
```

```text
                 Biomarker_Urate CAD_UKB
Biomarker_Urate           1.1120  0.0529
CAD_UKB                   0.0529  1.0296
```

Also reuse the harmonised multi-trait data from Step 1. We only keep the urate
and CAD_UKB columns needed by the no-confounder model:

```r
MRdat.no.conf.all <- data.frame(
  SNP = dat$SNP,
  A1 = dat$A1,
  A2 = dat$A2,
  b.exp = dat[[paste0("b_", exposure)]],
  b.out = dat[[paste0("b_", outcome)]],
  se.exp = dat[[paste0("se_", exposure)]],
  se.out = dat[[paste0("se_", outcome)]],
  pval.exp = dat[[paste0("pval_", exposure)]],
  pval.out = dat[[paste0("pval_", outcome)]],
  L2 = dat$L2
)
```

## 5.2. LD-Clump the Urate Instruments

Use the same IV threshold and local PLINK reference:

```r
MRdat.no.conf <- MRAPSS::clump(
  MRdat.no.conf.all,
  IV.Threshold = p_thresh,
  SNP_col = "SNP",
  pval_col = "pval.exp",
  clump_kb = 1000,
  clump_r2 = 0.001,
  bfile = bfile_prefix,
  plink_bin = plink_bin
)
```

The no-confounder clumping summary was:

| item | value |
|---|---:|
| candidate urate IVs at p <= 5e-5 | 14,090 |
| clumped IVs | 696 |
| removed due to LD/reference absence | 13,394 |

These are the same candidate and clumped IV counts as the MRMU run because both
models now reuse the same Step 1 harmonised multi-trait data.

## 5.3. Fit the No-Confounder Model

Fit the no-confounder model:

```r
fit.no.conf <- MRAPSS::MRAPSS(
  MRdat.no.conf,
  exposure = exposure,
  outcome = outcome,
  C = C.no.conf,
  Omega = Omega.no.conf,
  Cor.SelectionBias = TRUE
)
```

Extract a compact result table:

```r
pvalue.no.conf <- if (!is.null(fit.no.conf$pvalue)) {
  fit.no.conf$pvalue
} else {
  fit.no.conf$pval
}

res.no.conf <- data.frame(
  exposure = exposure,
  outcome = outcome,
  confounders = "None",
  threshold = p_thresh,
  beta = fit.no.conf$beta,
  se = fit.no.conf$beta.se,
  pvalue = pvalue.no.conf,
  nIV = nrow(MRdat.no.conf),
  neffIV = fit.no.conf$pi0 * nrow(MRdat.no.conf),
  pi0 = fit.no.conf$pi0,
  sigma_sq = fit.no.conf$sigma.sq,
  tau_sq = fit.no.conf$tau.sq,
  method = "MRAPSS(no confounders)"
)
```

## Step 5 Result: MRAPSS Without Confounders

The no-confounder MRAPSS result was:

| quantity | value |
|---|---:|
| beta | 0.04490893 |
| standard error | 0.01613808 |
| p-value | 5.3893e-03 |
| clumped IVs | 696 |
| effective IVs | 324.7741 |
| pi0 | 0.4666295 |
| sigma.sq | 4.59646e-06 |
| tau.sq | 3.29462e-07 |
| threshold | 5e-5 |
| method | MRAPSS without confounders |

Without SBP and DBP in the model, the estimated urate effect on CAD_UKB was
larger and statistically significant.

# MRMU vs MRAPSS Without Confounders

The side-by-side comparison is:

| method | confounders | beta | standard error | p-value | IVs | effective IVs | pi0 |
|---|---|---:|---:|---:|---:|---:|---:|
| MRMU | SBP + DBP | 0.01232753 | 0.01526020 | 0.41919 | 696 | 322.2309 | 0.4629754 |
| MRAPSS | none | 0.04490893 | 0.01613808 | 5.3893e-03 | 696 | 324.7741 | 0.4666295 |

The no-confounder analysis gives a larger and statistically significant effect
estimate. After including SBP and DBP in MRMU, the effect estimate is attenuated
and no longer significant, suggesting that blood pressure traits explain part of
the apparent urate-CAD_UKB association in the single-exposure analysis.

# A Complete Script

The complete script is saved as:

```text
package/scripts/run_urate_cad_mrmu_ldsc_example.R
```

Run it from the package root:

```bash
MRMU_EXAMPLE_DIR=/path/to/MRMU_R_github Rscript scripts/run_urate_cad_mrmu_ldsc_example.R
```

The script writes:

| output file | content |
|---|---|
| `urate_cad_mrmu_ldsc_result.csv` | final MRMU result table |
| `urate_cad_mrmu_ldsc_fit.rds` | fitted model object |
| `urate_cad_mrmu_ldsc_input.rds` | clumped MRMU input |
| `MRMU_ldsc_paras_urate_sbp_dbp_cadukb.rds` | estimated `Omega`, `C`, and genetic correlations |
| `ldsc_pairwise_results.csv` | pairwise LDSC summary table |
| `run.log` | full run log |

The no-confounder comparison script is saved as:

```text
package/scripts/run_urate_cad_mrapss_comparison.R
```

Run it from the package root with:

```bash
MRMU_EXAMPLE_DIR=/path/to/MRMU_R_github Rscript scripts/run_urate_cad_mrapss_comparison.R
```

It writes:

| output file | content |
|---|---|
| `urate_cad_mrapss_result.csv` | no-confounder MRAPSS result table |
| `urate_cad_mrmu_vs_mrapss_comparison.csv` | side-by-side MRMU vs MRAPSS comparison |
| `urate_cad_mrapss_fit.rds` | fitted no-confounder model |
| `urate_cad_mrapss_input.rds` | clumped no-confounder input |
| `urate_cad_mrapss_ldsc_paras_from_step2.rds` | no-confounder `Omega/C` submatrices extracted from Step 2 |

# Troubleshooting

If local clumping fails, check:

- `plink_bin` points to an executable macOS/Linux binary for your system.
- `bfile_prefix` omits the `.bed`, `.bim`, and `.fam` suffix.
- all three PLINK reference files exist in the same folder.

If LDSC parameter estimation fails, check:

- `ldscore_dir` contains all `1-22.l2.ldscore.gz` files.
- `ldscore_dir` also contains all `1-22.l2.M_5_50` files.
- the GWAS files contain valid `SNP`, `A1`, `A2`, `Z`, `N`, `chi2`, and `P`
  columns.

For large runs, save intermediate objects such as `dat`, `MRMUInput`, and
`MRMU_ldsc_paras`, so the full pipeline can resume without repeating the most
expensive steps.
