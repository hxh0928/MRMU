# MRMU

MRMU is an R package for Mendelian randomization analyses with multiple
exposures/confounders using GWAS summary statistics.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/MRMU")
```

The package depends on `MRAPSS`, `Matrix`, `expm`, and `data.table`.
`MRAPSS` is available from the MRC IEU r-universe repository.

## Basic workflow

```r
library(MRMU)

# Harmonise GWAS summary statistics across traits
dat <- MV_harmonise(data.files, ldscore)

# Build MR-MU input data after thresholding and LD clumping
mrdat <- MRMU_Input(
  dat = dat,
  exposure = "BMI",
  confounders = c("Trait1", "Trait2"),
  outcome = "T2D",
  p_thresh = 5e-5,
  bfile = NULL
)

# Fit MR-MU
fit <- MRMU(
  MRdat = mrdat,
  exposure = "BMI",
  confounders = c("Trait1", "Trait2"),
  outcome = "T2D",
  C = C,
  Omega = Omega
)
```

If PLINK-based clumping is needed, pass `bfile` and `plink_bin` to
`MRMU_Input()`, `SVMR_Input()`, or `MVMR_IV()`.

## Input summary-statistic format

`MV_harmonise()` expects each GWAS file listed in `data.files$file.dir` to contain
at least these columns:

- `SNP`
- `A1`
- `A2`
- `Z`
- `N`
- `P`

`data.files` should also contain `trait.name`, which is used to name the output
columns.
