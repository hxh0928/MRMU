#!/usr/bin/env Rscript

lib <- file.path(getwd(), ".Rlib")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(MRAPSS)
  library(data.table)
})

bundle_dir <- Sys.getenv("MRMU_EXAMPLE_DIR", unset = "path/to/MRMU_R_github")
if (!dir.exists(bundle_dir)) {
  stop("Set MRMU_EXAMPLE_DIR to the downloaded MRMU_R_github data folder.")
}
data_dir <- file.path(bundle_dir, "raw_gwas")
bfile_prefix <- file.path(bundle_dir, "plink_reference", "all_1000G_EUR_Phase3")
plink_bin <- file.path(bundle_dir, "tools", "plink_mac", "plink")

out_dir <- file.path(getwd(), "results", "urate_cad_mrapss_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exposure <- "Biomarker_Urate"
outcome <- "CAD_UKB"
p_thresh <- 5e-5

message("Loading Step 1 harmonised data and Step 2 LDSC parameters...")
mrmu_dir <- file.path(getwd(), "results", "urate_cad_mrmu_ldsc_clumped")
dat <- readRDS(file.path(mrmu_dir, "harmonised_urate_sbp_dbp_cadukb_with_l2.rds"))
MRMU_ldsc_paras <- readRDS(file.path(mrmu_dir, "MRMU_ldsc_paras_urate_sbp_dbp_cadukb.rds"))

message("Constructing no-confounder MRAPSS input from Step 1 output...")
MRdat_all <- data.frame(
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

message("Running local PLINK clumping...")
MRdat <- MRAPSS::clump(
  MRdat_all,
  IV.Threshold = p_thresh,
  SNP_col = "SNP",
  pval_col = "pval.exp",
  clump_kb = 1000,
  clump_r2 = 0.001,
  bfile = bfile_prefix,
  plink_bin = plink_bin
)

message("Fitting MRAPSS without confounders...")
fit <- MRAPSS::MRAPSS(
  MRdat,
  exposure = exposure,
  outcome = outcome,
  C = MRMU_ldsc_paras$C[c(exposure, outcome), c(exposure, outcome)],
  Omega = MRMU_ldsc_paras$Omega[c(exposure, outcome), c(exposure, outcome)],
  Cor.SelectionBias = TRUE
)

pvalue <- if (!is.null(fit$pvalue)) fit$pvalue else fit$pval
res <- data.frame(
  exposure = exposure,
  outcome = outcome,
  confounders = "None",
  threshold = p_thresh,
  beta = unname(fit$beta),
  se = unname(fit$beta.se),
  pvalue = pvalue,
  nIV = nrow(MRdat),
  neffIV = unname(fit$pi0 * nrow(MRdat)),
  pi0 = unname(fit$pi0),
  sigma_sq = unname(fit$sigma.sq),
  tau_sq = unname(fit$tau.sq),
  method = "MRAPSS(no confounders)",
  stringsAsFactors = FALSE
)

write.csv(res, file.path(out_dir, "urate_cad_mrapss_result.csv"), row.names = FALSE)
saveRDS(fit, file.path(out_dir, "urate_cad_mrapss_fit.rds"))
saveRDS(MRdat, file.path(out_dir, "urate_cad_mrapss_input.rds"))
saveRDS(
  list(
    C = MRMU_ldsc_paras$C[c(exposure, outcome), c(exposure, outcome)],
    Omega = MRMU_ldsc_paras$Omega[c(exposure, outcome), c(exposure, outcome)]
  ),
  file.path(out_dir, "urate_cad_mrapss_ldsc_paras_from_step2.rds")
)

mrmu_file <- file.path(getwd(), "results", "urate_cad_mrmu_ldsc_clumped",
                       "urate_cad_mrmu_ldsc_result.csv")
if (file.exists(mrmu_file)) {
  mrmu <- read.csv(mrmu_file)
  comparison <- rbind(
    data.frame(
      method = "MRMU(SBP+DBP confounders)",
      beta = mrmu$beta,
      se = mrmu$se,
      pvalue = mrmu$pvalue,
      nIV = mrmu$nIV,
      neffIV = mrmu$neffIV,
      pi0 = mrmu$pi0
    ),
    data.frame(
      method = "MRAPSS(no confounders)",
      beta = res$beta,
      se = res$se,
      pvalue = res$pvalue,
      nIV = res$nIV,
      neffIV = res$neffIV,
      pi0 = res$pi0
    )
  )
  write.csv(comparison, file.path(out_dir, "urate_cad_mrmu_vs_mrapss_comparison.csv"),
            row.names = FALSE)
  print(comparison)
}

message("Wrote MRAPSS comparison results to: ", out_dir)
print(res)
