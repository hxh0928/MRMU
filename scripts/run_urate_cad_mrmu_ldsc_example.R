#!/usr/bin/env Rscript

lib <- file.path(getwd(), ".Rlib")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(MRMU)
  library(MRAPSS)
  library(data.table)
})

bundle_dir <- Sys.getenv("MRMU_EXAMPLE_DIR", unset = "path/to/MRMU_R_github")
if (!dir.exists(bundle_dir)) {
  stop("Set MRMU_EXAMPLE_DIR to the downloaded MRMU_R_github data folder.")
}
data_dir <- file.path(bundle_dir, "raw_gwas")
ldscore_dir <- file.path(bundle_dir, "ldscore", "eur_w_ld_chr")
bfile_prefix <- file.path(bundle_dir, "plink_reference", "all_1000G_EUR_Phase3")
plink_bin <- file.path(bundle_dir, "tools", "plink_mac", "plink")

out_dir <- file.path(getwd(), "results", "urate_cad_mrmu_ldsc_clumped")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("Biomarker_Urate", "Metabolic_SBP", "Metabolic_DBP", "CAD_UKB")
data.files <- data.frame(
  trait.name = traits,
  file.dir = file.path(data_dir, traits),
  stringsAsFactors = FALSE
)

p_thresh <- 5e-5
exposure <- "Biomarker_Urate"
confounders <- c("Metabolic_SBP", "Metabolic_DBP")
outcome <- "CAD_UKB"

message("Reading LD scores...")
ld <- rbindlist(lapply(1:22, function(chr) {
  fread(file.path(ldscore_dir, paste0(chr, ".l2.ldscore.gz")))
}), use.names = TRUE)

M <- sum(unlist(lapply(1:22, function(chr) {
  scan(file.path(ldscore_dir, paste0(chr, ".l2.M_5_50")), quiet = TRUE)
})))

message("Harmonising four trait files with real L2...")
dat <- MV_harmonise(data.files, ld)
saveRDS(dat, file.path(out_dir, "harmonised_urate_sbp_dbp_cadukb_with_l2.rds"))

message("Estimating pairwise LDSC parameters...")
ldsc_cache <- file.path(out_dir, "MRMU_ldsc_paras_urate_sbp_dbp_cadukb.rds")
ldsc_pair_file <- file.path(out_dir, "ldsc_pairwise_results.csv")

if (file.exists(ldsc_cache)) {
  MRMU_ldsc_paras <- readRDS(ldsc_cache)
} else {
  Omega <- matrix(0, length(traits), length(traits), dimnames = list(traits, traits))
  C <- matrix(0, length(traits), length(traits), dimnames = list(traits, traits))
  Rg <- matrix(NA_real_, length(traits), length(traits), dimnames = list(traits, traits))
  pair_rows <- list()

  for (i in seq_along(traits)) {
    for (j in i:length(traits)) {
      trait_i <- traits[i]
      trait_j <- traits[j]
      message("LDSC: ", trait_i, " ~ ", trait_j)

      dat1 <- as.data.frame(fread(data.files$file.dir[i]))
      dat2 <- as.data.frame(fread(data.files$file.dir[j]))

      paras <- MRAPSS::est_paras(
        dat1 = dat1,
        dat2 = dat2,
        trait1.name = trait_i,
        trait2.name = trait_j,
        h2.fix.intercept = FALSE,
        LDSC = TRUE,
        ld = ld,
        M = M
      )

      Omega[i, j] <- Omega[j, i] <- paras$Omega[1, 2]
      C[i, j] <- C[j, i] <- paras$C[1, 2]
      Rg[i, j] <- Rg[j, i] <- paras$ldsc_res$rg

      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        trait1 = trait_i,
        trait2 = trait_j,
        rg = paras$ldsc_res$rg,
        rg_se = paras$ldsc_res$rg.se,
        omega = paras$Omega[1, 2],
        C = paras$C[1, 2],
        rho = paras$ldsc_res$I[1, 2],
        rho_se = paras$ldsc_res$I.se[1, 2],
        rhog = paras$ldsc_res$cov[1, 2],
        rhog_se = paras$ldsc_res$cov.se[1, 2]
      )
    }
  }

  diag(C) <- ifelse(diag(C) < 1, 1, diag(C))
  MRMU_ldsc_paras <- list(Omega = Omega, C = C, Rg = Rg)
  pair_results <- rbindlist(pair_rows)
  fwrite(pair_results, ldsc_pair_file)
  saveRDS(MRMU_ldsc_paras, ldsc_cache)
}

message("Running MRMU_Input with local PLINK clumping...")
MRMUInput <- MRMU_Input(
  dat = dat,
  exposure = exposure,
  confounders = confounders,
  outcome = outcome,
  p_thresh = p_thresh,
  bfile = bfile_prefix,
  plink_bin = plink_bin,
  clump_kb = 1000,
  clump_r2 = 0.001
)

traits_for_model <- c(exposure, confounders, outcome)
message("Fitting MRMU with LDSC Omega/C and real L2...")
fit <- MRMU(
  MRdat = MRMUInput,
  exposure = exposure,
  confounders = confounders,
  outcome = outcome,
  Omega = MRMU_ldsc_paras$Omega[traits_for_model, traits_for_model],
  C = MRMU_ldsc_paras$C[traits_for_model, traits_for_model],
  Cor.SelectionBias = TRUE,
  cut.confounders = FALSE,
  tol = 1e-8,
  a = 1,
  b = 1
)

res <- data.frame(
  exposure = exposure,
  outcome = outcome,
  confounders = paste(confounders, collapse = ", "),
  nconfs = fit$nexps,
  threshold = p_thresh,
  beta = unname(fit$beta1),
  se = unname(fit$beta1.se),
  pvalue = fit$beta1.pvalue,
  nIV = fit$nIV,
  neffIV = unname(fit$nvalid),
  pi0 = unname(fit$pi0),
  tau_sq = unname(fit$tau.sq),
  a = fit$a,
  b = fit$b,
  tol = fit$tol,
  method = "MR-MU(LDSC Omega/C, local PLINK clump)",
  stringsAsFactors = FALSE
)

write.csv(res, file.path(out_dir, "urate_cad_mrmu_ldsc_result.csv"), row.names = FALSE)
saveRDS(fit, file.path(out_dir, "urate_cad_mrmu_ldsc_fit.rds"))
saveRDS(MRMUInput, file.path(out_dir, "urate_cad_mrmu_ldsc_input.rds"))

metadata <- list(
  data_files = data.files,
  ldscore_dir = ldscore_dir,
  bfile_prefix = bfile_prefix,
  plink_bin = plink_bin,
  p_thresh = p_thresh,
  n_harmonised = nrow(dat),
  n_iv = fit$nIV,
  M = M
)
saveRDS(metadata, file.path(out_dir, "urate_cad_mrmu_ldsc_metadata.rds"))

message("Wrote LDSC-based results to: ", out_dir)
print(res)
