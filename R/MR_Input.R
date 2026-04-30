comple <- function(allele){
  # Return the Watson-Crick complement for standard nucleotide alleles.
  ifelse(allele == "A","T", ifelse(allele == "T","A", ifelse(allele == "G","C", ifelse(allele == "C","G", allele)) ))
  
}

harmonise <- function(dat, IDs){
  # Merge with the allele reference from the first trait.
  dat = merge(IDs, dat, by="SNP")
  
  # Flip z-scores when the study allele pair is reversed or strand-complemented.
  flip.index = which((dat$A1.x == dat$A2.y & dat$A1.y == dat$A2.x) |
                       (dat$A1.x ==comple(dat$A2.y) & dat$A1.y == comple(dat$A2.x)))
  
  dat[,"A1.y"] = dat[,"A1.x"]
  dat[,"A2.y"] = dat[,"A2.x"]
  dat[flip.index ,"Z"] = -dat[flip.index ,"Z"]
  
  dat = data.frame(SNP = dat$SNP,
                   A1 = dat$A1.x,
                   A2 = dat$A2.x,
                   b = dat$Z/sqrt(dat$N),
                   se = 1/sqrt(dat$N),
                   pval = dat$P)
}

#' Harmonise multiple GWAS summary statistic files
#'
#' Reads multiple GWAS summary statistic files, aligns all traits to the allele
#' coding of the first trait, converts signed z-scores to effect estimates, and
#' merges SNP-level LD scores.
#'
#' @param data.files A data frame with one row per trait. It must contain
#'   `trait.name`, the trait label used in output column names, and `file.dir`,
#'   the path to a GWAS summary-statistic file.
#' @param ldscore A data frame containing at least `SNP` and `L2` columns, where
#'   `L2` is the SNP-level LD score.
#'
#' @details Each GWAS file should contain columns `SNP`, `A1`, `A2`, `Z`, `N`,
#'   and `P`. The function uses the first trait as the allele reference and
#'   keeps SNPs available in all traits and in `ldscore`.
#'
#' @return A harmonised data frame containing `SNP`, `A1`, `A2`, trait-specific
#'   `b_`, `se_`, and `pval_` columns, plus `L2`.
#'
#' @examples
#' \dontrun{
#' data.files <- data.frame(
#'   trait.name = c("Biomarker_Urate", "CAD_UKB"),
#'   file.dir = c("path/to/Biomarker_Urate", "path/to/CAD_UKB")
#' )
#' ldscore <- data.frame(SNP = "rs1", L2 = 1)
#' dat <- MV_harmonise(data.files, ldscore)
#' }
#' @export
MV_harmonise <- function(data.files, ldscore){
  
   trait.names = data.files$trait.name
   
   trait1.dat = data.table::fread(data.files[1, "file.dir"])
   
   IDs = trait1.dat[, c("SNP", "A1", "A2")]
   
   dat = harmonise(trait1.dat , IDs)
   
   colnames(dat) = c("SNP", "A1", "A2", paste0("b_", trait.names[1]), paste0("se_", trait.names[1]), paste0("pval_", trait.names[1]))
   
   cat("Before hamonising, we have",nrow(dat), "SNPs \n")
   
  for (i in 2:length(trait.names)){
    
     trait.dat = data.table::fread(data.files[i, "file.dir"])
     
     dat0 = harmonise(trait.dat, IDs)
    
     colnames(dat0) = c("SNP", "A1", "A2", paste0("b_", trait.names[i]), 
                        paste0("se_", trait.names[i]), paste0("pval_", trait.names[i]))

     dat = merge(dat, dat0, by= c("SNP", "A1", "A2"))
     
     cat("Harmonising with",trait.names[i], ". Remaining", nrow(dat), "SNPs. \n")
   
  }
  
  dat = merge(dat, ldscore[, c("SNP", "L2")], by="SNP")
  return(dat)
  
}
#' Create single-variable MR input data
#'
#' Selects instruments for one exposure, performs LD clumping, and returns a
#' clumped data frame suitable for single-variable MR analyses.
#'
#' @param dat Harmonised data frame from [MV_harmonise()] or an equivalent data
#'   frame containing `SNP`, allele columns, and trait-specific `b_`, `se_`, and
#'   `pval_` columns.
#' @param exposure Character string naming the exposure trait.
#' @param outcome Character string naming the outcome trait.
#' @param p_thresh Instrument selection p-value threshold for the exposure.
#' @param clump_kb Window size in kilobases used by PLINK clumping.
#' @param clump_r2 LD r-squared threshold used by PLINK clumping.
#' @param bfile PLINK reference prefix without `.bed`, `.bim`, or `.fam`.
#'   Set to `NULL` to use the default clumping path from `MRAPSS::clump()`.
#' @param plink_bin Path to the PLINK executable.
#'
#' @return A clumped data frame with selected SNPs and exposure-outcome summary
#'   statistics.
#'
#' @examples
#' \dontrun{
#' svmr_dat <- SVMR_Input(dat, "Biomarker_Urate", "CAD_UKB", bfile = NULL)
#' }
#' @export
SVMR_Input <- function(dat, exposure, outcome, p_thresh=5e-08, clump_kb = 1000,clump_r2 = 0.001,
                       bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
                       plink_bin = "/import/home2/maxhu/plink/plink"){
  
  
  SVMRdat = dat[, names(dat) %in% c("SNP", "A1", "A2", 
                                     paste0("b_", exposure), paste0("se_", exposure), paste0("pval_", exposure),
                                     paste0("b_", outcome), paste0("se_", outcome), paste0("pval_", outcome))]
  
  SVMRdat$pval_selec = SVMRdat[, paste0("pval_", exposure)]
  
  # Select instruments by exposure p-value and remove correlated variants.
  if(is.null(bfile)){
    
    SVMR_Input = MRAPSS::clump(SVMRdat,
                               IV.Threshold = p_thresh,
                               SNP_col = "SNP",
                               pval_col = "pval_selec",
                               clump_kb = clump_kb,
                               clump_r2 = clump_r2) 
    
  }else{
    
    SVMR_Input = MRAPSS::clump(SVMRdat,
                               IV.Threshold = p_thresh,
                               SNP_col = "SNP",
                               pval_col = "pval_selec",
                               clump_kb = clump_kb,
                               clump_r2 = clump_r2,
                               bfile = bfile,
                               plink_bin = plink_bin) 
    
  }
  return(SVMR_Input)

}

#' Create MR-MU input data
#'
#' Selects instruments for the primary exposure, performs LD clumping, and
#' reshapes harmonised summary statistics into the list structure expected by
#' [MRMU()].
#'
#' @param dat Harmonised data frame from [MV_harmonise()] or an equivalent data
#'   frame containing `SNP`, `A1`, `A2`, `L2`, and trait-specific `b_`, `se_`,
#'   and `pval_` columns.
#' @param exposure Character string naming the primary exposure trait.
#' @param confounders Character vector naming confounder traits.
#' @param outcome Character string naming the outcome trait.
#' @param p_thresh Instrument selection p-value threshold for the exposure.
#' @param bfile PLINK reference prefix without `.bed`, `.bim`, or `.fam`.
#'   Set to `NULL` to use the default clumping path from `MRAPSS::clump()`.
#' @param plink_bin Path to the PLINK executable.
#' @param clump_kb Window size in kilobases used by PLINK clumping.
#' @param clump_r2 LD r-squared threshold used by PLINK clumping.
#'
#' @return A list containing exposure names, confounder names, outcome name,
#'   matrices of effect estimates and standard errors, LD scores, exposure
#'   p-values, and the selection threshold.
#'
#' @examples
#' \dontrun{
#' mrmu_dat <- MRMU_Input(
#'   dat,
#'   exposure = "Biomarker_Urate",
#'   confounders = c("Metabolic_SBP", "Metabolic_DBP"),
#'   outcome = "CAD_UKB",
#'   bfile = NULL
#' )
#' }
#' @export
MRMU_Input <- function(dat, exposure, confounders, outcome, p_thresh=5e-05,
                       bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
                       plink_bin = "/import/home2/maxhu/plink/plink",
                       clump_kb = 1000,clump_r2 = 0.001){
  
  dat$pval_selec = dat[, paste0("pval_", exposure)]
  
  # Use the primary exposure to select IVs, then LD-clump the candidates.
  if(is.null(bfile)){
    
    MRMU_IV = MRAPSS::clump(dat,
                            IV.Threshold = p_thresh,
                            SNP_col = "SNP",
                            pval_col = "pval_selec",
                            clump_kb = clump_kb,
                            clump_r2 = clump_r2) 
    
  }else{
    
    MRMU_IV = MRAPSS::clump(dat,
                            IV.Threshold = p_thresh,
                            SNP_col = "SNP",
                            pval_col = "pval_selec",
                            clump_kb = clump_kb,
                            clump_r2 = clump_r2,
                            bfile = bfile,
                            plink_bin = plink_bin) 
    
  }

  # Store exposure and confounder effects as matrices in a stable trait order.
  b.exp = as.matrix(MRMU_IV[, paste0("b_", c(exposure, confounders))])
  b.out = as.matrix(MRMU_IV[, paste0("b_", outcome)])
  se.exp = as.matrix(MRMU_IV[, paste0("se_", c(exposure, confounders))])
  se.out = as.matrix(MRMU_IV[, paste0("se_", outcome)])
  L2 = MRMU_IV$L2
  pval.exps = as.matrix(MRMU_IV[, paste0("pval_", c(exposure, confounders))])
  m = nrow(MRMU_IV)
  p = length(confounders) +1
  
  b.exp = matrix(b.exp, nrow = m, ncol=p)
  colnames(b.exp) = paste0("b_", c(exposure, confounders))
  
  b.out = matrix(b.out, nrow =m, ncol=1)
  colnames(b.out) = paste0("b_", outcome)
  
  se.exp= matrix(se.exp, nrow = m, ncol =p)
  colnames(se.exp) = paste0("se_", c(exposure, confounders))
  
  se.out = matrix(se.out, nrow = m, ncol =1)
  colnames(se.out) = paste0("se_", outcome)
  
  L2 = matrix(L2, ncol=1, nrow=m)
  
  pval.exps = matrix(pval.exps, nrow = m, ncol=p)
  colnames(pval.exps) = paste0("pval_", c(exposure, confounders))
  
  
  MRMU_Input = list(exposure=exposure, 
                    confounders=confounders,
                    outcome=outcome,
                    b.exp = b.exp,
                    b.out = b.out,
                    se.exp= se.exp,
                    se.out = se.out,
                    L2 = L2,
                    pval.exps = pval.exps,
                    Threshold = p_thresh) 
  return(MRMU_Input)
  
}

#' Select instruments for multivariable MR
#'
#' Selects and LD-clumps instruments across multiple exposure traits. The
#' function first identifies instruments for each exposure and then performs a
#' second clumping step to remove duplicates and correlated variants across the
#' pooled instrument set.
#'
#' @param dat Harmonised summary-statistic data frame containing `SNP`, allele
#'   columns, and trait-specific `b_`, `se_`, and `pval_` columns.
#' @param exposures Character vector of exposure trait names.
#' @param p_thresh Instrument selection p-value threshold.
#' @param bfile PLINK reference prefix without `.bed`, `.bim`, or `.fam`.
#'   Set to `NULL` to use the default clumping path from `MRAPSS::clump()`.
#' @param plink_bin Path to the PLINK executable.
#' @param clump_kb Window size in kilobases used by PLINK clumping.
#' @param clump_r2 LD r-squared threshold used by PLINK clumping.
#'
#' @return A clumped data frame containing the union of selected instruments
#'   across exposure traits.
#'
#' @examples
#' \dontrun{
#' ivs <- MVMR_IV(dat, c("Biomarker_Urate", "Metabolic_SBP"), 5e-5, bfile = NULL)
#' }
#' @export
MVMR_IV <- function(dat, exposures, p_thresh,
                    bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
                    plink_bin = "/import/home2/maxhu/plink/plink",
                    clump_kb = 1000,clump_r2 = 0.001){

  MVMRdat = NULL
  
  for (trait.name in exposures){
    cat(trait.name, "\n")
    
    dat$pval_selec = dat[, paste0("pval_", trait.name)]
    
    if(length(which(dat$pval_selec<=5e-08))==0){
      cat("No enough IVs \n")
      next
    } 
    # Select and clump instruments for each exposure separately.
    if(is.null(bfile)){
      
      MVMR_traitIV = MRAPSS::clump(dat,
                                   IV.Threshold = p_thresh,
                                   SNP_col = "SNP",
                                   pval_col = "pval_selec",
                                   clump_kb = clump_kb,
                                   clump_r2 = clump_r2) 
      
    }else{
      
      MVMR_traitIV = MRAPSS::clump(dat,
                              IV.Threshold = p_thresh,
                              SNP_col = "SNP",
                              pval_col = "pval_selec",
                              clump_kb = clump_kb,
                              clump_r2 = clump_r2,
                              bfile = bfile,
                              plink_bin = plink_bin) 
      
    }
    
    
    MVMRdat = rbind(MVMRdat,  MVMR_traitIV)
    
  }
  
  # Re-clump the pooled IV set to remove duplicates and correlated variants.
  
  if(is.null(bfile)){
    
    MVMR_IV = MRAPSS::clump(MVMRdat,
                            IV.Threshold = p_thresh,
                            SNP_col = "SNP",
                            pval_col = "pval_selec",
                            clump_kb = clump_kb,
                            clump_r2 = clump_r2)
    
  }else{
    
    MVMR_IV = MRAPSS::clump(MVMRdat,
                            IV.Threshold = p_thresh,
                            SNP_col = "SNP",
                            pval_col = "pval_selec",
                            clump_kb = clump_kb,
                            clump_r2 = clump_r2,
                            bfile = bfile,
                            plink_bin = plink_bin) 
    
  }
  return(MVMR_IV)
  
}
