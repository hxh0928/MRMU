#' Fit the MR-MU model
#'
#' Fits the MR-MU variational EM model using summary statistics for one primary
#' exposure, one outcome, and zero or more confounder traits. When no confounder
#' remains after optional filtering, the function falls back to the single-
#' exposure MRAPSS model.
#'
#' @param MRdat A list returned by [MRMU_Input()] or an equivalent object with
#'   matrices `b.exp`, `se.exp`, and `pval.exps`, vectors or one-column matrices
#'   `b.out`, `se.out`, `L2`, and an instrument selection `Threshold`.
#' @param exposure Character string naming the primary exposure trait.
#' @param confounders Character vector naming confounder traits. Use `NULL` for
#'   a single-exposure model.
#' @param outcome Character string naming the outcome trait.
#' @param pi0 Optional initial value for the proportion of selected instruments
#'   with foreground signal. If `NULL`, the algorithm initializes it internally.
#' @param SigmaX Optional initial covariance matrix for exposure and confounder
#'   foreground effects. If `NULL`, it is initialized from `MRdat$b.exp`.
#' @param tau.sq Optional initial foreground outcome-effect variance. If `NULL`,
#'   it is initialized from `MRdat$b.out`.
#' @param C Sample-structure matrix estimated from LDSC. Rows and columns should
#'   be named by `c(exposure, confounders, outcome)`.
#' @param Omega Polygenic-effect covariance matrix estimated from LDSC. Rows and
#'   columns should be named by `c(exposure, confounders, outcome)`.
#' @param Cor.SelectionBias Logical. If `TRUE`, corrects for IV selection using
#'   `MRdat$Threshold`. If `FALSE`, uses `Threshold = 1`.
#' @param tol Numeric convergence tolerance for the variational EM algorithm.
#' @param a,b Hyperparameters for the shrinkage prior on confounder effects.
#' @param cut.confounders Logical. If `TRUE`, removes confounder columns with
#'   very weak genome-wide evidence before fitting the model.
#'
#' @return A list containing the causal estimate `beta1`, standard error
#'   `beta1.se`, p-value `beta1.pvalue`, fitted variance parameters, posterior
#'   quantities, instrument counts, convergence settings, and model metadata.
#'
#' @examples
#' \dontrun{
#' library(MRMU)
#' data("urate_cad_mrmu_iv")
#' data("urate_cad_mrmu_background")
#'
#' fit <- MRMU(
#'   MRdat = urate_cad_mrmu_iv,
#'   exposure = "Biomarker_Urate",
#'   confounders = c("Metabolic_SBP", "Metabolic_DBP"),
#'   outcome = "CAD_UKB",
#'   C = urate_cad_mrmu_background$C,
#'   Omega = urate_cad_mrmu_background$Omega
#' )
#' fit$beta1
#' }
#' @export
MRMU <- function(MRdat = NULL,
                 exposure = NULL,
                 confounders = NULL,
                 outcome = NULL,
                 pi0 = NULL,
                 SigmaX = NULL,
                 tau.sq = NULL,
                 C = NULL,
                 Omega = NULL,
                 Cor.SelectionBias = TRUE,
                 tol = 1e-8,
                 a = 1, b = 1,
                 cut.confounders = TRUE){
  
  if(is.null(MRdat)){
    cat("No MRdat for MR testing")
    return(NULL)
  }
  
  # Selection threshold controls the truncation correction for selected IVs.
  if(!Cor.SelectionBias){
    
    Threshold = 1
    message("Threshold = 1, the model will not correct for selection bias")
    
  }else{
    
    Threshold = unique(stats::na.omit(MRdat$Threshold))
    
    if(length(Threshold) == 0) Threshold = max(MRdat$pval.exp)
    if(length(Threshold) > 1) stop("MRdat$Threshold must contain a single threshold value.")
    
  }
  
  if(nrow(MRdat$b.exp) <= 10) stop(" Not enough IVs.")
  m = nrow(MRdat$b.exp)
  p0 = ncol(MRdat$b.exp)
  cols = 1:p0
  
  # Provide neutral defaults when LDSC background matrices are not supplied.
  if(is.null(C)){
    C = diag((p0+1))}
  
  
  if(is.null(Omega)){
    Omega = matrix(0, (p0+1), (p0+1))}
  

  if(is.null(exposure)| is.null(confounders)){
      exp.names = c(paste0("X", 1:p0))
    }else{
      exp.names = c(exposure, confounders)
    }
  
  if(is.null(outcome)) outcome ="Y"
  
  # Align matrix names before subsetting, allowing either row or column names.
  if(!is.null(rownames(Omega))) colnames(Omega) = row.names(Omega) 

  if(!is.null(colnames(Omega))) rownames(Omega) =  colnames(Omega) 
  
  if(is.null(colnames(Omega)) & is.null(rownames(Omega))) colnames(Omega) = rownames(Omega) = c(exp.names, outcome)
  
  if(!is.null(rownames(C))) colnames(C) = rownames(C) 
  
  if(!is.null(colnames(C))) rownames(C) =  colnames(C) 
  
  if(is.null(colnames(C)) & is.null(rownames(C))) colnames(C) = row.names(C) = c(exp.names, outcome)
  
  if(!all(c(exp.names, outcome) %in%  colnames(Omega)))   stop("Check the columns of your Omega matrix")
  if(!all(c(exp.names, outcome) %in%  colnames(C)))   stop("Check the columns of your C matrix")
  
  
  if(cut.confounders & ! is.null(confounders)){
    # Keep the exposure and confounders with enough strong instruments.
    counts = colSums(MRdat$pval.exps < 5e-08)
    cols = sort(unique(c(1, which(counts > 5))))
  }
  exp.names = exp.names[cols]
  p = length(cols)
  p2 = p - 1
  print(c(exp.names, outcome))
  
  MRdat$b.exp = matrix(MRdat$b.exp[, cols], nrow = m, ncol = p)
  MRdat$se.exp = matrix(MRdat$se.exp[, cols], nrow = m, ncol = p)
  MRdat$pval.exps = matrix(MRdat$pval.exps[, cols], nrow = m, ncol = p)
  Omega = Omega[c(exp.names, outcome), c(exp.names, outcome)]
  C = C[c(exp.names, outcome), c(exp.names, outcome)]

  
  if(!is.null(SigmaX)){
    SigmaX = SigmaX[exp.names, exp.names]
  }  
  
  
  if(p2==0){
    # With no remaining confounder, fit the single-exposure MRAPSS model.
    MRdat = data.frame(b.exp = MRdat$b.exp[,1],
                       b.out = MRdat$b.out,
                       se.exp = MRdat$se.exp[,1],
                       se.out = MRdat$se.out,
                       pval.exp = MRdat$pval.exps[,1],
                       L2 = MRdat$L2,
                       Threshold = MRdat$Threshold)
    colnames(MRdat) = c("b.exp", "b.out", "se.exp", "se.out", "pval.exp", "L2", "Threshold")
    
    fit.MRAPSS = MRAPSS::MRAPSS(MRdat,
                                exposure = exposure,
                                outcome=outcome,
                                Omega = Omega,
                                C=C,
                                Cor.SelectionBias = TRUE)
    
    
    
    
    beta1 = fit.MRAPSS$beta
    beta1.se = fit.MRAPSS$beta.se
    beta1.pvalue =  if(!is.null(fit.MRAPSS$pvalue)) fit.MRAPSS$pvalue else fit.MRAPSS$pval
    tau.sq = fit.MRAPSS$tau.sq
    SigmaX = fit.MRAPSS$sigma.sq
    nIV = m
    nvalid = fit.MRAPSS$pi0*m
    nexps = 1
    pi0 = fit.MRAPSS$pi0
    post = fit.MRAPSS$post
    
    
  }else{
    # Stage 1 fits the null model with the primary causal effect fixed to zero.
    fit_s1 = MRMU_vEMfunc(MRdat = MRdat,
                          beta1 = 0,
                          SigmaX = SigmaX,
                          tau.sq = tau.sq,
                          pi0 = pi0,
                          fix.beta1 = T,
                          fix.tau=F,
                          fix.SigmaX = F,
                          C = C,
                          Omega = Omega,
                          tol = tol,
                          Threshold = Threshold,
                          a=a,b=b )
    
    if(fit_s1$pi0<=1e-04){
      cat("Change hyperparameters \n")
      a=0.5
      b=0.1
      fit_s1 = MRMU_vEMfunc(MRdat = MRdat,
                            beta1 = 0,
                            SigmaX = SigmaX,
                            tau.sq = tau.sq,
                            pi0 = pi0,
                            fix.beta1 = T,
                            fix.tau=F,
                            fix.SigmaX = F,
                            C = C,
                            Omega = Omega,
                            tol = tol,
                            Threshold = Threshold,
                            a=a, b=b)
    }
    
    # Stage 2 releases the primary causal effect and reuses stage 1 posterior
    # values as stable starting points.
    fit_s2 = MRMU_vEMfunc(MRdat = MRdat,
                          beta1 = 0,
                          mu.beta2 = fit_s1$post$mu.beta2,
                          Sigma.beta2 = fit_s1$post$Sigma.beta2,
                          SigmaX = fit_s1$SigmaX,
                          tau.sq = fit_s1$tau.sq,
                          pi0 = fit_s1$pi0,
                          fix.beta1 = F,
                          fix.tau=F,
                          fix.SigmaX = F,
                          C = C,
                          Omega = Omega,
                          tol = tol,
                          Threshold = Threshold,
                          a=a,b=b)
    
    #cat("stage1: likelihood ", fit_s1$log_elbo, "stage2: likelihood ", fit_s2$log_elbo, "\n")
    
    # Likelihood-ratio inference compares stage 2 with the stage 1 null fit.
    LR1 = 2*(fit_s2$log_elbo - fit_s1$log_elbo)
    pvalue1 = pchisq(LR1, 1, lower.tail = F)
    pvalue1 = formatC(pvalue1, format = "e", digits = 4)
    beta1.se = suppressWarnings(abs(fit_s2$beta1/sqrt(LR1)))
    SigmaX = fit_s2$SigmaX
    
    
    cat("***********************************************************\n")
    cat("MR test results of ", exposure , " on ", outcome, ": \n")
    cat("beta1 = ", round(fit_s2$beta1,4), "beta1.se = ", round(beta1.se, 4), "beta1.pvalue = ", pvalue1,  "\n")
    cat("Total No.of IVs:", nrow(MRdat$b.exp), "Effective NO. of IVs:", fit_s2$pi0 * nrow(MRdat$b.exp), "\n")
    cat("***********************************************************\n")
    
    beta1 = fit_s2$beta1
    beta1.se = beta1.se
    beta1.pvalue = pvalue1
    tau.sq = fit_s2$tau.sq
    SigmaX = fit_s2$SigmaX
    nIV = m
    nvalid = fit_s2$pi0*m
    nexps = p2
    pi0 = fit_s2$pi0
    post = fit_s2$post
    colnames(SigmaX) = rownames(SigmaX) = exp.names
    tol=fit_s2$tol
  }
  
  
  return( list(exposure = exposure,
               outcome = outcome,
               beta1 = beta1,
               beta1.se = beta1.se,
               beta1.pvalue = beta1.pvalue,
               tau.sq = tau.sq,
               SigmaX = SigmaX,
               nIV = m,
               nvalid = nvalid,
               nexps = p2,
               pi0 = pi0,
               post = post,
               a=a,
               b=b,
               tol=tol,
               Threshold = Threshold,
               method = "MR-MU"))
  
}
# 
