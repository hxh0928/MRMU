#' #' @title  A function for implementing MR-APSS.
#' #' @description MR-APSS: a unified approach to Mendelian Randomization accounting for pleiotropy and sample structure using genome-wide summary statistics.
#' #' MA-APSS uses a variantional EM algorithm for estimation of parameters.
#' #' MR-APSS uses likelihood ratio test for inference.
#' #' 
#' #' @param MRdat  MRdat frame at least contain the following varaibles: b.exp b.out se.exp se.out L2 Threshold. L2:LD score, Threshold: modified IV selection threshold for correction of selection bias
#' #' @param exposure exposure name
#' #' @param outcome   outcome name
#' #' @param pi0 initial value for pi0, default `NULL` will use the default initialize procedure.
#' #' @param sigma.sq initial value for sigma.sq , default `NULL`will use the default initialize procedure.
#' #' @param tau.sq initial value for tau.sq , default `NULL` will use the default initialize procedure.
#' #' @param C  the estimated C matrix capturing the effects of sample structure. default `diag(2)`.
#' #' @param Omega  the estimated variance-covariance matrix of polygenic effects. default `matrix(0,2,2)`.
#' #' @param tol     tolerence, default '1e-08'
#' #' @param Cor.SelectionBias   Whether use the selection Threshold for correction of selection bias. If FALSE, the model won't correct for selection bias.
#' #' @param ELBO     Whether check the evidence lower bound or not, if `FALSE`, check the maximum likelihood instead. default `FALSE`.
#' #' 
#' #' @return a list with the following elements:
#' #' \describe{
#' #' \item{MRdat: }{Input MRdat frame}
#' #' \item{exposure: }{exposure of interest}
#' #' \item{outcome: }{outcome of interest}
#' #' \item{beta: }{causal effect estimate}
#' #' \item{beta.se: }{standard error}
#' #' \item{pval: }{p-value}
#' #' \item{sigma.sq: }{variance of forground exposure effect}
#' #' \item{tau.sq: }{variance of forground outcome effect}
#' #' \item{pi0: }{The probability of a SNP with forground signal after selection}
#' #' \item{post: }{Posterior estimates of latent varaibles}
#' #' \item{method: }{"MR-APSS"}
#' #' }
#' #' 
#' #' @examples
#' #' library(MRAPSS)
#' #' exposure = "BMI"
#' #' outcome = "T2D"
#' #' Threshold = 5e-05  # IV selection Threshold
#' #' MRdat(C)
#' #' MRdat(Omega)
#' #' MRdat(MRdat)
#' #' MRres = MRAPSS(MRdat,
#' #'                exposure = "BMI",
#' #'                outcome = "T2D",
#' #'                C = C,
#' #'                Omega =  Omega ,
#' #'                Cor.SelectionBias = TRUE)
#' #' MRplot(MRres, exposure = "BMI", outcome = "T2D")
#' #' @export
#' #' 
#' 

MRMU <- function(MRdat = NULL,
                 exposure = NULL,
                 confounders=NULL,
                 outcome = NULL,
                 pi0 = NULL,
                 SigmaX = NULL,
                 tau.sq = NULL,
                 C = NULL,
                 Omega = NULL,
                 Cor.SelectionBias = TRUE,
                 tol = 1e-8,
                 a =1, b=1,
                 cut.confounders = TRUE){
  
  if(is.null(MRdat)){
    cat("No MRdat for MR testing")
    return(NULL)
  }
  
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
  
  if(!is.null(rownames(Omega))) colnames(Omega) = row.names(Omega) 

  if(!is.null(colnames(Omega))) rownames(Omega) =  colnames(Omega) 
  
  if(is.null(colnames(Omega)) & is.null(rownames(Omega))) colnames(Omega) = rownames(Omega) = c(exp.names, outcome)
  
  if(!is.null(rownames(C))) colnames(C) = rownames(C) 
  
  if(!is.null(colnames(C))) rownames(C) =  colnames(C) 
  
  if(is.null(colnames(C)) & is.null(rownames(C))) colnames(C) = row.names(C) = c(exp.names, outcome)
  
  if(!all(c(exp.names, outcome) %in%  colnames(Omega)))   stop("Check the columns of your Omega matrix")
  if(!all(c(exp.names, outcome) %in%  colnames(C)))   stop("Check the columns of your C matrix")
  
  
  if(cut.confounders & ! is.null(confounders)){
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
    # run MR-APSS
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
    ## stage 1
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
    
    # Inference
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
