#' Estimate MR-MU LDSC parameters
#'
#' Estimates pairwise LDSC-based nuisance parameters for all traits listed in
#' `data.files`. The resulting `Omega` and `C` matrices can be passed directly
#' to [MRMU()] after harmonising and clumping instruments.
#'
#' @param data.files A data frame with columns `trait.name` and `file.dir`.
#'   `trait.name` gives the trait label and `file.dir` points to a GWAS
#'   summary-statistic file with `SNP`, `A1`, `A2`, `Z`, `N`, and `P` columns.
#' @param ldscore.files Directory containing LD score reference files in the
#'   format expected by `MRAPSS::est_paras()`.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{`LDSC.res`}{Pairwise LDSC summary table.}
#'   \item{`gc`}{Matrix of pairwise genetic correlations.}
#'   \item{`Omega`}{Matrix of polygenic-effect covariance parameters.}
#'   \item{`C`}{Matrix of sample-structure parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' paras <- est_MRMU_ldsc_paras(data.files, "path/to/eur_w_ld_chr")
#' paras$Omega
#' paras$C
#' }
#' @export
est_MRMU_ldsc_paras <- function(data.files, ldscore.files="/import/home/share/xhu/database/1KG/eur_w_ld_chr"){
  
  traits = data.files$trait.name
  
  # Allocate square matrices using the trait order supplied by data.files.
  Omega = matrix(0,length(traits), length(traits))
  C = matrix(0,length(traits), length(traits))
  gc =  matrix(0,length(traits), length(traits))
  
  LDSC.res = NULL
  for(i in 1:length(traits)){
    
    for(j in i: length(traits)){
      
      trait1 = data.table::fread(data.files[i, "file.dir"])
      trait2 = data.table::fread(data.files[j, "file.dir"])
      
      # Estimate pairwise LDSC parameters for traits i and j.
      paras = MRAPSS::est_paras(dat1 = trait1,
                                dat2 = trait2,
                                trait1.name = traits[i],
                                trait2.name = traits[j],
                                h2.fix.intercept = F,
                                LDSC = T,
                                ldscore.dir = ldscore.files)
      
      Omega[i,j] = Omega[j,i] = paras$Omega[1,2]
      
      C[i,j] = C[j,i] = paras$C[1,2]
      
      gc[i,j] =  gc[j,i] = paras$ldsc_res$rg
      
      LDSC.res = rbind(LDSC.res, data.frame(traits[i],
                                            traits[j],
                                            rg = paras$ldsc_res$rg,
                                            rg.se = paras$ldsc_res$rg.se,
                                            c = paras$ldsc_res$I[1,2],
                                            c.se = paras$ldsc_res$I.se[1,2],
                                            rhog = paras$ldsc_res$cov[1,2],
                                            rhog.se = paras$ldsc_res$cov.se[1,2]))
    }
  }
  
  return(list(LDSC.res=LDSC.res, gc=gc, Omega = Omega, C=C))
}
