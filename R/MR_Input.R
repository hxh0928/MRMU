comple <- function(allele){
  
  ifelse(allele == "A","T", ifelse(allele == "T","A", ifelse(allele == "G","C", ifelse(allele == "C","G", allele)) ))
  
}

harmonise <- function(dat, IDs){
  # merge datasets
  dat = merge(IDs, dat, by="SNP")
  
  #"Harmonize the direction of SNP effects of exposure and outcome"
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



SVMR_Input <- function(dat, exposure, outcome, p_thresh=5e-08, clump_kb = 1000,clump_r2 = 0.001,
                       bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
                       plink_bin = "/import/home2/maxhu/plink/plink"){
  
  
  SVMRdat = dat[, names(dat) %in% c("SNP", "A1", "A2", 
                                     paste0("b_", exposure), paste0("se_", exposure), paste0("pval_", exposure),
                                     paste0("b_", outcome), paste0("se_", outcome), paste0("pval_", outcome))]
  
  SVMRdat$pval_selec = SVMRdat[, paste0("pval_", exposure)]
  
  # p value thresholding and LD clumping
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


MRMU_Input <- function(dat, exposure, confounders, outcome, p_thresh=5e-05,
                       bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
                       plink_bin = "/import/home2/maxhu/plink/plink",
                       clump_kb = 1000,clump_r2 = 0.001){
  
  dat$pval_selec = dat[, paste0("pval_", exposure)]
  
  # p value thresholding and LD clumping
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
    # p value thresholding and LD clumping for each exposure
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
  
  # p value thresholding and LD clumping to remove duplicate IVs
  
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
