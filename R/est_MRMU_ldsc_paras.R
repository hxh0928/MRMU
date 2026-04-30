est_MRMU_ldsc_paras <- function(data.files, ldscore.files="/import/home/share/xhu/database/1KG/eur_w_ld_chr"){
  
  traits = data.files$trait.name
  
  Omega = matrix(0,length(traits), length(traits))
  C = matrix(0,length(traits), length(traits))
  gc =  matrix(0,length(traits), length(traits))
  
  LDSC.res = NULL
  for(i in 1:length(traits)){
    
    for(j in i: length(traits)){
      
      trait1 = data.table::fread(data.files[i, "file.dir"])
      trait2 = data.table::fread(data.files[j, "file.dir"])
      
      # estimate parameters using LDSC
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

