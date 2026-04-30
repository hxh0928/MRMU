\name{est_MRMU_ldsc_paras}
\alias{est_MRMU_ldsc_paras}
\title{Estimate MR-MU LDSC parameters}
\usage{
est_MRMU_ldsc_paras(data.files,
  ldscore.files = "/import/home/share/xhu/database/1KG/eur_w_ld_chr")
}
\description{
Estimates pairwise LDSC-based nuisance parameters used by MR-MU.
}
\arguments{
\item{data.files}{Data frame with \code{trait.name} and \code{file.dir}
columns.}
\item{ldscore.files}{Directory containing LD score reference files.}
}
\value{A list containing LDSC results, genetic correlations, \code{Omega}, and
\code{C}.}
