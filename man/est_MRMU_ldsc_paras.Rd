\name{est_MRMU_ldsc_paras}
\alias{est_MRMU_ldsc_paras}
\title{Estimate MR-MU LDSC parameters}
\usage{
est_MRMU_ldsc_paras(
  data.files,
  ldscore.files = "/import/home/share/xhu/database/1KG/eur_w_ld_chr"
)
}
\arguments{
\item{data.files}{A data frame with columns \code{trait.name} and
\code{file.dir}. \code{trait.name} gives the trait label and \code{file.dir}
points to a GWAS summary-statistic file with \code{SNP}, \code{A1},
\code{A2}, \code{Z}, \code{N}, and \code{P} columns.}

\item{ldscore.files}{Directory containing LD score reference files in the
format expected by \code{MRAPSS::est_paras()}.}
}
\value{
A list with four elements:
\describe{
  \item{\code{LDSC.res}}{Pairwise LDSC summary table.}
  \item{\code{gc}}{Matrix of pairwise genetic correlations.}
  \item{\code{Omega}}{Matrix of polygenic-effect covariance parameters.}
  \item{\code{C}}{Matrix of sample-structure parameters.}
}
}
\description{
Estimates pairwise LDSC-based nuisance parameters for all traits listed in
\code{data.files}. The resulting \code{Omega} and \code{C} matrices can be
passed directly to \code{\link{MRMU}} after harmonising and clumping
instruments.
}
\examples{
\dontrun{
paras <- est_MRMU_ldsc_paras(data.files, "path/to/eur_w_ld_chr")
paras$Omega
paras$C
}
}
