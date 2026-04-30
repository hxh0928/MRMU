\name{MVMR_IV}
\alias{MVMR_IV}
\title{Select instruments for multivariable MR}
\usage{
MVMR_IV(
  dat,
  exposures,
  p_thresh,
  bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
  plink_bin = "/import/home2/maxhu/plink/plink",
  clump_kb = 1000,
  clump_r2 = 0.001
)
}
\arguments{
\item{dat}{Harmonised summary-statistic data frame containing \code{SNP},
allele columns, and trait-specific \code{b_}, \code{se_}, and \code{pval_}
columns.}

\item{exposures}{Character vector of exposure trait names.}

\item{p_thresh}{Instrument selection p-value threshold.}

\item{bfile}{PLINK reference prefix without \code{.bed}, \code{.bim}, or
\code{.fam}. Set to \code{NULL} to use the default clumping path from
\code{MRAPSS::clump()}.}

\item{plink_bin}{Path to the PLINK executable.}

\item{clump_kb}{Window size in kilobases used by PLINK clumping.}

\item{clump_r2}{LD r-squared threshold used by PLINK clumping.}
}
\value{
A clumped data frame containing the union of selected instruments across
exposure traits.
}
\description{
Selects and LD-clumps instruments across multiple exposure traits. The function
first identifies instruments for each exposure and then performs a second
clumping step to remove duplicates and correlated variants across the pooled
instrument set.
}
\examples{
\dontrun{
ivs <- MVMR_IV(dat, c("Biomarker_Urate", "Metabolic_SBP"), 5e-5, bfile = NULL)
}
}
