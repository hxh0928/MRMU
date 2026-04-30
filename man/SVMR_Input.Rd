\name{SVMR_Input}
\alias{SVMR_Input}
\title{Create single-variable MR input data}
\usage{
SVMR_Input(
  dat,
  exposure,
  outcome,
  p_thresh = 5e-08,
  clump_kb = 1000,
  clump_r2 = 0.001,
  bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
  plink_bin = "/import/home2/maxhu/plink/plink"
)
}
\arguments{
\item{dat}{Harmonised data frame from \code{\link{MV_harmonise}} or an
equivalent data frame containing \code{SNP}, allele columns, and trait-specific
\code{b_}, \code{se_}, and \code{pval_} columns.}

\item{exposure}{Character string naming the exposure trait.}

\item{outcome}{Character string naming the outcome trait.}

\item{p_thresh}{Instrument selection p-value threshold for the exposure.}

\item{clump_kb}{Window size in kilobases used by PLINK clumping.}

\item{clump_r2}{LD r-squared threshold used by PLINK clumping.}

\item{bfile}{PLINK reference prefix without \code{.bed}, \code{.bim}, or
\code{.fam}. Set to \code{NULL} to use the default clumping path from
\code{MRAPSS::clump()}.}

\item{plink_bin}{Path to the PLINK executable.}
}
\value{
A clumped data frame with selected SNPs and exposure-outcome summary
statistics.
}
\description{
Selects instruments for one exposure, performs LD clumping, and returns a
clumped data frame suitable for single-variable MR analyses.
}
\examples{
\dontrun{
svmr_dat <- SVMR_Input(dat, "Biomarker_Urate", "CAD_UKB", bfile = NULL)
}
}
