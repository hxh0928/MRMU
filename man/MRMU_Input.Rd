\name{MRMU_Input}
\alias{MRMU_Input}
\title{Create MR-MU input data}
\usage{
MRMU_Input(dat, exposure, confounders, outcome, p_thresh = 5e-5,
  bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
  plink_bin = "/import/home2/maxhu/plink/plink", clump_kb = 1000,
  clump_r2 = 0.001)
}
\description{
Selects and LD-clumps instruments, then formats harmonised summary statistics
for \code{MRMU()}.
}
\arguments{
\item{dat}{Harmonised data frame containing \code{b_}, \code{se_}, and
\code{pval_} columns for each trait.}
\item{exposure}{Primary exposure trait name.}
\item{confounders}{Character vector of confounder trait names.}
\item{outcome}{Outcome trait name.}
\item{p_thresh}{Instrument selection p-value threshold.}
\item{bfile}{PLINK reference prefix; use \code{NULL} for API-based clumping.}
\item{plink_bin}{Path to the PLINK executable.}
\item{clump_kb, clump_r2}{LD clumping parameters.}
}
\value{A list suitable for \code{MRMU()}.}
