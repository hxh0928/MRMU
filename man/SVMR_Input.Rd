\name{SVMR_Input}
\alias{SVMR_Input}
\title{Create single-variable MR input data}
\usage{
SVMR_Input(dat, exposure, outcome, p_thresh = 5e-8, clump_kb = 1000,
  clump_r2 = 0.001,
  bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
  plink_bin = "/import/home2/maxhu/plink/plink")
}
\description{
Formats and LD-clumps instruments for a single exposure-outcome MR analysis.
}
\arguments{
\item{dat}{Harmonised summary-statistic data frame.}
\item{exposure}{Exposure trait name.}
\item{outcome}{Outcome trait name.}
\item{p_thresh}{Instrument selection p-value threshold.}
\item{clump_kb, clump_r2}{LD clumping parameters.}
\item{bfile}{PLINK reference prefix; use \code{NULL} for API-based clumping.}
\item{plink_bin}{Path to the PLINK executable.}
}
\value{A clumped data frame.}
