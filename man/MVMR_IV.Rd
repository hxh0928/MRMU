\name{MVMR_IV}
\alias{MVMR_IV}
\title{Select instruments for multivariable MR}
\usage{
MVMR_IV(dat, exposures, p_thresh,
  bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
  plink_bin = "/import/home2/maxhu/plink/plink", clump_kb = 1000,
  clump_r2 = 0.001)
}
\description{
Selects and LD-clumps instruments across multiple exposure traits.
}
\arguments{
\item{dat}{Harmonised summary-statistic data frame.}
\item{exposures}{Character vector of exposure trait names.}
\item{p_thresh}{Instrument selection p-value threshold.}
\item{bfile}{PLINK reference prefix; use \code{NULL} for API-based clumping.}
\item{plink_bin}{Path to the PLINK executable.}
\item{clump_kb, clump_r2}{LD clumping parameters.}
}
\value{A clumped data frame.}
