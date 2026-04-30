\docType{data}
\name{urate_cad_mrmu_iv}
\alias{urate_cad_mrmu_iv}
\title{Clumped IV Data for the Urate-CAD MRMU Example}
\format{
A list in the format expected by \code{\link{MRMU}}:
\describe{
  \item{\code{exposure}}{Exposure trait name, \code{"Biomarker_Urate"}.}
  \item{\code{confounders}}{Confounder trait names, \code{"Metabolic_SBP"} and \code{"Metabolic_DBP"}.}
  \item{\code{outcome}}{Outcome trait name, \code{"CAD_UKB"}.}
  \item{\code{b.exp}}{Matrix of SNP-exposure and SNP-confounder effect estimates.}
  \item{\code{b.out}}{Vector of SNP-outcome effect estimates.}
  \item{\code{se.exp}}{Matrix of standard errors for exposure and confounder effects.}
  \item{\code{se.out}}{Vector of standard errors for outcome effects.}
  \item{\code{L2}}{SNP-level LD scores.}
  \item{\code{pval.exps}}{Matrix of exposure and confounder p-values.}
  \item{\code{Threshold}}{Instrument selection threshold, \code{5e-5}.}
}
}
\usage{
data(urate_cad_mrmu_iv)
}
\description{
Clumped instrumental-variable input data for the package demonstration. The
example estimates the causal effect of urate on coronary artery disease while
adjusting for systolic and diastolic blood pressure traits as confounders.
}
\details{
The instruments were selected using \code{p <= 5e-5} for \code{Biomarker_Urate}
and LD clumped with a 1000 Genomes EUR reference panel. The object contains
696 clumped instruments and is intentionally small enough to ship with the
package.
}
\keyword{datasets}
