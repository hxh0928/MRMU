\docType{data}
\name{urate_cad_ldsc_parameters}
\alias{urate_cad_ldsc_parameters}
\title{LDSC Parameters for the Urate-CAD MRMU Example}
\format{
A list with three named matrices:
\describe{
  \item{\code{Omega}}{Estimated variance-covariance matrix of polygenic effects.}
  \item{\code{C}}{Estimated sample-structure matrix.}
  \item{\code{Rg}}{Estimated genetic-correlation matrix.}
}
All matrices are named by the four traits \code{Biomarker_Urate},
\code{Metabolic_SBP}, \code{Metabolic_DBP}, and \code{CAD_UKB}.
}
\usage{
data(urate_cad_ldsc_parameters)
}
\description{
LDSC-estimated parameters used by the package demonstration for the
urate-CAD MRMU analysis with SBP and DBP as confounders.
}
\details{
The parameters were estimated pairwise from EUR LD score files and aligned to
the same trait order used by \code{urate_cad_mrmu_iv}. Pass
\code{urate_cad_ldsc_parameters$C} and
\code{urate_cad_ldsc_parameters$Omega} directly to \code{\link{MRMU}}.
}
\keyword{datasets}
