\name{MRMU}
\alias{MRMU}
\title{Fit the MR-MU model}
\usage{
MRMU(
  MRdat = NULL,
  exposure = NULL,
  confounders = NULL,
  outcome = NULL,
  pi0 = NULL,
  SigmaX = NULL,
  tau.sq = NULL,
  C = NULL,
  Omega = NULL,
  Cor.SelectionBias = TRUE,
  tol = 1e-08,
  a = 1,
  b = 1,
  cut.confounders = TRUE
)
}
\arguments{
\item{MRdat}{A list returned by \code{\link{MRMU_Input}} or an equivalent
object with matrices \code{b.exp}, \code{se.exp}, and \code{pval.exps},
vectors or one-column matrices \code{b.out}, \code{se.out}, \code{L2}, and an
instrument selection \code{Threshold}.}

\item{exposure}{Character string naming the primary exposure trait.}

\item{confounders}{Character vector naming confounder traits. Use \code{NULL}
for a single-exposure model.}

\item{outcome}{Character string naming the outcome trait.}

\item{pi0}{Optional initial value for the proportion of selected instruments
with foreground signal. If \code{NULL}, the algorithm initializes it
internally.}

\item{SigmaX}{Optional initial covariance matrix for exposure and confounder
foreground effects. If \code{NULL}, it is initialized from
\code{MRdat$b.exp}.}

\item{tau.sq}{Optional initial foreground outcome-effect variance. If
\code{NULL}, it is initialized from \code{MRdat$b.out}.}

\item{C}{Sample-structure matrix estimated from LDSC. Rows and columns should
be named by \code{c(exposure, confounders, outcome)}.}

\item{Omega}{Polygenic-effect covariance matrix estimated from LDSC. Rows and
columns should be named by \code{c(exposure, confounders, outcome)}.}

\item{Cor.SelectionBias}{Logical. If \code{TRUE}, corrects for IV selection
using \code{MRdat$Threshold}. If \code{FALSE}, uses \code{Threshold = 1}.}

\item{tol}{Numeric convergence tolerance for the variational EM algorithm.}

\item{a, b}{Hyperparameters for the shrinkage prior on confounder effects.}

\item{cut.confounders}{Logical. If \code{TRUE}, removes confounder columns with
very weak genome-wide evidence before fitting the model.}
}
\value{
A list with the following elements:
\describe{
  \item{\code{exposure}}{Primary exposure name.}
  \item{\code{outcome}}{Outcome name.}
  \item{\code{beta1}}{Estimated primary causal effect.}
  \item{\code{beta1.se}}{Standard error of \code{beta1}.}
  \item{\code{beta1.pvalue}}{Likelihood-ratio p-value for the primary effect.}
  \item{\code{tau.sq}}{Estimated foreground outcome-effect variance.}
  \item{\code{SigmaX}}{Estimated covariance matrix for exposure and confounder
  foreground effects.}
  \item{\code{nIV}}{Number of selected instruments.}
  \item{\code{nvalid}}{Effective number of foreground instruments,
  \code{pi0 * nIV}.}
  \item{\code{nexps}}{Number of confounder effects included in the fitted
  model.}
  \item{\code{pi0}}{Estimated foreground mixture proportion.}
  \item{\code{post}}{Posterior quantities from the variational EM algorithm.}
  \item{\code{a}, \code{b}}{Hyperparameters used by the model.}
  \item{\code{tol}}{Convergence tolerance used by the model.}
  \item{\code{Threshold}}{Instrument selection threshold used for selection
  correction.}
  \item{\code{method}}{Method label, \code{"MR-MU"}.}
}
}
\description{
Fits the MR-MU variational EM model using summary statistics for one primary
exposure, one outcome, and zero or more confounder traits. When no confounder
remains after optional filtering, the function falls back to the single-
exposure MRAPSS model.
}
\examples{
\dontrun{
library(MRMU)
data("urate_cad_mrmu_iv")
data("urate_cad_mrmu_background")

fit <- MRMU(
  MRdat = urate_cad_mrmu_iv,
  exposure = "Biomarker_Urate",
  confounders = c("Metabolic_SBP", "Metabolic_DBP"),
  outcome = "CAD_UKB",
  C = urate_cad_mrmu_background$C,
  Omega = urate_cad_mrmu_background$Omega
)
fit$beta1
}
}
