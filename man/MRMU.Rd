\name{MRMU}
\alias{MRMU}
\title{Fit the MR-MU model}
\usage{
MRMU(MRdat = NULL, exposure = NULL, confounders = NULL, outcome = NULL,
  pi0 = NULL, SigmaX = NULL, tau.sq = NULL, C = NULL, Omega = NULL,
  Cor.SelectionBias = TRUE, tol = 1e-8, a = 1, b = 1,
  cut.confounders = TRUE)
}
\description{
Fits the MR-MU variational EM model for Mendelian randomization with one
primary exposure and optional confounder traits.
}
\arguments{
\item{MRdat}{Input list from \code{MRMU_Input()} or an equivalent object.}
\item{exposure}{Name of the primary exposure.}
\item{confounders}{Character vector of confounder trait names.}
\item{outcome}{Name of the outcome trait.}
\item{pi0, SigmaX, tau.sq}{Optional initial parameter values.}
\item{C, Omega}{Sample-structure and polygenic-effect covariance matrices.}
\item{Cor.SelectionBias}{Whether to correct for IV selection bias.}
\item{tol}{Convergence tolerance.}
\item{a, b}{Hyperparameters for the prior on confounder effects.}
\item{cut.confounders}{Whether to remove weak confounder columns.}
}
\value{A list containing effect estimate, standard error, p-value, fitted
variance parameters, posterior quantities, and model metadata.}
