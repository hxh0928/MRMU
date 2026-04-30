\name{MRMU_Input}
\alias{MRMU_Input}
\title{Create MR-MU input data}
\usage{
MRMU_Input(
  dat,
  exposure,
  confounders,
  outcome,
  p_thresh = 5e-05,
  bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
  plink_bin = "/import/home2/maxhu/plink/plink",
  clump_kb = 1000,
  clump_r2 = 0.001
)
}
\arguments{
\item{dat}{Harmonised data frame from \code{\link{MV_harmonise}} or an
equivalent data frame containing \code{SNP}, \code{A1}, \code{A2}, \code{L2},
and trait-specific \code{b_}, \code{se_}, and \code{pval_} columns.}

\item{exposure}{Character string naming the primary exposure trait.}

\item{confounders}{Character vector naming confounder traits.}

\item{outcome}{Character string naming the outcome trait.}

\item{p_thresh}{Instrument selection p-value threshold for the exposure.}

\item{bfile}{PLINK reference prefix without \code{.bed}, \code{.bim}, or
\code{.fam}. Set to \code{NULL} to use the default clumping path from
\code{MRAPSS::clump()}.}

\item{plink_bin}{Path to the PLINK executable.}

\item{clump_kb}{Window size in kilobases used by PLINK clumping.}

\item{clump_r2}{LD r-squared threshold used by PLINK clumping.}
}
\value{
A list containing:
\describe{
  \item{\code{exposure}}{Primary exposure name.}
  \item{\code{confounders}}{Confounder trait names.}
  \item{\code{outcome}}{Outcome trait name.}
  \item{\code{b.exp}}{Matrix of SNP-exposure and SNP-confounder effects.}
  \item{\code{b.out}}{Matrix of SNP-outcome effects.}
  \item{\code{se.exp}}{Matrix of exposure/confounder standard errors.}
  \item{\code{se.out}}{Matrix of outcome standard errors.}
  \item{\code{L2}}{Matrix of SNP-level LD scores.}
  \item{\code{pval.exps}}{Matrix of exposure/confounder p-values.}
  \item{\code{Threshold}}{Instrument selection threshold.}
}
}
\description{
Selects instruments for the primary exposure, performs LD clumping, and
reshapes harmonised summary statistics into the list structure expected by
\code{\link{MRMU}}.
}
\examples{
\dontrun{
mrmu_dat <- MRMU_Input(
  dat,
  exposure = "Biomarker_Urate",
  confounders = c("Metabolic_SBP", "Metabolic_DBP"),
  outcome = "CAD_UKB",
  bfile = NULL
)
}
}
