\name{MV_harmonise}
\alias{MV_harmonise}
\title{Harmonise multiple GWAS summary statistic files}
\usage{
MV_harmonise(data.files, ldscore)
}
\arguments{
\item{data.files}{A data frame with one row per trait. It must contain
\code{trait.name}, the trait label used in output column names, and
\code{file.dir}, the path to a GWAS summary-statistic file.}

\item{ldscore}{A data frame containing at least \code{SNP} and \code{L2}
columns, where \code{L2} is the SNP-level LD score.}
}
\value{
A harmonised data frame containing \code{SNP}, \code{A1}, \code{A2},
trait-specific \code{b_}, \code{se_}, and \code{pval_} columns, plus
\code{L2}.
}
\description{
Reads multiple GWAS summary statistic files, aligns all traits to the allele
coding of the first trait, converts signed z-scores to effect estimates, and
merges SNP-level LD scores.
}
\details{
Each GWAS file should contain columns \code{SNP}, \code{A1}, \code{A2},
\code{Z}, \code{N}, and \code{P}. The function uses the first trait as the
allele reference and keeps SNPs available in all traits and in \code{ldscore}.
}
\examples{
\dontrun{
data.files <- data.frame(
  trait.name = c("Biomarker_Urate", "CAD_UKB"),
  file.dir = c("path/to/Biomarker_Urate", "path/to/CAD_UKB")
)
ldscore <- data.frame(SNP = "rs1", L2 = 1)
dat <- MV_harmonise(data.files, ldscore)
}
}
