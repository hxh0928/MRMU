\name{MV_harmonise}
\alias{MV_harmonise}
\title{Harmonise multiple GWAS summary statistic files}
\usage{
MV_harmonise(data.files, ldscore)
}
\description{
Harmonises alleles across multiple GWAS summary statistic files and merges LD
scores.
}
\arguments{
\item{data.files}{Data frame with \code{trait.name} and \code{file.dir}
columns.}
\item{ldscore}{Data frame containing \code{SNP} and \code{L2} columns.}
}
\value{A harmonised data frame.}
