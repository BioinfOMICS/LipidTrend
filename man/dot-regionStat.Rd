% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyzeLipidRegion.R
\name{.regionStat}
\alias{.regionStat}
\title{Region Statistics Calculation}
\usage{
.regionStat(X, Y, test = "t.test", ...)
}
\arguments{
\item{X}{Matrix. Lipid abundance data matrix.}

\item{Y}{Matrix. Group information matrix.}

\item{test}{Character. Statistical test to use ("t.test" or "Wilcoxon").
Default is "t.test".}
}
\value{
Numeric vector of statistical test results.
}
\description{
Calculate region statistics using either t-test or Wilcoxon test
}
\keyword{internal}
