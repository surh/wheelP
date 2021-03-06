% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm_binaryP_enrich.r
\name{perm_binaryP_enrich}
\alias{perm_binaryP_enrich}
\title{Permutation enrichment test of phosphate modulators}
\usage{
perm_binaryP_enrich(test, N = 999, qval.thres = 0.1, max.val = Inf,
  min.val = -Inf)
}
\arguments{
\item{test}{result from test_all_strains}

\item{N}{number of permutations}

\item{qval.thres}{q-value threshold}

\item{max.val}{Maximum value for log fold chance}

\item{min.val}{maximum value for log fold change}
}
\value{
A vector of p-values
}
\description{
Tests whether the number of strains that modulate
plant shoot phosphate accumulation is signigicantly
enriched for positive or negative modulators.
}
\seealso{
\code{\link{test_all_strains}}
}
\author{
Sur Herrera Paredes
}
\keyword{binP}
