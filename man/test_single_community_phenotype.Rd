% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_single_community_phenotype.r
\name{test_single_community_phenotype}
\alias{test_single_community_phenotype}
\title{Test the effect of all single communities}
\usage{
test_single_community_phenotype(Dat, dir, var.name, bacteria.col = "Bacteria",
  ref.level = "none", plot = FALSE, f1.extra = "+ Experiment + Plate")
}
\arguments{
\item{Dat}{A data.frame with columns 'StartP', 'EndP', as well as columns
defined by the next parameters}

\item{dir}{Where to store output}

\item{var.name}{Column name containing the phenptype to test}

\item{bacteria.col}{Column name containing the bacterial treatment, this
is the variable of interest.}

\item{ref.level}{Reference level for the bacteria.col factor}

\item{plot}{Logical, produce or not plots?}

\item{f1.extra}{Part of the formula to be included. Normally confounders}
}
\description{
Test the effect of all single communities
}
\author{
Sur Herrera Paredes
}
