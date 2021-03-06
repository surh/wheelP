% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_binary.r
\name{plot_binary}
\alias{plot_binary}
\title{Plot results from binary association}
\usage{
plot_binary(Dat, xbreaks = c(1, 2, 5, 10, 20, 50))
}
\arguments{
\item{Dat}{A data.frame. Requires columns StartP and EndP defining
conditions and used for faceting. Column Pi_content indicates the
phosphate concentration. Column bacteria indicates whether bacteria
was added or not.}

\item{xbreaks}{vector of values to be indicated in x-axis.}
}
\value{
A ggplot2 plot object.
}
\description{
Generates density plot of Pi content in four conditions with
and without bacteria
}
\author{
Sur Herrera Paredes
}
\keyword{binP}
\keyword{plots}
