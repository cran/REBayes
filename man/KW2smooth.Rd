% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{KW2smooth}
\alias{KW2smooth}
\title{Smooth a bivariate Kiefer-Wolfowitz NPMLE}
\usage{
KW2smooth(f, bw = NULL, k = 2)
}
\arguments{
\item{f}{bivariate KW fitted object as from GLVmix}

\item{bw}{bandwidth defaults to bwKW2(f),}

\item{k}{kernel 1 for Gaussian, 2 for biweight, 3 for triweight}
}
\description{
Smooth a bivariate Kiefer-Wolfowitz NPMLE
}
\author{
R. Koenker
}
\keyword{utility}
