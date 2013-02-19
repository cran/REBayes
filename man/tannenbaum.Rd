\name{tannenbaum}
\alias{tannenbaum}
\docType{data}
\title{Gaussian Mixture data }
\description{
Gaussian Location Mixture data to illustrate Mosek tolerance problem
}
\usage{data(tannenbaum)}
\format{
  5000 iid Gaussians
}
\details{
This data set was generated in the course of trying to understand some
anomalies in estimating Gaussian location mixture problems with GLmix.
It is used by the demo tannenbaum to illustrate that sometimes it is
worthwhile reducing the convergence tolerance for mosek.
See \code{demo(tannenbaum)} for further details.
}
\keyword{datasets}
