\name{GVmix}
\alias{GVmix}
\title{
GVmix:  Generalized Maximum Likelihood for Empirical Bayes Estimation of 
Gamma Variances
}
\description{
A Kiefer-Wolfowitz procedure for ML estimation of a Gaussian model with
independent variance components with longitudinal data.
}
\usage{
GVmix(y, id, v, pv = 300, eps = 1e-06, rtol = 1.0e-6, verb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ A vector of observations }
  \item{id}{ A strata indicator vector of the same length }
  \item{v}{ A vector of bin boundaries for the variance effects}
  \item{pv}{The number of variance effect bins, if u is missing}
  \item{eps}{A tolerance for determining the support of the bins}
  \item{rtol}{ relative tolerance for dual gap convergence criterion }
  \item{verb}{A flag indicating how verbose the Mosek output should be}
}
\details{
See Gu and Koenker (2012?)
}
\value{
An object of class \code{density}  consisting of the following components:
\item{x}{the variance bin boundaries}
\item{y}{the function values of the mixing density for the variances. }
\item{g}{the function values of the mixture density at the observed y's.}
\item{logLik}{the value of the log likelihood at the solution}
\item{status}{the mosek convergence status.}
}
\references{
Gu Y. and R. Koenker (2013)  
}
\author{R. Koenker}
\keyword{nonparametric}
