\name{TLVmix}
\alias{TLVmix}
\title{
TLVmix:  Maximum Likelihood for Empirical Bayes Estimation of 
Gaussian Means and Variances 
}
\description{
A Kiefer-Wolfowitz procedure for ML estimation of a Gaussian model with
independent mean and variance components with longitudinal data.
This version employs Student t likelihood formulation
}
\usage{
TLVmix(y, id, u, v, pu = 300, pv = 300, eps = 1e-06, rtol = 1.0e-6, verb = 0)
}
\arguments{
  \item{y}{ A vector of observations }
  \item{id}{ A strata indicator vector of the same length }
  \item{u}{ A vector of bin boundaries for the mean effects}
  \item{v}{ A vector of bin boundaries for the variance effects}
  \item{pu}{The number of mean effect bins, if u is missing}
  \item{pv}{The number of variance effect bins, if u is missing}
  \item{eps}{A tolerance for determining the support of the bins}
  \item{rtol}{ relative tolerance for dual gap convergence criterion }
  \item{verb}{A flag indicating how verbose the Mosek output should be}
}
\details{
The vector \code{v} can be specified as a scalar; in this case no mixture
distribution for variances is estimated, instead the variance distribution
is set to be point mass one at the value specified by the scalar.  This is
(potentially) useful for inference about the need for a variance mixture.
See Gu and Koenker (2013)
}
\value{
A list consisting of the following components:
\item{u}{the mean bin boundaries}
\item{fu}{the function values of the mixing density for the means }
\item{v}{the variance bin boundaries}
\item{fv}{the function values of the mixing density for the variances. }
\item{logLik}{log likelihood value for mean problem}
\item{status}{Mosek convergence status}

}
\references{
Gu, J. and R. Koenker (2013)  Empirical Bayesball 
}
\author{R. Koenker}
\keyword{nonparametric}
