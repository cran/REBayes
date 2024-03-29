% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLVmix.R
\name{GLVmix}
\alias{GLVmix}
\title{NPMLE of Gaussian Location-Scale Mixture Model}
\usage{
GLVmix(t, s, m, u = 30, v = 30, ...)
}
\arguments{
\item{t}{A vector of location estimates}

\item{s}{A vector of variance estimates}

\item{m}{A vector of sample sizes of the same length as t and s, or if scalar
a common sample size length}

\item{u}{A vector of bin boundaries for the location effects}

\item{v}{A vector of bin boundaries for the variance effects}

\item{...}{optional parameters to be passed to KWDual to control optimization}
}
\value{
A list consisting of the following components:
	\item{u}{midpoints of mean bin boundaries}
	\item{v}{midpoints of variance bin boundaries}
	\item{fuv}{the function values of the mixing density.}
	\item{logLik}{log likelihood value for mean problem}
	\item{du}{Bayes rule estimate of the mixing density means.}
	\item{dv}{Bayes rule estimate of the mixing density variances.}
	\item{A}{Constraint matrix}
	\item{status}{Mosek convergence status}
}
\description{
A Kiefer-Wolfowitz procedure for ML estimation of a Gaussian model with
possibly dependent mean and variance components. This version differs from
\code{WGLVmix} in that it doesn't assume the data is in longitudinal form.
This version assumes a general bivariate distribution for the mixing
distribution. The defaults use a rather coarse bivariate gridding.
}
\references{
Gu, J. and R. Koenker (2014) Heterogeneous Income Dynamics: An
Empirical Bayes Perspective, \emph{JBES},35, 1-16.

Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
\emph{Journal of Statistical Software}, 82, 1--26.
}
\seealso{
WTLVmix for an implementation assuming independent heterogeneity, and WGLVmix
	for a version that requires access to a full longitudinal data structure.
}
\author{
R. Koenker and J. Gu
}
\keyword{nonparametric}
