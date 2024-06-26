% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/B2mix.R
\name{B2mix}
\alias{B2mix}
\title{Bivariate Binomial mixture estimation via Kiefer Wolfowitz MLE}
\usage{
B2mix(x, k, u = 40, v = 40, weights = NULL, ...)
}
\arguments{
\item{x}{n by 2 matrix of counts of "successes" for binomial observations}

\item{k}{n by 2 matrix of Number of trials for binomial observations}

\item{u}{Grid Values for the mixing distribution defaults to equal
spacing of length u on [eps, 1- eps], if u is scalar.}

\item{v}{Grid Values for the mixing distribution defaults to equal
spacing of length v on [eps, 1- eps], if v is scalar.}

\item{weights}{replicate weights for x obervations, should sum to 1}

\item{...}{Other arguments to be passed to KWDual to control optimization}
}
\value{
An object of class density with components: 
	\item{u}{grid of evaluation points of the mixing density} 
	\item{v}{grid of evaluation points of the mixing density} 
	\item{y}{function values of the mixing density at x} 
	\item{g}{estimates of the mixture density at the distinct data values} 
	\item{logLik}{Log Likelihood value at the estimate}
	\item{dy}{Bayes rule estimates of binomial probabilities for distinct data values}
	\item{status}{exit code from the optimizer}
}
\description{
Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of bivariate binomials
}
\details{
This function was inspired by a paper by Kline and Walters (2019) on evaluation of audit
experiments for employment discrimination.  An example of its usage is available
with `demo(B2mix1)`.  There can be identification issues particularly when the
numbers of trials are modest as described in Koenker and Gu (2024).  Caveat emptor!
The predict method for B2mix objects will compute posterior means,
}
\references{
Kiefer, J. and J. Wolfowitz Consistency of the Maximum
Likelihood Estimator in the Presence of Infinitely Many Incidental
Parameters \emph{Ann. Math. Statist}. 27, (1956), 887-906.

Kline, P. and C. Walters, (2019) Audits as Evidence: Experiments, Ensembles
and Enforcement, preprint.

Koenker, R. and Gu, J. (2024) Empirical Bayes: Some Tools, Rules and Duals,
Cambridge University Press.
}
\seealso{
`Bmix` for univariate binomial mixtures.
}
\author{
R. Koenker
}
\keyword{nonparametric}
