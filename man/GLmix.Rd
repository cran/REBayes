\name{GLmix}
\alias{GLmix}
\title{MLE (Kiefer-Wolfowitz)  Density Estimation}
\description{Density estimation based on Kiefer Wolfowitz method}
\usage{ GLmix(x, v, m = 300, eps = 1e-06, hist = FALSE, rtol = 1.0e-6, verb = 0) }
\arguments{
  \item{x}{Data:  Sample Observations} 
  \item{v}{Undata:  Grid Values defaults equal spacing of length m } 
  \item{m}{Number of (default) grid points}
  \item{eps}{tolerance parameter}
  \item{hist}{If TRUE then aggregate x to histogram weights} 
  \item{rtol}{ relative tolerance for dual gap convergence criterion }
  \item{verb}{integer determining how verbose the output should be}
}
\details{
Kiefer Wolfowitz MLE density estimation as proposed by Jiang and Zhang for the
Gaussian compound decision problem.  The histogram option is intended for large
problems, say n > 1000, where reducing the sample size dimension is desirable.
By default the grid for the binning is equally spaced on the support of the data.
This function does the normal convolution problem, for gamma mixtures of variances
see \code{GVmix}, or for mixtures of both means and variances \code{GLVmix}. 
The optimization is carried out by calls to the function \code{mosek} in the
required package \pkg{Rmosek}.
}
\value{An object of class density with components
\item{x}{points of evaluation on the domain of the density}
\item{y}{estimated function values at the points x, the mixing density}
\item{logLik}{Log likelihood value at the proposed solution}
\item{flag}{exit code from the optimizer}
}
\references{
Kiefer, J. and J. Wolfowitz
Consistency of the Maximum Likelihood Estimator in the Presence of 
Infinitely Many Incidental Parameters \emph{Ann. Math. Statist}. 
Volume 27, Number 4 (1956), 887-906.

Jiang, Wenhua and Cun-Hui Zhang
General maximum likelihood empirical Bayes estimation of normal means
\emph{Ann. Statist.}, Volume 37, Number 4 (2009), 1647-1684.
}

\author{Roger Koenker}
\keyword{nonparametric}
