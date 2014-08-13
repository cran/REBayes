\name{Gompertzmix}
\alias{Gompertzmix}
\title{
NPMLE for Gompertz Mixtures
}
\description{
Kiefer-Wolfowitz NPMLE for Gompertz Mixtures of scale parameter
}
\usage{
Gompertzmix(x, v, alpha, theta, m = 300, weight = rep(1, length(x)), 
	eps = 1e-06, hist = FALSE, rtol = 1e-06, verb = 0, control = NULL)
}
\arguments{ \item{x}{ Survival times }
  \item{v}{ Grid values for mixing distribution }
  \item{alpha}{ Shape parameter for Gompertz distribution }
  \item{theta}{ Scale parameter for Gompertz Distribution }
  \item{m}{ Number of v grid points }
  \item{weight}{ Frequency (count) for survival times }
  \item{eps}{ tolerance parameter }
  \item{hist}{ If TRUE aggregate to histogram counts }
  \item{rtol}{ convergence tolerance for optimization }
  \item{verb}{ verboseness of optimization }
  \item{control}{ Mosek control parameters see KWDual documentation}
}
\details{
Kiefer Wolfowitz NPMLE density estimation for Gompertz scale mixtures.
The histogram option is intended for relatively large
problems, say n > 1000, where reducing the sample size dimension is desirable.
By default the grid for the binning is equally spaced on the support of the data.
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
}
\author{ Roger Koenker and Jiaying Gu } 

\seealso{ \code{Weibullmix} }

\keyword{nonparametric}
