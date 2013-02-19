\name{Bmix}
\alias{Bmix}
\title{ Binomial mixture estimation via Kiefer Wolfowitz MLE }
\description{ Binomial mixture estimation via Kiefer Wolfowitz MLE }
\usage{
Bmix(x, k, v, m = 300, eps = 1e-06, rtol = 1e-09, collapse = TRUE, verb = 0)
}
\arguments{
  \item{x}{ count of "successes" for binomial observations }
  \item{k}{ number of trials for binomial observations }
  \item{v}{ Undata:  Grid Values for the mixing distribution defaults 
	to equal spacing of length m on [eps, 1- eps].  }
  \item{m}{Number of equally spaced grid points desired.}
  \item{eps}{ Tolerance parameter to determine support of grid.  }
  \item{rtol}{ Convergence tolerance for duality gap in Mosek.  }
  \item{collapse}{Collapse observations into cell counts. }
  \item{verb}{ Integer determining verbosity of Mosek output }
}
\details{Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of binomials }

\value{ 
An object of class density with components:
\item{x}{points of evaluation of the mixing density} 
\item{y}{function values of the mixing density at x}
\item{logLik}{Log Likelihood value at the estimate}
\item{flag}{exit code from the optimizer} 
}
\references{
Kiefer, J. and J. Wolfowitz
Consistency of the Maximum Likelihood Estimator in the Presence of
Infinitely Many Incidental Parameters \emph{Ann. Math. Statist}.
Volume 27, Number 4 (1956), 887-906.
}
\author{R. Koenker}
\keyword{nonparametric}
