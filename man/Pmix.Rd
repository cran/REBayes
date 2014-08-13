\name{Pmix}
\alias{Pmix}
\title{
Poisson mixture estimation via Kiefer Wolfowitz MLE
}
\description{
Poisson mixture estimation via Kiefer Wolfowitz MLE
}
\usage{
Pmix(X, v, m = 300, eps = 1e-06, rtol = 1e-06, verb = 0, control = NULL)
}
\arguments{
  \item{X}{
Data: Sample observations (integer valued)
}
  \item{v}{
Undata:  Grid Values for the mixing distribution defaults to equal spacing
of length m
}
  \item{m}{
Number of default grid points
}
  \item{eps}{
Tolerance parameter to determine support for grid
}
  \item{rtol}{
Tolerance parameter for Mosek converence
}
  \item{verb}{
integer determining how verbose the Mosek output should be
}
\item{control}{ Mosek control list see KWDual documentation}
}
\details{
Kiefer Wolfowitz NPMLE estimation for Poisson mixtures.
}
\value{ An object of class density with components:
\item{x}{points of evaluation of the mixing density} 
\item{y}{function values of the mixing density at x}
\item{g}{function values of the mixture density on 0, 1, ... max(x)+1,
for computing Bayes rule}
\item{logLik}{Log Likelihood value at the estimate}
\item{flag}{exit code from the optimizer} 
}
\references{
Kiefer, J. and J. Wolfowitz
Consistency of the Maximum Likelihood Estimator in the Presence of
Infinitely Many Incidental Parameters \emph{Ann. Math. Statist}.
Volume 27, Number 4 (1956), 887-906.
}
\author{ Roger Koenker }
\keyword{nonparametric}
