\name{KWDual}
\alias{KWDual}
\title{
Dual optimization for Kiefer-Wolfowitz problems
}
\description{
Interface function for calls to Mosek from various REBayesfunctions
}
\usage{
KWDual(x, w, d, A, rtol = 1.0e-6, verb = 0)
}
\arguments{
  \item{x}{
Observed data
}
  \item{w}{
weights for \code{x}
}
  \item{d}{
constraint vector
}
  \item{A}{
Linear constraint matrix
}
  \item{rtol}{
relative tolerance for dual gap convergence criterion
}
  \item{verb}{
verbosity desired from mosek
}
}
\details{
See userguide for Rmosek for further details.
}
\value{
Returns a list with components:
\item{f}{dual solution vector, the mixing density}
\item{g}{primal solution vector, the mixture density evaluated at the data points}
\item{logLik}{log likelihood}
\item{status}{return status from Mosek}
}
\author{
R. Koenker
}
\keyword{ nonparametrics }
