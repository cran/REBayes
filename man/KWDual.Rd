\name{KWDual}
\alias{KWDual}
\title{
Dual optimization for Kiefer-Wolfowitz problems
}
\description{
Interface function for calls to Mosek from various REBayes functions
}
\usage{
KWDual(x, w, d, A, rtol = 1.0e-6, verb = 0, control = NULL)
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
  \item{control}{
control list consisting of sublists \code{iparam}, \code{dparam}, and \code{sparam},
containing elements of various mosek control parameters.  See the Rmosek and Mosek
manuals for further details.  A prime example is \code{rtol} which should eventually
be deprecated and folded into \code{control}, but will persist for a while for
compatibility reasons.  Another example, the one that motivated the introduction
of \code{control} would be \code{control = list(iparam = list(num_threads = 1))},
which forces Mosek to use a single threaded process.  The default allows Mosek
to uses multiple threads (cores) if available, which is generally desirable, but
may have unintended consequences when running simulations on clusters.
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
