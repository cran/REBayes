\name{L1norm}
\alias{L1norm}
\title{ L1norm for piecewise linear functions }
\description{
Intended to compute the L1norm of the difference between two distribution functions.
}
\usage{ L1norm(F, G, eps = 1e-06) }
\arguments{
  \item{F}{ A stepfunction }
  \item{G}{ Another stepfunction }
  \item{eps}{ A tolerance parameter }
}
\details{
Both F and G should be of class \code{stepfun}, and they should be non-defective 
distribution functions.  There are some tolerance issues in checking whether both
functions are proper distribution functions at the extremes of their support.  For
simulations it may be prudent to wrap \code{L1norm} in \code{try}.
}
\value{ A real number.  }
\author{ R. Koenker }
\examples{
# Make a random step (distribution) function with Gaussian knots
rstep <- function(n){
        x <- sort(rnorm(n))
        y <- runif(n)
        y <- c(0,cumsum(y/sum(y)))
        stepfun(x,y)
        }
F <- rstep(20)
G <- rstep(10)
S <- L1norm(F,G)
plot(F,main = paste("||F - G|| = ", round(S,4)))
lines(G,col = 2)
}
\keyword{ utility }
