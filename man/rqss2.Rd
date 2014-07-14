\name{rqss2}
\alias{rqss2}
\title{
Quantile Regression L2 Penalized Smoothing Spline 
}
\description{
    Estimation of nonparametric conditional quantile functions using the classical (Wahba)
    L2 roughness penalty.
}
\usage{
rqss2(x, y, tau = 0.5, lambda = 1, dual = TRUE, rtol = 1e-06, verb = 0)
}
\arguments{
  \item{x}{
      the design points
}
  \item{y}{
      the response values
}
  \item{tau}{
      the quantile of interest
}
  \item{lambda}{
      the penalty parameter
}
  \item{dual}{
      the computational method, dual if TRUE, primal otherwise 
}
  \item{rtol}{
      tolerance for the solution
}
  \item{verb}{
      verbosity of the output, use verb = 5 for details of the solution iterations
}
}
\details{
    Calls the Mosek primal-dual interior point solver to compute solutions to the quadratic
    program defining the spline estimator.  
}
\value{
    Returns the fitted values at the (ordered) design points.
}
\references{
Henrik Alsing Friberg (2012). Rmosek: The R-to-MOSEK Optimization
  Interface. R package version 1.2.4.
  http://CRAN.R-project.org/package=Rmosek
}
\author{
    Roger Koenker
}
\note{
    Only does bivariate scatterplot smoothing, see \code{rqss} in the quantreg package for 
    more general additive modeling features.
}

\seealso{
\code{\link[quantreg]{rqss}}.
}
\examples{
# Test Problem (as usual) from Ruppert, Wand, and Carroll

n <- 50
x <- sort(runif(n))
g0 <- function(x, j = 4)
    sqrt(x *(1-x)) * sin(2*pi*(1 + 2^((9-4*j)/5))/(x + 2^((9-4*j)/5)))
dgp <- function(x,f=rnorm,j=2,sigma=.2,gamma=0)
        f(length(x))*sigma*(1+gamma*x)
y <- g0(x) + dgp(x)
plot(x,y, cex = .5, col = "grey")
lines(x, g0(x))
fd <- rqss2(x,y, lambda = 0.0001)
if(fd$status == "OPTIMAL") lines(x, fd$fit, col = 2)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{smooth}
