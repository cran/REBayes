#' NPMLE for Gamma Mixtures 
#' 
#' A Kiefer-Wolfowitz MLE for Gamma mixture models 
#' 
#' @param x vector of observed variances
#' @param shape vector of shape parameters corresponding to x
#' @param v A vector of bin boundaries, if scalar then v equally spaced bins
#'	are constructed 
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param eps tolerance for default gridding 
#' @param ...  optional parameters passed to KWDual to control optimization
#' @return An object of class \code{density} with components:
#' 	\item{x}{midpoints of the bin boundaries} 
#' 	\item{y}{estimated function values of the mixing density} 
#' 	\item{g}{function values of the mixture density at the observed x's.} 
#' 	\item{logLik}{the value of the log likelihood at the solution} 
#' 	\item{dy}{Bayes rule estimates of } 
#' 	\item{status}{the Mosek convergence status.}
#' @author J. Gu and R. Koenker
#' @seealso Gammamix for a general implementation for Gamma mixtures
#' @references Gu J. and R. Koenker (2014) Unobserved heterogeneity in 
#' income dynamics: an empirical Bayes perspective, \emph{JBES}, 35, 1-16.
#'
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @keywords nonparametric
#' @importFrom stats dgamma
#' @export
Gammamix <- function(x, v = 300, shape = 1, weights = NULL, eps = 1e-10, ...){

    n = length(x)
    if (length(v) == 1)
	v <- seq(min(x/shape) - eps, max(x/shape) + eps, length = v)
    p <- length(v)
    d <- rep(1,p)
    if(length(weights)) w <- weights
    else w <- rep(1, n)/n
    A <- outer(x, 1/v, FUN = dgamma, shape = shape)
    f <- KWDual(A, d, w, ...)
    y <- f$f
    g <- f$g
    logLik <- n * sum(w * log(g))
    dy <- as.vector((A %*% (y * d * v))/g)
    z <- list(x = v, y = y, g = g, logLik = logLik, dy = dy, status = f$status)
    class(z) <- c("Gammamix", "density")
    return(z)
    }
