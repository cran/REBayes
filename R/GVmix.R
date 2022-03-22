#' NPMLE for Gaussian Variance Heterogeneity
#' 
#' A Kiefer-Wolfowitz MLE for Gaussian models with independent variances.  This
#' can be viewed as a general form for \eqn{\chi^2} mixtures, see \code{Gammamix}
#' for a more general form for Gamma mixtures.
#' 
#' @param x vector of observed variances
#' @param m vector of sample sizes corresponding to x
#' @param v A vector of bin boundaries, if scalar then v equally spaced bins
#'	are constructed 
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ...  optional parameters passed to KWDual to control optimization
#' @return An object of class \code{density} with components:
#' 	\item{x}{midpoints of the bin boundaries} 
#' 	\item{y}{estimated function values of the mixing density} 
#' 	\item{g}{function values of the mixture density at the observed x's.} 
#' 	\item{logLik}{the value of the log likelihood at the solution} 
#' 	\item{dy}{Bayes rule estimates of } 
#' 	\item{status}{the Mosek convergence status.}
#' @author R. Koenker
#' @seealso Gammamix for a general implementation for Gamma mixtures
#' @references 
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' \emph{JASA}, 109, 674--685.
#'
#' Gu J. and R. Koenker (2014) Unobserved heterogeneity in 
#' income dynamics: an empirical Bayes perspective, \emph{JBES}, 35, 1-16. 
#'
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @keywords nonparametric
#' @export
GVmix <- function(x, m, v = 300, weights = NULL, ...){
    n = length(x)
    r <- (m - 1)/2
    eps <- 1e-8
    if (length(v) == 1)
        v <- seq(min(x) - eps, max(x) + eps, length = v)
    v <- v[v > 0]
    p <- length(v)
    d <- rep(1,p)
    p <- length(v)
    if(length(weights)) w <- weights
    else w <- rep(1, n)/n
    R <- outer(r * x, v, "/")
    A <- matrix(NA, n, p)
    for(i in 1:n){
	for(j in 1:p){
	    A[i,j] <- dgamma(x[i], r[i], scale = v[j]/r[i])
	}
    }
    f <- KWDual(A, d, w, ...)
    y <- f$f
    g <- f$g
    logLik <- n * sum(w * log(g))
    dy <- as.vector((A %*% (y * d * v))/g)
    z <- list(x = v, y = y, g = g, logLik = logLik, dy = dy, status = f$status)
    class(z) <- c("GVmix", "density")
    return(z)
    }
