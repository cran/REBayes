#' Poisson mixture estimation via Kiefer Wolfowitz MLE
#' 
#' Poisson mixture estimation via Kiefer Wolfowitz MLE
#' 
#' The predict method for \code{Pmix} objects will compute means, medians or
#' modes of the posterior according to whether the \code{Loss} argument is 2, 1
#' or 0, or posterior quantiles if \code{Loss} is in (0,1).
#'
#' @param x Data: Sample observations (integer valued)
#' @param v Grid Values for the mixing distribution defaults to equal
#' spacing of length v when v is specified as a scalar
#' @param support a 2-vector containing the lower and upper support points
#' of sample observations to account for possible truncation.
#' @param exposure observation specific exposures to risk see details
#' @param ... other parameters passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{x}{points of evaluation of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#' 	\item{g}{function values of the mixture density on \eqn{0, 1, ... max(x)+1}} 
#' 	\item{logLik}{Log Likelihood value at the estimate} 
#' 	\item{dy}{Bayes rule estimate of Poisson rate parameter at each x}  
#'	\item{status}{exit code from the optimizer}
#' @details  In the default case \code{exposure = 1} it is assumed that
#' \code{x} contains individual observations that are aggregated into
#' count bins via \code{table}.  When \code{exposure} has the same length as
#' \code{x} then it is presumed to be individual specific risk exposure and
#' the Poisson mixture is taken to be \eqn{x | v ~ Poi(v * exposure)} and the
#' is not aggregated.  See for example the analysis of the Norberg data in
#' Koenker and Gu (2016).
#' @author Roger Koenker and Jiaying Gu
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. Volume 27, Number 4 (1956), 887-906.
#'
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @keywords nonparametric
#' @export
#' @importFrom stats ppois dpois
Pmix <- function (x, v = 300, support = NULL, exposure = NULL, ...) 
{
    n <- length(x)
    eps <- 1e-04
    if (!length(exposure)) 
	exposure <- 1
    if (length(v) == 1) 
        v <- seq(max(2 * eps, min(x/exposure)) - eps, max(x/exposure) + eps, length = v)
    dtpois <- function(x, v)
	dpois(x, v)/(ppois(support[2], v) - ppois(support[1], v))
    m <- length(v)
    d <- rep(1, length(v))
    if (length(exposure) == 1) {
        y <- table(x)
        w <- y/sum(y)
        x <- as.integer(unlist(dimnames(y)))
        if(length(support))
	    A <- outer(x, v, "dtpois")
	else
	    A <- outer(x, v, "dpois")
    }
    else if (length(exposure) == n) {
        w <- rep(1, n)/n
        A <- matrix(0, n, m)
        for (i in 1:n) A[i, ] = sapply(x[i], v * exposure[i], FUN = dpois)
    }
    else
	stop("length(exposure) must be 1 or length(x)")
    f <- KWDual(A, d, w, ...)
    logLik <- n * sum(w * log(f$g))
    dy <- as.vector((A %*% (f$f * d * v))/f$g)
    z <- list(x = v, y = f$f, g = f$g, logLik = logLik, dy = dy, 
        status = f$status)
    class(z) <- c("Pmix", "density")
    return(z)
}
