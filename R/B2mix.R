#' Bivariate Binomial mixture estimation via Kiefer Wolfowitz MLE
#' 
#' Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of bivariate binomials
#' 
#' This function was inspired by a paper by Kline and Walters (2019) on evaluation of audit
#' experiments for employment discrimination.  An example of its usage is available
#' with `demo(B2mix1)`.  There can be identification issues particularly when the
#' numbers of trials are modest as described in Gu (2020).  Caveat emptor!
#' The predict method for B2mix objects will compute posterior means, 
#'
#' @param x n by 2 matrix of counts of "successes" for binomial observations
#' @param k n by 2 matrix of Number of trials for binomial observations
#' @param u Grid Values for the mixing distribution defaults to equal
#' spacing of length u on [eps, 1- eps], if u is scalar.
#' @param v Grid Values for the mixing distribution defaults to equal
#' spacing of length v on [eps, 1- eps], if v is scalar.
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... Other arguments to be passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{u}{grid of evaluation points of the mixing density} 
#' 	\item{v}{grid of evaluation points of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#' 	\item{g}{estimates of the mixture density at the distinct data values} 
#' 	\item{logLik}{Log Likelihood value at the estimate}
#' 	\item{dy}{Bayes rule estimates of binomial probabilities for distinct data values}
#' 	\item{status}{exit code from the optimizer}
#' @author R. Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. 27, (1956), 887-906.
#'
#' Kline, P. and C. Walters, (2019) Audits as Evidence: Experiments, Ensembles
#' and Enforcement, preprint.
#'
#' Gu, J. (2020) A Collection of Notes on Binomial Mixtures, preprint.
#'
#' @seealso `Bmix` for univariate binomial mixtures.
#' @keywords nonparametric
#' @importFrom stats dbinom
#' @export
B2mix = function (x, k, u = 40, v = 40, weights = NULL, ...) 
{
    n <- NROW(x)
    w <- weights
    if(!length(w))
	w = rep(1,n)/n
    eps <- 1e-04
    if (length(u) == 1) 
        u <- seq(eps, 1 - eps, length = u)
    if (length(v) == 1) 
        v <- seq(eps, 1 - eps, length = v)
    mu <- length(u)
    mv <- length(v)
    uv <- c(outer(u,v))
    d <- rep(1, mu*mv)
    makeA = function(x,k,u,v)
	outer(outer(x, u, function(x,u,k) dbinom(x,k,u), k = k), rep(1,length(v)))
    A1 = makeA(x[,1], k[,1], u, v)
    A2 = makeA(x[,2], k[,2], v, u)
    A = A1 * aperm(A2, c(1,3,2))
    B = NULL
    for(j in 1:mv) B = cbind(B, A[,,j])
    z <- KWDual(B, d, w, ...)
    g <- z$g
    logLik <- n * sum(w * log(g))
    dy <- c((B %*% (z$f * d * uv))/g)
    z <- list(u = u, v = v, y = z$f, g = g, logLik = logLik, dy = dy, status = z$status)
    class(z) <- c("B2mix", "density")
    return(z)
}

