#' Binomial mixtures with Poisson Trials via Kiefer Wolfowitz NPMLE
#' 
#' Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of Poisson Binomials
#' 
#' The joint distribution of the probabilities of success and the number of trials
#' is estimated.  The grid specification for success probabilities is as for \code{Bmix}
#' whereas the grid for the Poisson rate parameters is currently the support of the 
#' observed trials.  There is no predict method as yet.  See \code{demo(BPmix1)}.
#'
#' @param x Count of "successes" for binomial observations
#' @param m Number of trials for binomial observations
#' @param v Grid Values for the mixing distribution defaults to equal
#' spacing of length v on [eps, 1- eps], if v is scalar.
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... Other arguments to be passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{v}{grid points of evaluation of the success probabilities} 
#' 	\item{u}{grid points of evaluation of the Poisson rate for number of trials} 
#' 	\item{y}{function values of the mixing density at (v,u)} 
#' 	\item{g}{estimates of the mixture density at the distinct data values} 
#' 	\item{logLik}{Log Likelihood value at the estimate}
#' 	\item{status}{exit code from the optimizer}
#' @author R. Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. 27, (1956), 887-906.
#'
#' @keywords nonparametric
#' @importFrom stats dbinom
#' @export
BPmix <-  function(x, m, v = 50, weights = NULL, ...){
# Binomial mixture with possibly dependent Poisson sample sizes
    n <- length(x)
    w <- weights
    if(!length(w))
	w <- rep(1,n)/n
    M <- table(m)
    u <- as.numeric(names(M))
    eps <- 1e-04
    if(length(v) == 1)
	v <- seq(eps,1-eps, length = v)
    J <- length(v)
    K <- length(u)
    d <- rep(1,J*K)
    A <- array(NA,c(n, J, K))
    for(k in 1:K)
	A[,,k] <- outer(x,v,function(x,v,k) dbinom(x, size = k, prob = v), k = m) * 
	    dpois(m, u[k])
    A <- matrix(A,n,J*K)
    z <- KWDual(A, d, w, ...)
    g <- z$g
    logLik <- n * sum(w * log(g))
    z <- list(v = v, u = u, y = z$f, g = g, logLik = logLik, status = z$status)
    class(z) <- c("BPmix", "density")
    return(z)
}
