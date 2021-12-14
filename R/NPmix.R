#' Normal mixture with Poisson sample size via Kiefer Wolfowitz NPMLE
#' 
#' Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of Normal/Poissons 
#' 
#' The joint distribution of the means and the number of trials determining sample standard
#' deviations is estimated.  The grid specification for means is as for \code{GLmix}
#' whereas the grid for the Poisson rate parameters by default depends on the support of the 
#' observed trials.  There is no predict method as yet.  See \code{demo(NPmix1)}.
#'
#' @param x observed response for Gaussian observations
#' @param m Number of trials for Poisson observations
#' @param v Grid Values for the Gaussian means mixing distribution defaults to equal
#' spacing of length v on [min(x) + eps, max(x) - eps], if v is scalar.
#' @param u Grid Values for the Poisson rate mixing distribution defaults to equal
#' spacing of length u on [min(m) + eps, max(m) - eps], if u is scalar.
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... Other arguments to be passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{v}{grid points of evaluation of the success probabilities} 
#' 	\item{u}{grid points of evaluation of the Poisson rate for number of trials} 
#' 	\item{y}{function values of the mixing density at (v,u)} 
#' 	\item{g}{estimates of the mixture density at the distinct data values} 
#' 	\item{logLik}{Log Likelihood value at the estimate}
#' 	\item{status}{exit code from the optimizer}
#' @author R. Koenker and J. Gu
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. 27, (1956), 887-906.
#'
#' @keywords nonparametric
#' @importFrom stats dnorm
#' @export
NPmix = function(x, m, v = 50, u = 50, weights = NULL, ...){
# Normal mixture with possibly dependent Poisson sample sizes
    n <- length(x)
    w <- weights
    if(!length(w))
	w <- rep(1,n)/n
    #M <- table(m)
    #u <- as.numeric(names(M))
    eps <- 1e-04
    if(length(u) == 1)
	u <- seq(min(m) + eps, max(m) - eps, length = u)
    if(length(v) == 1)
	v <- seq(min(x) + eps,max(x) - eps, length = v)
    J <- length(v)
    K <- length(u)
    d <- rep(1,J*K)
    A <- array(NA,c(n, J, K))
    for(k in 1:K)
	A[,,k] <- outer(x,v,function(x,v,s) dnorm(x - v, sd = 1/sqrt(s)), s = m) * 
	dpois(m, u[k])
    A <- matrix(A,n,J*K)
    z <- KWDual(A, d, w, ...)
    g <- z$g
    logLik <- n * sum(w * log(g))
    z <- list(v = v, u = u, y = z$f, g = g, logLik = logLik, status = z$status)
    class(z) <- c("BPmix", "density")
    return(z)
}
