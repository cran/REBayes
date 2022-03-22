#' Weighted NPMLE of Longitudinal Gaussian Mean and Variances Model
#' 
#' A Kiefer-Wolfowitz procedure for ML estimation of a Gaussian model with
#' dependent mean and variance components and weighted longitudinal data.
#' This version assumes a general bivariate distribution for the mixing
#' distribution. The defaults use a rather coarse bivariate gridding.
#' In contrast to the function \code{GLVmix} the full longitudinal data
#' structure is required for this function and the likelihood evaluation
#' reflects this difference.
#' 
#' @param y A vector of observations
#' @param id A strata indicator vector of the same length
#' @param w A vector of weights
#' @param u A vector of bin boundaries for the mean effects
#' @param v A vector of bin boundaries for the variance effects
#' @param ... optional parameters to be passed to KWDual to control optimization
#' @return A list consisting of the following components: 
#' 	\item{u}{midpoints of mean bin boundaries} 
#' 	\item{v}{midpoints of variance bin boundaries} 
#' 	\item{fuv}{the function values of the mixing density.} 
#' 	\item{logLik}{log likelihood value for mean problem} 
#' 	\item{du}{Bayes rule estimate of the mixing density means.} 
#' 	\item{dv}{Bayes rule estimate of the mixing density variances.} 
#' 	\item{A}{Constraint matrix} 
#' 	\item{status}{Mosek convergence status}
#' @author R. Koenker and J. Gu
#' @references Gu, J. and R. Koenker (2014) Heterogeneous Income Dynamics: An
#' Empirical Bayes Perspective, \emph{JBES},35, 1-16. 
#'
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @seealso WTLVmix for an implementation assuming independent heterogeneity,
#' GLVmix for an implementation that assumes the availability of only the summary 
#' statistics but not the full longitudinal data structure.
#' @keywords nonparametric
#' @export
WGLVmix <- function(y, id, w, u = 30, v = 30, ...){

    n <- length(y)
    eps <- 1e-4
    if(missing(w)) w <- rep(1,n)
    wsum <- tapply(w, id, "sum")
    t <- tapply(w * y, id, "sum")/wsum
    m <- tapply(y, id, "length")
    r <- (m - 1)/2
    s <- (tapply(w * y^2, id, "sum") - t^2 * wsum)/(m - 1)
    n <- length(s)
    if(length(u) == 1) u <- seq(min(t) - eps, max(t) + eps, length = u)
    if(length(v) == 1) v <- seq(min(s) - eps, max(s) + eps, length = v)
    v <- v[v > 0]
    pu <- length(u)
    du <- rep(1,pu)
    wu <- rep(1,n)/n
    pv <- length(v)
    dv <- rep(1,pv)
    Av <- matrix(NA, n, pv)
    for(i in 1:n){
	for(j in 1:pv){
	    Av[i,j] <- dgamma(s[i], r[i], scale = v[j]/r[i])
	}
    }
    Av <- outer(Av, rep(1,pu))
    Av <- aperm(Av, c(1,3,2))
    Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum), rep(1, pu)), sqrt(v), "/"))
    Auv <- Av * Au
    A <- NULL
    for(j in 1:pv)
	A <- cbind(A, Auv[,,j])
    duv = as.vector(kronecker(du, dv))
    uv <- expand.grid(alpha = u, theta = v)
    f <- KWDual(A, duv, wu, ...)
    fuv <- f$f
    g = f$g
    status <- f$status
    r <- (m-1)/2
    logK <- log(gamma(r)) - r * log(r) - 0.5 * log(wsum) - 
	r * log(2*pi) - log(s^(r-1)) + 0.5 * tapply(log(w), id, "sum")
    logLik <- sum(log(g)) + sum(logK)
    du <- A%*%(uv[,1] * duv * fuv)/g  #Bayes rule for u: E(u|t,s)
    dv <- A %*% (uv[,2] * duv * fuv)/g  # Bayes rule for v: E(v|t,s)
    z <- list(u = u, v = v, fuv = fuv, logLik = logLik, du = du, dv = dv,
	A = A, status = status)
    class(z) <- c("WGLVmix", "GLVmix")
    z
}
