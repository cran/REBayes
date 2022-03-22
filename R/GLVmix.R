#' NPMLE of Gaussian Location-Scale Mixture Model
#'
#' A Kiefer-Wolfowitz procedure for ML estimation of a Gaussian model with
#' possibly dependent mean and variance components. This version differs from
#' \code{WGLVmix} in that it doesn't assume the data is in longitudinal form.
#' This version assumes a general bivariate distribution for the mixing
#' distribution. The defaults use a rather coarse bivariate gridding.
#'
#' @param t A vector of location estimates
#' @param s A vector of variance estimates
#' @param m A vector of sample sizes of the same length as t and s, or if scalar
#' 	a common sample size length
#' @param u A vector of bin boundaries for the location effects
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
#' @seealso WTLVmix for an implementation assuming independent heterogeneity, and WGLVmix
#' 	for a version that requires access to a full longitudinal data structure.
#' @keywords nonparametric
#' @export
GLVmix <- function (t, s,  m, u = 30, v = 30, ...) 
{
    n <- length(t)
    w <- rep(1, n)/n
    eps <- 1e-04
    if(length(m) == 1) m <- rep(m, length(t))
    r <- (m - 1)/2
    if (length(u) == 1) 
        u <- seq(min(t) - eps, max(t) + eps, length = u)
    if (length(v) == 1) 
        v <- seq(min(s) - eps, max(s) + eps, length = v)
    v <- v[v > 0]
    pu <- length(u)
    du <- rep(1, pu)
    pv <- length(v)
    dv <- rep(1, pv)
    Av <- matrix(NA, n, pv)
    for (i in 1:n){
    	for (j in 1:pv){
    		Av[i,j] = dgamma(s[i], r[i], scale = v[j]/r[i])
	}
    }
    Av <- outer(Av, rep(1, pu))
    Av <- aperm(Av,c(1,3,2))
    Au <- dnorm(outer(outer(t, u, "-") * 
	      outer(sqrt(m), rep(1, pu)), sqrt(v), "/"))
    Au <- Au/outer(outer(1/sqrt(m), rep(1, pu)), sqrt(v))
    Auv <- Av * Au
    A <- NULL
    for (j in 1:pv) A <- cbind(A, Auv[, , j])
    duv = as.vector(kronecker(du, dv))
    f <- KWDual(A, duv, w, ...)
    fuv <- f$f
    uv <- expand.grid(alpha = u, theta = v)
    g <- as.vector(A %*% (duv * fuv))
    logLik <- n * sum(w * log(f$g))
    du <- A%*%(uv[,1] * duv * fuv)/g  #Bayes rule for u: E(u|t,s)
    dv <- A %*% (uv[,2] * duv * fuv)/g  # Bayes rule for v: E(v|t,s)
    z <- list(u = u, v = v, fuv = fuv, logLik = logLik, 
	      du = du, dv = dv, A = A, status = f$status)
    class(z) <- "GLVmix"
    z
}
