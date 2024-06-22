#' Efron Bayesian Deconvolution Estimator for Gaussian Mixtures
#'
#' Efron (2016, 2019) penalized logspline density estimator for Gaussian
#' mixture model g-modeling.  Returns an object of class GLmix to facilitate
#' prediction compatible with Kiefer-Wolfowitz GLmix  estimation.  In particular
#' percentile confidence intervals can be constructed based on posterior quantiles.
#' Assumes homoscedastic standard Gaussian noise, for the moment.
#'
#' @param y Data: Sample Observations
#' @param T Undata: Grid Values defaults equal spacing of with T bins, when T is
#' a scalar
#' @param sigma scale parameter of the Gaussian noise, may take vector value of 
#' length(y)
#' @param df degrees of freedom of the natural spline basis
#' @param c0 penalty parameter for the Euclidean norm penalty.
#' @return An object of class GLmix, density with components: 
#'      \item{x}{points of  evaluation on the domain of the density} 
#'      \item{y}{estimated function values at these points of the mixing density} 
#'      \item{sigma}{returns a sigma = 1 for compatibility with GLmix} 
#' @author Adapted from a similar implementation in the R package deconvolveR of
#' Narasimhan and Efron. 
#' @references Efron, B. (2016) Empirical Bayes deconvolution estimates,
#' Biometrika, 103, 1–20, 
#' Efron, B. (2019) Bayes, Oracle Bayes and Empirical Bayes,  
#' Statistical Science, 34, 177-201.
#'
#' @keywords nonparametric
#' @importFrom stats dnorm nlm
#' @importFrom splines ns
#' @export
#'
BDGLmix <- function(y, T = 300, sigma = 1, df = 5, c0 = 0.1){
    # Bayesian Deconvolution Estimator: Efron (B'ka, 2016)
    eps <- 1e-04
    if(length(T) == 1) T <- seq(min(y)-eps, max(y)+eps, length = T)
    X <- ns(T, df = df)
    a0 <- rep(0, ncol(X))
    A <- dnorm(outer(y,T,"-"), sd = sigma)
    qmle <- function(a, X, A, c0){
        g <- exp(X %*% a)
        g <- g/sum(g)
        f <- A %*% g
        -sum(log(f)) + c0 * sum(a^2)^.5
    }
    ahat <- nlm(qmle, a0, X=X, A=A, c0 = c0)$estimate
    g <- exp(X %*% ahat)
    g <- c(g/sum(g * diff(T)[1]))
    z <- list(x = T,y = g, sigma = sigma)
    class(z) <- c("GLmix", "density")
    z
}
