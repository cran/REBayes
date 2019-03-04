#' NPMLE for Longitudinal Gaussian Means and Variances Model with Independent Prior
#' 
#' A Kiefer-Wolfowitz NPMLE procedure for estimation of a Gaussian model with
#' independent mean and variance prior components with weighted longitudinal data.
#' This version iterates back and forth from Gamma and Gaussian forms of the likelihood.
#' 
#' @param y A vector of observations
#' @param id A strata indicator vector indicating grouping of y
#' @param w A vector of weights corresponding to y
#' @param u A vector of bin boundaries for the mean effects
#' @param v A vector of bin boundaries for the variance effects
#' @param eps  Convergence tolerance for iterations
#' @param maxit A limit on the number of allowed iterations
#' @param ... optional parameters to be passed to KWDual to control optimization
#' @return A list consisting of the following components: 
#'      \item{u}{midpoints of the mean bin boundaries} 
#'      \item{fu}{the function values of the mixing density of the means } 
#'      \item{v}{midpoints of the variance bin boundaries} 
#'      \item{fv}{the function values of the mixing density of the variances.} 
#'      \item{logLik}{vector of log likelihood values for each iteration} 
#'      \item{du}{Bayes rule estimate of the mixing density means.} 
#'      \item{dv}{Bayes rule estimate of the mixing density variances.} 
#'      \item{status}{Mosek convergence status for each iteration}
#' @author J. Gu and R. Koenker
#' @seealso WGLVmix for a more general bivariate mixing distribution version and
#' 	WTLVmix for an alternative estimator exploiting a Student/Gamma decomposition
#' @references Gu, J. and R. Koenker (2015)  Empirical Bayesball Remixed:  Empirical
#' Bayes Methods for Longitudinal Data, \emph{J. Applied Econometrics}, 32, 575-599.
#'
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @keywords nonparametric
#' @export
WLVmix <- function (y, id, w, u = 300, v = 300, eps = 1e-4, maxit = 2, ...) {
    n <- length(y)
    if (missing(w)) 
        w <- rep(1, n)
    wsum <- tapply(w, id, "sum")
    t <- tapply(w * y, id, "sum")/wsum
    m <- tapply(y, id, "length")
    r <- (m - 1)/2
    s <- (tapply(w * y^2, id, "sum") - t^2 * wsum)/(m - 1)
    n <- length(s)
    logK <- log(gamma(r)) - r * log(r) - 0.5 * log(wsum) - r * 
        log(2 * pi) - log(s^(r - 1)) + 0.5 * tapply(log(w), id, "sum")
    statit <- matrix(0, 2, maxit)
    likit <- rep(0, maxit)
    if (length(u) == 1) 
        u <- seq(min(t) - eps, max(t) + eps, length = u)
    if (length(v) == 1) 
        v <- seq(min(s) - eps, max(s) + eps, length = v)
    pu <- length(u)
    du <- rep(1,pu)
    pu <- length(u)
    wu <- rep(1, n)/n
    FV0 <- GVmix(s, m, v = v, ...)
    statit[1, 1] <- FV0$status
    v <- FV0$x
    pv <- length(v)
    R <- outer(r * s, v, "/")
    G <- outer(s * gamma(r), rep(1, pv))
    r <- outer((m - 1)/2, rep(1, pv))
    Av <- outer((exp(-R) * R^r)/G, rep(1, pu))
    Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum), rep(1, 
        pu)), sqrt(v), "/"))
    Au <- Au/outer(outer(1/sqrt(wsum), rep(1, pu)), sqrt(v))
    Au <- aperm(Au, c(1, 3, 2))
    A <- Av * Au
    B <- matrix(0, n, pu)
    for (i in 1:n) B[i, ] <- t(A[i, , ]) %*% FV0$y/sum(FV0$y)
    FU0 <- KWDual(B, du, wu, ...)
    statit[2, 1] <- FU0$status
    FU0 <- FU0$f
    likit[1] <- sum(log(B %*% FU0/sum(FU0))) + sum(logK)
    FV <- FV0
    FU <- FU0
    it <- 1
    dlik <- Inf
    while ((dlik > eps) && (it < maxit)) {
	it <- it + 1
        C <- matrix(0, n, pv)
        for (i in 1:n) C[i, ] <- (A[i, , ]) %*% as.vector(FU)/sum(FU)
        FV1 <- KWDual(C, dv, wu, ...)
        statit[1, it] <- FV1$status
        FV1 <- FV1$f
        D <- matrix(0, n, pu)
        for (i in 1:n) D[i, ] <- t(A[i, , ]) %*% FV1/sum(FV1)
        FU1 <- KWDual(D, du, wu, ...)
        statit[1, it] <- FU1$status
        FU1 <- FU1$f
        likit[it] <- sum(log(D %*% FU1/sum(FU1))) + sum(logK)
	dlik <- likit[it] - likit[it - 1]
        FV <- FV1
        FU <- FU1
    }
    au <- matrix(0, n, pu)
    for (i in 1:pu) au[, i] <- A[, , i] %*% (dv * FV)
    av <- matrix(0, n, pv)
    for (i in 1:pv) av[, i] <- A[, i, ] %*% (du * FU)
    g <- au %*% (du * FU)
    du <- au %*% (u * du * FU)/g
    dv <- av %*% (v * dv * FV)/g
    list(u = u, fu = FU, v = v, fv = FV, logLik = likit[1:it], du = du, 
        dv = dv, status = statit[,1:it])
}
