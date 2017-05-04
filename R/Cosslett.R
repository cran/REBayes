#' Kiefer-Wolfowitz estimator for Cosslett (1983) estimator
#'
#' Kiefer-Wolfowitz-Cosslett  estimator for binary response model.  
#'
#' In the primal form of the problem the pseudo log likelihood is:
#' 
#' 	\deqn{l(f|y) =  sum_i [ y_i \log \sum_j (I(v_j <= x_i) * f_j) + 
#' 		(1 - y_i) \log \sum_j (I(v_j > x_i) * f_j) ]}
#'
#' as usual the implementation used here solves the corresponding dual problem.
#' Cumsum of the output y gives the CDF of the unobserved utility difference.
#' See the \code{demo(Cosslett1)}and  \code{demo(Cosslett2)} for illustrations
#' without any covariate, and \code{demo(Cosslett3)} for an illustration with a
#' covariate using profile likelihood.  This model is also known as current
#' status linear regression in the biostatistics literature, see e.g. Groeneboom
#' and Hendrickx (2016) for recent results and references.
#'
#' @param y is the binary outcome
#' @param x is the observed utility difference between two choices, it would be
#' possible to extend this to make x a linear (index) function of some parameters
#' @param v the unobserved utility difference taking values on a grid, by default
#' this grid is equally spaced with 300 distinct points, however it is known that
#' the mass points for the problem are located at the data points, x, so users may
#' wish to set \code{v = sort(x)} although if the sample size is large this can be
#' slow.
#' @param weights  replicate weights for x observations, should sum to 1 
#' @param ... optional parameters to be passed to KWDual to control optimization
#' @return an object of class density with the components:
#' 	\item{x}{points of evaluation of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#' 	\item{logL}{log likelihood of estimated model} 
#'	\item{status}{exit code from the optimizer}
#' @author Jiaying Gu and Roger Koenker
#' @references Kiefer, J. and J. Wolfowitz (1956) Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters, \emph{Ann. Math. Statist}, 27, 887-906.
#'
#' Cosslett, S. (1983) Distribution Free Maximum Likelihood Estimator of the 
#' Binary Choice Model, \emph{Econometrica}, 51, 765-782.
#'
#' Groeneboom, P. and K. Hendrickx (2016) Current Status Linear Regression,
#' preprint available from \url{https://arxiv.org/abs/1601.00202}.
#' @keywords nonparametric
#' @export

Cosslett <- function(x, y, v = 300, weights = NULL, ...) {
    n <- length(x)
    eps <- 1e-4
    if (length(v) == 1) 
        v <- seq(min(x)-eps, max(x) + eps, length = v)
    if(length(weights)) w <- weights
    else w <- rep(1, n)/n
    d <- rep(1, length(v))
    A <- outer(x, v, ">=")
    A <- (y == 1) * A + (y == 0) * (1 - A)
    f <- KWDual(A, d, w, ...)
    logL <- sum(log(A %*% f$f))
    z <- list(x = v, y = f$f, logL = logL, status = f$status)
    class(z) <- c("Cosslett", "density")
    return(z)
}
