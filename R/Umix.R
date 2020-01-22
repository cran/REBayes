#' NPMLE for Uniform Scale Mixtures
#' 
#' Kiefer-Wolfowitz Nonparametric MLE for Uniform Scale Mixtures
#' 
#' Kiefer-Wolfowitz MLE for the mixture model \eqn{Y \sim U[0,T], \; T \sim G} 
#' No gridding is required since mass points of the mixing distribution, \eqn{G},
#' must occur at the data points.  This formalism is equivalent, as noted by
#' Groeneboom and Jongbloed (2014) to the Grenander estimator of a monotone
#' density in the sense that the estimated mixture density, i.e. the marginal
#' density of \eqn{Y}, is the Grenander estimate,  see the remark at the end
#' of their Section 2.2.  See also \code{demo(Grenander)}.  Note that this
#' refers to the decreasing version of the Grenander estimator, for the 
#' increasing version try standing on your head.
#' 
#' @param x Data: Sample Observations
#' @param ... other parameters to pass to KWDual to control optimization
#' @return An object of class density with components: 
#'      \item{x}{points of  evaluation on the domain of the density} 
#'      \item{y}{estimated mass at the points x of the mixing density} 
#'      \item{g}{the estimated mixture density function values at x} 
#'      \item{logLik}{Log likelihood value at the proposed solution} 
#'      \item{status}{exit code from the optimizer}
#' @author Jiaying Gu and Roger Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}.  Volume 27, Number 4 (1956), 887-906.
#' 
#' Groeneboom, P. and G. Jongbloed, \emph{Nonparametric Estimation under 
#' Shape Constraints}, 2014, Cambridge U. Press.
#'
#' @keywords nonparametric
#' @export
Umix <- function (x, ...){
    n = length(x)
    v = sort(x)
    A = outer(x, v, "<=")/outer(rep(1,n),v,"*")
    d = rep(1,n)
    w = rep(1,n)/n
    f = KWDual(A,d,w)
    y = f$f
    g = f$g
    logLik = n * sum(w * log(g))
    z = list(x = v, y = y, g = g, logLik = logLik, status = f$status)
    class(z) = c("Umix", "density")
    return(z)
}
