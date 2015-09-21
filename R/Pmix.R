#' Poisson mixture estimation via Kiefer Wolfowitz MLE
#' 
#' Poisson mixture estimation via Kiefer Wolfowitz MLE
#' 
#' Kiefer Wolfowitz NPMLE estimation for Poisson mixtures.
#' 
#' @param x Data: Sample observations (integer valued)
#' @param v Grid Values for the mixing distribution defaults to equal
#' spacing of length v when v is specified as a scalar
#' @param ... other parameters passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{x}{points of evaluation of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#' 	\item{g}{function values of the mixture density on \eqn{0, 1, ... max(x)+1}} 
#' 	\item{logLik}{Log Likelihood value at the estimate} 
#' 	\item{dy}{Bayes rule estimate of Poisson rate parameter at each x}  
#'	\item{status}{exit code from the optimizer}
#' @author Roger Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. Volume 27, Number 4 (1956), 887-906.
#' @keywords nonparametric
#' @export
Pmix <- function(x, v = 300, ...){
   n <- length(x)
   eps <- 1e-4
   if(length(v) == 1) 
      v <- seq(max(2*eps, min(x))-eps,max(x)+eps,length = v)
   y <- table(x)
   w <- y/sum(y)
   x <- as.integer(unlist(dimnames(y)))
   m <- length(v)
   d <- diff(v)
   v <- (v[-1] + v[-m])/2
   A <- outer(x,v,"dpois")
   f <- KWDual(A, d, w, ...)
   logLik <- n * sum(w * log(f$g))
   dy <- as.vector((A%*%(f$f * d * v))/f$g)
z <- list(x = v, y = f$f, g = f$g, logLik = logLik, dy = dy, status = f$status)
class(z) <- c("Pmix", "density")
return(z)
}
