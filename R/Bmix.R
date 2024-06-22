#' Binomial mixture estimation via Kiefer Wolfowitz MLE
#' 
#' Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of binomials
#' 
#' The predict method for \code{Bmix} objects will compute means, medians or
#' modes of the posterior according to whether the \code{Loss} argument is 2, 1
#' or 0, or posterior quantiles if \code{Loss} is in (0,1).
#' When the number of trials is small the NPMLE may be non-unique.  This happens
#' when there exists a vector \eqn{v} in the unit simplex of \eqn{R^m}  
#' such that \eqn{Av = f} where \eqn{f = (n_0/n , ... , n_k/n)} the observed frequencies, 
#' and A is the k by m matrix with typical element \deqn{C(k,x) p_j^x (1-p_j)^{k - x}.}
#' If there exists such a solution, it follows that the maximal likelihood value is attained by any Ghat
#' such that \deqn{p_j = \int C(k,j) p^j (1-p)^{k-j} dGhat (p) = n_j/n,} for j = 0, ... , k.
#' There will be many such solutions, but by the Caratheodory theorem any one of them can be expressed
#' as a linear combination of no more than k extreme points of the constraint set.
#' In contrast, when there are no solutions
#' inside the simplex satisfying the equation, then the NPMLE is the unique projection onto the boundary
#' of that set.  To facilitate checking this condition if the \code{check} parameter is \code{TRUE}, the 
#' linear program is feasible and the \code{unique} component is returned as \code{TRUE} if
#' the program is infeasible, and \code{FALSE} is returned otherwise.  This check is restricted to
#' settings in which k is fixed, and \code{collapse} is \code{TRUE}.  See Robbins (1956, p 161) for
#' some further discussion of the binomial mixture model and a very clever alternative approach to
#' prediction.
#'
#' @param x Count of "successes" for binomial observations
#' @param k Number of trials for binomial observations
#' @param v Grid Values for the mixing distribution defaults to equal
#' spacing of length v on [eps, 1- eps], if v is scalar.
#' @param collapse Collapse observations into cell counts.
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param unique option to check unique of reported solution
#' @param ... Other arguments to be passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{x}{grid midpoints of evaluation of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#' 	\item{g}{estimates of the mixture density at the distinct data values} 
#' 	\item{logLik}{Log Likelihood value at the estimate}
#' 	\item{dy}{Bayes rule estimates of binomial probabilities for distinct data values}
#' 	\item{unique}{Flag indicating whether the solution is unique}
#' 	\item{status}{exit code from the optimizer}
#' @author R. Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. 27, (1956), 887-906.
#'
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' \emph{JASA}, 109, 674--685.
#'
#' Robbins, H. (1956) An Empirical Bayes Approach to Statistics, 3rd Berkeley
#' Symposium.
#'
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @keywords nonparametric
#' @importFrom stats dbinom
#' @export
Bmix <- function(x, k, v = 300, collapse = TRUE, weights = NULL, unique = FALSE, ...){

    n <- length(x)
    w <- weights
    if(collapse){ #collapse observations into cell counts
      T <- table(x,k)
      x <- rep(as.numeric(dimnames(T)[[1]]), NCOL(T))
      k <- rep(as.numeric(dimnames(T)[[2]]), each = NROW(T))
      y <- c(T)
      s <- y > 0
      y <- y[s]
      x <- x[s]
      k <- k[s]
      w <- y/sum(y)
      }
   if(!length(w)) w <-  rep(1,n)/n
   eps <- 1e-4
   if(length(v) == 1) v <- seq(eps, 1 - eps, length = v)
   m <- length(v)
   d <- rep(1,m)
   A <- outer(x,v,function(x, v, k) dbinom(x,size = k, prob = v), k = k)
   z <- KWDual(A, d, w, ...)
   g <- z$g
   logLik <- n * sum(w * log(g))
   dy <- as.vector((A%*%(z$f * d * v))/g)
   if(unique){ # Check feasibility LP
       if(var(k) > 0)
	  stop("No uniqueness checking when k is heterogeneous") 
       P <- list(sense = "min")
       P$c <- rep(0, m)
       P$A <- rbind(A,1)
       P$bc <- rbind(blc = c(w,1), buc = c(w,1))
       P$bx <- rbind(blx = rep(0L, m), bux = rep(1L, m))
       result <- Rmosek::mosek(P, ...)
       unique <- grep("INFEAS", result$sol$bas$solsta)>0
   }
   z <- list(x = v, y = z$f, g = g, logLik = logLik, dy = dy, unique = unique, status= z$status)
class(z) <- c("Bmix", "density")
return(z)
}
