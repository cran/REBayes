#' Primal optimization for Kiefer-Wolfowitz problems
#' 
#' Interface function for calls to optimizer from various REBayes functions
#' There is currently only one option for the optimization that based on  Mosek. 
#' It relies on the \pkg{Rmosek} interface to R see installation instructions in
#' the Readme file in the inst directory of this package.  This version of the function
#' works only with versions of Mosek 9.0.  This is an experimental alternative to the
#' main KWDual which is the usual interface from fitting functions to Mosek, caveat emptor..
#' 
#' @param A Linear constraint matrix
#' @param d constraint vector
#' @param w weights for \code{x} should sum to one.
#' @param ...  other parameters passed to control optimization:  These may
#' include \code{rtol} the relative tolerance for dual gap convergence criterion,
#' \code{verb} to control verbosity desired from mosek, \code{verb = 0} is quiet,
#' \code{verb = 5} produces a fairly detailed iteration log,
#' \code{control} is a control list consisting of sublists \code{iparam},
#' \code{dparam}, and \code{sparam}, containing elements of various mosek
#' control parameters.  See the Rmosek and Mosek manuals for further details.
#' A prime example is \code{rtol} which should eventually be deprecated and
#' folded into \code{control}, but will persist for a while for compatibility
#' reasons.  The default for \code{rtol} is 1e-6, but in some cases it is
#' desirable to tighten this, say to 1e-10.  Another example that motivated the introduction of
#' \code{control} would be \code{control = list(iparam = list(num_threads =
#' 1))}, which forces Mosek to use a single threaded process.  The default
#' allows Mosek to uses multiple threads (cores) if available, which is
#' generally desirable, but may have unintended (undesirable) consequences when running
#' simulations on clusters.
#' @return Returns a list with components: \item{f}{primal solution vector, the
#' mixing density} \item{g}{the mixture density
#' evaluated at the data points} \item{logLik}{log likelihood}
#' \item{status}{return status from Mosek}.  Mosek termination messages are
#' treated as warnings from an R perspective since solutions producing, for example,
#' MSK_RES_TRM_STALL: The optimizer is terminated due to slow progress, may still
#' provide a satisfactory solution, especially when the return status variable is
#' "optimal".
#' @author R. Koenker
#' @references
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' \emph{JASA}, 109, 674--685.
#'
#' Mosek Aps (2015) Users Guide to the R-to-Mosek Optimization Interface, 
#' \url{https://docs.mosek.com/8.1/rmosek/index.html}.  
#' 
#' Koenker, R. and J. Gu, (2017) REBayes: An {R} Package for Empirical Bayes Mixture Methods,
#' \emph{Journal of Statistical Software}, 82, 1--26.
#' @keywords nonparametrics
#' @importFrom methods as new
#' @export
KWPrimal <- function(A, d, w, ...){
#
# First version:  15 May 20 20 

n <- nrow(A)
m <- ncol(A)
A <- Matrix::Matrix(A, sparse = TRUE)

dots <- list(...)

# Default mosek method
rtol <- ifelse(length(dots$rtol), dots$rtol, 1e-6)
verb <- ifelse(length(dots$verb), dots$verb, 0)
if(length(dots$control)) control <- dots$control
else control <- NULL

if(utils::packageVersion("Rmosek") < "9")
    stop("Version Mosek 9 only")
else { #Mosek Version => 9
    P <- list(sense = "min")
    A0 <- Matrix::Diagonal(n)
    A1 <- Matrix::Matrix(0, n, n)
    P$c <- c(rep(0, m + n), -w)
    P$A <- rbind(cbind(A, -A0, A1),c(d,rep(0,2*n)))
    P$bc <- rbind(c(rep(0, n),1), c(rep(0,n), 1))
    P$bx <- rbind(c(rep(0, m + n), rep(-Inf,n)), rep(Inf, m + 2*n))
    P$F <- Matrix::sparseMatrix(c(seq(1,3*n, by = 3), seq(3, 3*n, by = 3)),
		 m + c(1:n, (n+1):(2*n)), x = rep(1,2*n))
    P$g <- rep(c(0,1,0), n)
    P$cones <- matrix(list("PEXP", 3, NULL), nrow = 3, ncol = n)
    rownames(P$cones) <- c("type", "dim", "conepar")
    P$dparam$intpnt_co_tol_rel_gap <- rtol
}
if(length(control)){
    P$iparam <- control$iparam
    P$dparam <- control$dparam
    P$sparam <- control$sparam
}
z <- Rmosek::mosek(P, opts = list(verbose = verb))
status <- z$sol$itr$solsta
if (status != "OPTIMAL")
    warning(paste("Solution status = ", status))
f <- z$sol$itr$xx[1:m]
if(min(f) < -rtol) 
    warning("estimated mixing distribution has some negative values: consider reducing rtol")
else f[f < 0] <- 0
g <- as.vector(A %*% (f * d))
list(f = f, g = g, status = status)
}
