#' Dual optimization for Kiefer-Wolfowitz problems
#' 
#' Interface function for calls to optimizer from various REBayes functions
#' There is currently only one options for the optimization that based on  Mosek. 
#' It relies on the \pkg{Rmosek} interface to R see installation instructions in
#' the Readme file in the inst directory of this package.  This version of the function
#' is intended to work with versions of Mosek after 7.0.  A more experimental option
#' employing the \pkg{pogs} package available from \url{https://github.com/foges/pogs}
#' and employing an ADMM (Alternating Direction Method of Multipliers) approach has
#' been deprecated, those interested could try installing version 1.4 of REBayes, and
#' following the instructions provided there.
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
#' @return Returns a list with components: \item{f}{dual solution vector, the
#' mixing density} \item{g}{primal solution vector, the mixture density
#' evaluated at the data points} \item{logLik}{log likelihood}
#' \item{status}{return status from Mosek}
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
KWDual <- function(A, d, w, ...){
# Dual Kiefer-Wolfowitz MLE for Mixture Problems
#
#       This version implements a class of density estimators solving:
#
#       min_x  {F(x) := sum -log (x_i)}  s.t. A' x <= d, 0 <= x,  
#
#
#	where e.g.  A = phi(outer(Y,g,"fun")), with Y data and g a grid on the support of Y,
#	and "fun"  is some function representing the dependence of the base distribution.
#
#-------------------------------------------------------------------------------------
#
# Roger Koenker 
#
# First version:24 Feb 2012  
# Revised:	10 Jun 2015 # Simplified signature
# Revised:	 2 Jul 2015 # Added pogs method
# Revised	30 Jan 2019 # Removed pogs method, added Mosek V9 option

n <- nrow(A)
m <- ncol(A)
A <- t(A) 
A <- Matrix::Matrix(A, sparse = TRUE)

dots <- list(...)

# Default mosek method
rtol <- ifelse(length(dots$rtol), dots$rtol, 1e-6)
verb <- ifelse(length(dots$verb), dots$verb, 0)
if(length(dots$control)) control <- dots$control
else control <- NULL

if(utils::packageVersion("Rmosek") < 9){
    P <- list(sense = "min")
    P$c <- rep(0, n)
    P$A <- A
    P$bc <- rbind(rep(0,m),d)
    P$bx <- rbind(rep(0,n),rep(Inf,n))
    opro <- matrix ( list (), nrow =5, ncol = n)
    rownames ( opro ) <- c(" type ","j","f","g","h")
    
    opro[1,] <-  as.list(rep('log',n))
    opro[2,] <-  as.list(1:n)
    opro[3,] <-  as.list(-w)
    opro[4,] <-  as.list(rep(1,n))
    opro[5,] <-  as.list(rep(0,n))
    P$scopt<- list(opro = opro)
    P$dparam$intpnt_nl_tol_rel_gap <- rtol
}
else { #Mosek Version > 9
    P <- list(sense = "min")
    A0 <- Matrix::Matrix(0, m, n)
    P$c <- c(rep(0,n), -w)
    P$A <- cbind(A, A0)
    P$bc <- rbind(rep(0, m), d)
    P$bx <- rbind(c(rep(0, n), rep(-Inf,n)), rep(Inf, 2*n))
    P$F <- sparseMatrix(c(seq(1,3*n, by = 3), seq(3, 3*n, by = 3)),
		 c(1:n, (n+1):(2*n)), x = rep(1,2*n))
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
if(z$response$code != 0)
    stop(paste("Mosek error: ", z$response$msg))
status <- z$sol$itr$solsta
if (status != "OPTIMAL")
        warning(paste("Solution status = ", status))
f <- z$sol$itr$suc
if(min(f) < -rtol) warning("estimated mixing distribution has some negative values: consider reducing rtol")
else f[f < 0] <- 0
g <- as.vector(t(A) %*% (f * d))
list(f = f, g = g, status = status)
}
