
#' Regularized Logistic Regression 
#' 
#' Logistic Regression with lasso like penalties
#' 
#' In some logistic regression problems, especially those with a large number of fixed effects
#' like the Bradley-Terry rating model, it may be plausible to consider groups of effects that 
#' would be considered equivalence classes.  One way to implement such prior information is to
#' impose some form of regularization penalty.  In the general formulation we are trying to
#' solve the problem:
#' \deqn{ \min \ell (\theta | X, y) + \| D \theta \|_1 }.
#' For example in the Bradley-Terry rating model, we may consider penalties of the form, 
#' \deqn{ \| D \theta \|_1 = \sum_{i < j} |\theta_i - \theta_j | }   
#' so differences in all pairs of ratings are pulled together.  This form of the penalty
#' has been used by Hocking et al (2011) for clustering, by Masarotto and Varin (2012)
#' for estimation of the Bradley Terry model and by Gu and Volgushev (2019) for grouping
#' fixed effects in panel data models.  This is an implementation in
#' Mosek, so the package \pkg{Rmosek} and Mosek must be available at run time.
#' The \code{demo(RLR1)} illustrates use with the conventional lasso penalty and produces a 
#' lasso shrinkage plot.  The \code{demo(RLR2)} illustrates use with the ranking/grouping
#' lasso penalty and produces a plot of how the number of groups is reduced as lambda rises.
#' 
#' @param X a design matrix for the unconstrained logistic regression model
#' @param Y a response vector of Boolean values, or n by 2 matrix of binomials as in \code{glm}
#' @param D is a matrix specifying the penalty, \code{diag(ncol(X))} for the conventional 
#' lasso penalty 
#' @param lambda a scalar specifying the intensity of one's belief in the prior.  No  
#' provision for automatic selection has been made (yet).
#' @param ...  other parameters passed to control optimization:  These may
#' include \code{rtol} the relative tolerance for dual gap convergence criterion,
#' \code{verb} to control verbosity desired from mosek, \code{verb = 0} is quiet,
#' \code{verb = 5} produces a fairly detailed iteration log.  See the documentation for 
#' \code{KWDual} for further details.
#' @return A list with components: \item{coef}{vector of coefficients}\item{logLik}{log likelihood
#' value at the solution}\item{status}{return status from the Mosek optimizer}.
#' @author Roger Koenker with crucial help from Michal Adamaszek of Mosek ApS
#' @references  Gu, J. and Volgushev, S. (2019), `Panel data quantile regression with grouped 
#' fixed effects', \emph{Journal of Econometrics}, 213, 68--91.
#'
#' Hocking, T. D., Joulin, A., Bach, F. and Vert, J.-P. (2011), `Clusterpath: an algorithm for
#' clustering using convex fusion penalties', Proceedings of the 28th International Conference
#' on International Conference on Machine Learning, 745--752.
#'
#' Masarotto, G. and Varin, C. (2012), `The ranking lasso and its application to sport 
#' tournaments', \emph{The Annals of Applied Statistics}, 6, 1949--1970.
#'
#' @keywords regression
#' @export
#'
RLR <- function(X, Y, D, lambda, ...)
{
    P <- list(sense="min")
    n <- nrow(X)
    p <- ncol(X)
    k <- NCOL(Y)
    m <- nrow(D)

    if(NROW(Y) != n) stop("X and Y don't match n")
    if(ncol(D) != p) stop("X and D don't match p")

    dots <- list(...)

    if(utils::packageVersion("Rmosek") < "9")
	stop("RLR requires Mosek Version >= 9")

    # Default mosek method
    rtol <- ifelse(length(dots$rtol), dots$rtol, 1e-6)
    verb <- ifelse(length(dots$verb), dots$verb, 0)
    if(length(dots$control)) control <- dots$control
    else control <- NULL

    if(k == 2){
	# Variables: r(1), theta(d), u(m), t(2n), z1(2n), z2(2n)
	P$c <- c(lambda, rep(0,p+m), Y[,1], Y[,2],rep(0,2*n), rep(0,2*n))
	P$bx <-rbind(rep(-Inf,1+p+m+6*n), rep(Inf,1+p+m+6*n))
    
	# l1 constraints
	A1 <- rbind(cbind(0, D, diag(m), Matrix(0, m, 6*n)),
		    cbind(0, -D, diag(m), Matrix(0, m, 6*n)),
		    c(1,rep(0,p),rep(-1,m),rep(0, 6*n)))
	# z1 + z2 <= 1
	A2 <- sparseMatrix( rep(1:n, 4), 
		    c((1:n)+1+p+m+2*n, (1:n)+1+p+m+3*n,
		    (1:n)+1+p+m+4*n, (1:n)+1+p+m+5*n),
		    x = rep(1, 4*n))
	P$A <- rbind(A1,A2)
	P$bc <- rbind(c(rep(0,1+2*m), rep(-Inf, n)), 
			c(rep(Inf,1+2*m),rep(1, n)))
    
	F0 <- Matrix(nrow=0, ncol = 1+p+m+6*n)
	F1 <- Matrix(nrow=0, ncol = 1+p+m+6*n)
	for(i in 1:n) {
	    xi <- X[i,]
	    F0 <- rbind(F0, sparseMatrix( c(1, 3, 4, rep(6, p), 6),
		    c(1+p+m+2*n+i, 1+p+m+i, 1+p+m+3*n+i, 2:(p+1), 1+p+m+i),
		    x = c(1, -1, 1, -xi, -1), dims = c(6, 1+p+m+6*n) ) )
	    F1 <- rbind(F1, sparseMatrix( c(1, 3, 4, rep(6, p), 6),
		    c(1+p+m+4*n+i, 1+p+m+n+i, 1+p+m+5*n+i, 2:(p+1), 1+p+m+n+i),
		    x = c(1, -1, 1, xi, -1), dims = c(6, 1+p+m+6*n) ) )
	}
	g <- rep(rep(c(0, 1, 0), 4), n)
    
	P$F <- rbind(F0,F1)
	P$g <- g
	P$cones <- matrix(list("PEXP", 3, NULL), nrow=3, ncol=4*n)
	rownames(P$cones) <- c("type","dim","conepar")
    }

    else{

	# Variables: r(1), theta(d), u(m), t(n), z1(n), z2(n)
	P$c <- c(lambda, rep(0,p+m), rep(1, n), rep(0,n), rep(0,n))
	P$bx <-rbind(rep(-Inf,1+p+m+3*n), rep(Inf,1+p+m+3*n))
    
	# l1 constraints
	A1 <- rbind(cbind(0, D, diag(m), Matrix(0, m, 3*n)),
		    cbind(0, -D, diag(m), Matrix(0, m, 3*n)),
		    c(1,rep(0,p),rep(-1,m),rep(0, 3*n)))
	# z1 + z2 <= 1
	A2 <- sparseMatrix( rep(1:n, 2), 
		    c((1:n)+1+p+m+n, (1:n)+1+p+m+2*n),
		    x = rep(1, 2*n))
	P$A <- rbind(A1,A2)
	P$bc <- rbind(c(rep(0,1+2*m), rep(-Inf, n)), 
			c(rep(Inf,1+2*m),rep(1, n)))
    
	# (z1(i), 1, -t(i)) \in \KEXP, 
	# (z2(i), 1, (1-2y(i))*X(i,) - t(i)) \in \KEXP
	FE <- Matrix(nrow=0, ncol = 1+p+m+3*n)
	for(i in 1:n) {
	    yx <- (1 - 2*Y[i])*X[i,] 
	    FE <- rbind(FE, sparseMatrix( c(1, 3, 4, rep(6, p), 6),
		    c(1+p+m+n+i, 1+p+m+i, 1+p+m+2*n+i, 2:(p+1), 1+p+m+i),
		    x = c(1, -1, 1, yx, -1), dims = c(6, 1+p+m+3*n) ) )
	}
	gE <- rep(rep(c(0, 1, 0), 2), n)
    
	P$F <- FE
	P$g <- gE
	P$cones <- matrix(list("PEXP", 3, NULL), nrow=3, ncol=2*n)
	rownames(P$cones) <- c("type","dim","conepar")
    }
    P$dparam$intpnt_co_tol_rel_gap <- rtol

    if(length(control)){
	P$iparam <- control$iparam
	P$dparam <- control$dparam
	P$sparam <- control$sparam
    }
    z <- Rmosek::mosek(P, opts = list(verbose = verb))
    status <- z$sol$itr$solsta
    if (status != "OPTIMAL")
	warning(paste("Solution status = ", status))
    coef <- z$sol$itr$xx[2:(p+1)]
    topt <- z$sol$itr$xx[(1+p+m+1):(1+p+m+k*n)]
    logLik <- -sum(topt)
    z <- list(coef = coef, logLik = logLik, status = status)
    class(z)
    z
}
