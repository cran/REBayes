#' Maximum Entropy [De]Regularized Density Estimation
#' 
#' Density estimation based on maximum entropy methods
#' 
#' See the references for further details. And also Mosek "Manuals". The
#' acronym, according to the urban dictionary has a nice connection to
#' a term used in Bahamian dialect, mostly on the Family Islands like Eleuthera
#' and Cat Island meaning "mess with" "get involved," "get entangled," "fool
#' around," "bother:"
#' "I don't like to medder up with all kinda people"
#' "Don't medder with people (chirren)"
#' "Why you think she medderin up in their business."
#' 
#' This version implements a class of penalized density estimators solving:
#
#' \deqn{\min_x \phi(x_1) | A_1 x_1 - A_2 x_2 = b,  0 \le x_1, -\lambda \le x_2 \le \lambda }{
#' min_x \phi(x_1) | A_1 x_1 - A_2 x_2 = b,  0 \le x_1, -\lambda \le x_2 \le \lambda }
#
#' where \eqn{x} is a vector with two component subvectors: \eqn{x_1} is a 
#' vector of function values of the density \eqn{x_2} is a vector of dual values,
#' \eqn{\lambda} is typically positive, and controls the fluctuation of the Dorder
#' derivative of some transform of the density. When alpha = 1 this transform is
#' simply the logarithm of the density, and Dorder = 1 yields a piecewise exponential
#' estimate; when Dorder = 2 we obtain a variant of Silverman's (1982) estimator
#' that shrinks the fitted density toward the Gaussian, i.e. with total variation
#' of the second derivative of \eqn{log f} equal to zero.  See demo(Silverman) for
#' an illustration of this case.  If \eqn{\lambda} is in \eqn{(-1,0]} then the 
#' \eqn{x_2} TV constraint is replaced by \eqn{x_2 \geq 0}, which for \eqn{\alpha = 1}, 
#' constrains the fitted density to be log-concave; for \eqn{\alpha = 0.5},  \eqn{-1/\sqrt f}
#' is constrained to be concave; and for \eqn{\alpha \le 0}, \eqn{1/f^{\alpha -1}} is
#' constrained to be concave.  In these cases no further regularization of the smoothness
#' of density is required as the concavity constraint acts as  regularizer.
#' As explained further in Koenker and Mizera (2010) and
#' Han and Wellner (2016) decreasing \eqn{\alpha} constrains the fitted density to lie 
#' in a larger class of quasi-concave
#' densities.  See \code{demo(velo)} for an illustration of these options, but be aware
#' that more extreme \eqn{\alpha} pose more challenges from an numerical optimization
#' perspective.  Fitting for \eqn{\alpha < 1} employs a fidelity criterion closely 
#' related to Renyi entropy that is more suitable than likelihood for very peaked, or very heavy
#' tailed target densities.  For \eqn{\lambda < 0}  fitting for \code{Dorder != 1}
#' proceed at your own risk.  A closely related problem is illustrated in the demo
#' Brown which imposes a convexity constraint 
#' on \eqn{0.5 x^2 + log f(x)}. This ensures that the resulting Bayes rule,
#' aka Tweedie formula, is monotone in \eqn{x}, as described further in Koenker and
#' Mizera (2013).  
#'

#' @param x Data: either univariate or bivariate, the latter is highly experimental 
#' @param v Undata: either univariate or bivariate, univariate default is an
#' equally spaced grid of 300 values, for bivariate data there is not (yet) a default.
#' @param lambda total variation penalty parameter, if lambda is in [-1,0], a
#' concavity constraint is imposed. see Koenker and Mizera (2010) for
#' further details on the concavity constrained options.  
#' @param alpha Renyi entropy parameter characterizing fidelity criterion
#' by default 1 is log-concave and 0.5 is Hellinger.
#' @param Dorder Order of the derivative operator for the penalty
#' default is Dorder = 1, corresponding to TV norm constraint on the first derivative,
#' or a concavity constraint on some transform of the density.
#' @param w weights associated with x,
#' @param mass  normalizing constant for fitted density,
#' @param rtol Convergence tolerance for Mosek algorithm,
#' @param verb Parameter controlling verbosity of solution, 0 for silent, 5
#' gives rather detailed iteration log.
#' @param control Mosek control list see KWDual documentation
#' @return An object of class "medde" with components \item{x}{points of
#' evaluation on the domain of the density} \item{y}{estimated function values
#' at the evaluation points x}  \item{status}{exit status from Mosek}
#' @author Roger Koenker and Ivan Mizera
#' @seealso This function is based on an earlier function of the same name in
#' the deprecated package MeddeR that was based on an R-Matlab interface.
#' A plotting method is available, or medde estimates can be added to plots
#' with the usual \code{lines(meddefit, ...} invocation.  For log concave
#' estimates there is also a quantile function \code{qmedde} and a random
#' number generation function \code{rmedde}, eventually there should be
#' corresponding functionality for other alphas.
#' @references  Chen, Y. and R.J. Samworth, (2013) "Smoothed log-concave
#' maximum likelihood estimation with applications", \emph{Statistica Sinica},
#' 23, 1373--1398.
#' 
#' Han, Qiyang and Jon Wellner (2016) ``Approximation and estimation of s-concave 
#' densities via Renyi divergences, \emph{Annals of Statistics}, 44, 1332-1359.
#'
#' Koenker, R and I. Mizera, (2007) ``Density Estimation by Total Variation
#' Regularization,'' \emph{Advances in Statistical Modeling and Inference:
#' Essays in Honor of Kjell Doksum}, V.N. Nair (ed.), 613-634.
#' 
#' Koenker, R and I. Mizera, (2006) ``The alter egos of the regularized maximum
#' likelihood density estimators: deregularized maximum-entropy, Shannon,
#' Renyi, Simpson, Gini, and stretched strings,'' \emph{ Proceedings of the 7th
#' Prague Symposium on Asymptotic Statistics}.
#' 
#' Koenker, R and I. Mizera, (2010) ``Quasi-Concave Density Estimation''
#' \emph{Annals of Statistics}, 38, 2998-3027.
#' 
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' JASA, 109, 674--685.
#' 
#' Koenker, R and I. Mizera, (2014) ``Convex Optimization in R.'',
#' \emph{Journal of Statistical Software}, 60, 1-23.
#' @keywords nonparametric
#' @export
#' @import Matrix
#' @importFrom stats dpois
#' @examples
#' 
#' #Maximum Likelihood Estimation of a Log-Concave Density
#' set.seed(1968)
#' x <- rgamma(50,10)
#' m <- medde(x, v = 50, lambda = -.5, verb = 5)
#' plot(m, type = "l", xlab = "x", ylab = "f(x)")
#' lines(m$x,dgamma(m$x,10),col = 2)
#' title("Log-concave Constraint")
#' 
#' #Maximum Likelihood Estimation of a Gamma Density with TV constraint
#' set.seed(1968)
#' x <- rgamma(50,5)
#' f <- medde(x, v = 50, lambda = 0.2, verb = 5)
#' plot(f, type = "l", xlab = "x", ylab = "f(x)")
#' lines(f$x,dgamma(f$x,5),col = 2)
#' legend(10,.15,c("ghat","true"),lty = 1, col = 1:2)
#' title("Total Variation Norm Constraint")
#' 
#' 

medde <- function (x, v = 300, lambda = 0.5, alpha = 1, Dorder = 1, 
    w = NULL, mass = 1, rtol = 1e-06, verb = 0, control = NULL) 
{
    n <- length(x)
    if (length(v) == 1) 
        v <- seq(min(x), max(x), length = v)
    p <- length(v)
    d <- (c(diff(v),0) + c(0,diff(v)))/2
    A <- switch(Dorder + 1,
	t(diff(Diagonal(length(v)))),
	t(diff(diff(Diagonal(length(v))))),
	t(diff(diff(diff(Diagonal(length(v)))))))
    q <- dim(A)[2]
    if(length(w)) 
	e <- c(w %*% apply(diag(p),2,function(e, v, x) 
		     approx(v, e, x)$y, v = v, x = x))
    else
	e <- apply(apply(diag(p),2,function(e, v, x) 
		     approx(v, e, x)$y, v = v, x = x),2,mean)
    if(is.finite(mass)) dv <- mass * d/sum(d) else dv <- d
    beta <- alpha/(alpha - 1)
    if (lambda > 0) {
        LX <- c(rep(0, p), -rep(lambda, q))
        UX <- c(rep(Inf, p), rep(lambda, q))
    }
    else {
        LX <- rep(0, p + q)
        UX <- rep(Inf, p + q)
    }
    opro <- matrix(list(), nrow = 5, ncol = p)
    opro[2,] <- 1:p
    opro[5,] <- rep(0,p)
    if (alpha == 1) {
        opro[1, ] <- "ent"
        opro[3, ] <- -dv
        opro[4, ] <- rep(0, p)
    }
    else if (alpha == 0) {
        opro[1, ] <- "log"
        opro[3, ] <- dv
        opro[4, ] <- as.list(rep(1, p))
    }
    else {
        opro[1, ] <- "pow"
        opro[3, ] <- -sign(beta) * dv
        opro[4, ] <- rep(alpha, p)
    }
    P <- list(sense = "max")
    P$c <- rep(0, p + q)
    P$A <- cBind(Diagonal(p, x = dv), A)
    P$bx <- rbind(LX, UX)
    P$bc <- rbind(e, e)
    P$scopt <- list(opro = opro)
    P$dparam$intpnt_nl_tol_rel_gap <- rtol
    if (length(control)) {
        P$iparam <- control$iparam
        P$dparam <- control$dparam
        P$sparam <- control$sparam
    }
    z <- Rmosek::mosek(P, opts = list(verbose = verb))
    if(z$response$code != 0)
	stop(paste("Mosek error: ", z$response$msg))
    status = z$sol$itr$solsta
    if (status != "OPTIMAL") 
        warning(paste("Solution status = ", status))
    f <- z$sol$itr$xx[1:p]
    if(is.finite(mass)) f <- mass * f/sum(d) 
    #if (alpha == 1) 
    #    logLik <- sum(log(g))
    #else logLik <- NULL
    z <- list(x = v, y = f, status = status)
    class(z) <- "medde"
    z
}


#' Plotting method for medde objects
#'
#' @param x object obtained from medde fitting
#' @param ... other parameters to be passed to plot method
#' @importFrom stats  approx
#' @importFrom graphics contour
#' @importFrom graphics plot
plot.medde <- function(x, ...){
    plot.default(x$x, x$y, type = "l", xlab = "x", ylab = "f(x)", ...)
}


#' Quantile function for medde estimate
#'
#' Slightly modified version  borrowed from the package logcondens 
#' Todo:  extend this to cases with \eqn{\alpha != 1}.
#' @param p vector of probabilities at which to evaluate the quantiles
#' @param medde fitted object from medde
#' @keywords nonparametric
#' @export
qmedde <- function (p, medde) { 
    if (any(p < 0 | p > 1)) stop("All p's must be in [0, 1]!\n")
    x <- medde$phi[,1]
    phi <- medde$phi[,2]

    Fhat <- function (x, phi) {
       Jexp <- function (x, y, v = 1) {
          m <- length(x)
          z <- exp(x)
          d <- y - x
          k <- (1:m)[abs(d) > 0.005]
          z[k] <- z[k] * (exp(v * d[k]) - 1)/d[k]
          k <- (1:m)[abs(d) <= 0.005]
          z[k] <- z[k] * 
              (v + d[k] * (v/2 + d[k] * (v/6 + d[k] * (v/24 + d[k] * v/120))))
          return(z)
          }
    n <- length(x)
    Fhat <- 1:n * 0
    dx <- diff(x)
    Fhat[2:n] <- cumsum(dx * Jexp(phi[1:(n - 1)], phi[2:n]))
    Fhat[n] <- max(Fhat[n], 1)
    Fhat
    }
    qloglin <- function (u, v) 
        ifelse(abs(v) > 1e-06, log(1 + ((exp(v) - 1) * u))/v, u + v * u * (1 - u)/2)

    Fhat <- Fhat(x,phi)
    n <- length(x)
    m <- length(p)
    qs <- rep(0,m)
    for (i in 1:m) {
        p0 <- p[i]
        if (p0 == 0) {
            q <- -Inf
        }
        if (p0 == 1) {
            q <- x[n]
        }
        if ((p0 > 0) && (p0 < 1)) {
            xj <- max(x[Fhat <= p0])
            j <- length(x[x <= xj])
            u <- (p0 - Fhat[j])/(Fhat[j + 1] - Fhat[j]) 
            v <- (x[j + 1] - x[j]) * (phi[j + 1] - phi[j])
            q <- xj + (x[j + 1] - x[j]) * qloglin(u,v)
        }
        qs[i] <- as.numeric(q)
    }
    qs
}

#' Random number generation from a medde estimate
#'
#' @param n number of observations desired in calls to rmedde
#' @param medde fitted medde object for calls in qmedde and rmedde
#' @param smooth option to draw random meddes from the smoothed density
#' @keywords nonparametric
#' @export
rmedde <- function(n, medde, smooth = TRUE) {
	z <- qmedde(sort(runif(n)), medde) 
	if(smooth){
	   x <- medde$x
	   fx <- medde$y
	   X <- medde$phi[,1]
	   dx <- diff(medde$x)[1]
	   v <- var(X) - sum(((x - mean(X))^2)*fx)*dx
	   z <- z + rnorm(n,sd = sqrt(v))
	   }
	z
}
