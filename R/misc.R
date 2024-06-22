#' Integration by Trapezoidal Rule
#' 
#' Integration by Trapezoidal Rule
#' 
#' Crude Riemann sum approximation.
#' 
#' @param x points of evaluation
#' @param y function values
#' @return A real number.
#' @author R. Koenker
#' @keywords utility
#' @export
traprule <- function(x,y) sum(diff(x) * (y[-1] + y[-length(y)]))/2
#'
#' Smooth a bivariate Kiefer-Wolfowitz NPMLE 
#'
#' @param f bivariate KW fitted object as from GLVmix
#' @param bw bandwidth defaults to bwKW2(f), 
#' @param k kernel 1 for Gaussian, 2 for biweight, 3 for triweight
#' @author R. Koenker
#' @keywords utility
#' @export
#'
KW2smooth <- function(f, bw = NULL, k = 2){
    kernel <- function(x0, x, bw){
	t <- (x - x0)/bw
	switch(k-1,
	(1-t^2)^2*((t> -1 & t<1)-0) * 15/16,
	(1-t^2)^3*((t> -1 & t<1)-0) * 35/32)
    }
    kernel2 = function(u,v,bw)
	kernel(0,u,bw[1]) * kernel(0,v,bw[2])
    fs <- f
    if(max(abs(diff(diff(f$u))), abs(diff(diff(f$v)))) > 1e-10)
        stop("Only equally spaced grids allowed")
    if(!length(bw)) bw = bwKW2(f)
    mu <- length(f$u)
    mv <- length(f$v)
    fsuv <- matrix(0,mu,mv)
    uv <- expand.grid(f$u, f$v)
    for(i in 1:mu){ 
	for(j in 1:mv){
	    fsuv[i,j] <- crossprod(kernel2(f$u[i] - uv[,1], f$v[j] - uv[,2], bw), f$fuv)
	}
    }
    fsuv = c(fsuv)/sum(fsuv)
    fs$fuv <- fsuv
    fs$bw <- bw
    fs
}

#'
#' Smooth a Kiefer-Wolfowitz NPMLE 
#'
#' @param f KW fitted object
#' @param bw bandwidth defaults to 2 * mad
#' @param k kernel 2 for biweight, 3 for triweight
#' @author R. Koenker
#' @keywords utility
#' @export
#'
KWsmooth <- function(f, bw = NULL, k = 2){
    kernel <- function(x0, x, bw){
	t <- (x - x0)/bw
	switch(k, dnorm(t),
	(1-t^2)^k*((t> -1 & t<1)-0) * 15/16,
	(1-t^2)^k*((t> -1 & t<1)-0) * 35/32)
    }
    fs <- f
    if(max(abs(diff(diff(f$x)))) > 1e-10) 
	stop("Only equally spaced grids allowed")
    if(!length(bw)) bw = bwKW(f)
    for(i in 1:length(f$x)) # Would FFT help here?
	fs$y[i] = sum(kernel(f$x[i], f$x, bw)*f$y)
    fs$y <- fs$y/sum(fs$y * diff(f$x)[1])
    fs
}
#' Random sample from KW object
#'
#' @param n sample size
#' @param g KW object
#' @author R. Koenker
#' @keywords utility
#' @export
#'
rKW <- function(n, g){
    sample(g$x, n, prob = g$y, replace = TRUE)
}
#' Quantiles of KW fit
#'
#' @param g KW fitted object
#' @param q vector of quantiles to be computed
#' @author R. Koenker
#' @keywords utility
#' @export
#'
qKW <- function(g, q){
    n <- length(g$y)
    G <- cumsum(c(0, g$y[-n])/sum(g$y))
    g$x[apply(outer(G, q, "<="), 2, sum)]
}
#' Quantiles of bivariate KW fit
#'
#' @param g KW fitted object
#' @param q vector of quantiles to be computed
#' @author R. Koenker
#' @keywords utility
#' @export
#'
qKW2 <- function(g, q){
    fuvm = matrix(g$fuv, length(g$u), length(g$v))
    fum = rowSums(fuvm)
    fum = c(fum)/sum(fum)
    n <- length(fum)
    G = cumsum(c(0, fum[-n]))
    g$u[apply(outer(G,q,"<="),2,sum)]
}
#' Function inversion
#'
#' Given a function, F(x, ...), and a scalar y, find
#' x such that F(x, ...) = y.  Note that there is no
#' checking for the monotonicity of F wrt to x, or that
#' the interval specified is appropriate to the problem.
#' Such fine points are entirely the responsibility of
#' the user/abuser.  If the interval specified doesn't
#' contain a root some automatic attempt to expand the
#' interval will be made.  Originally intended for use 
#' with F as \code{ThreshFDR}.
#'
#' @param y the scalar at which to evaluate the inverse
#' @param F the function 
#' @param interval  the domain within which to begin looking
#' @param ...  other arguments for the function F
#' @author R. Koenker
#' @keywords utility
#' @importFrom stats uniroot
#' @export
#'
Finv = function(y, F, interval = c(0,1), ...) {
    uniroot(function(x) F(x, ...) - y, interval, extendInt = "yes", maxiter = 50)$root 
}
#' Thresholding for False Discovery Rate
#' 
#' This function approximates FDR for various values of lambda
#' and is usually employed in conjunction with \code{Finv} to
#' find an appropriate cutoff value lambda.
#'
#' @param lambda is the proposed threshold
#' @param stat is the statistic used for ranking 
#' @param v is the local false discovery statistic
#' @keywords utility
#' @export
#'
ThreshFDR = function(lambda, stat, v){
    mean((1-v) * (stat > lambda))/mean(stat > lambda)
}   
#'
#' Bandwidth selection for KW smoothing
#'
#' @param g KW fitted object
#' @param k multiplicative fudge factor
#' @param minbw minimum allowed value of bandwidth
#' @author R. Koenker
#' @keywords utility
#' @export
#'
bwKW <- function(g, k = 1, minbw = 0.1){
    m = qKW(g, 0.5)
    max(k * sum(abs(g$x - m) * g$y), minbw) 
}
#' Bandwidth selection for bivariate KW smoothing
#'
#' @param g bivariate KW fitted object
#' @param k multiplicative fudge factor
#' @author R. Koenker
#' @keywords utility
#' @export
#'
bwKW2 <- function(g, k = 1){
    nu = length(g$u)
    nv = length(g$v)
    bw = c(0,0)
    A = matrix(g$fuv, nu, nv)
    fu = list(x = g$u, y = apply(A,1,sum))
    bw[1] = bwKW(fu, k = k)
    fv = list(x = g$v, y = apply(A,2,sum))
    bw[2] = bwKW(fv, k = k)
    bw
}



#' L1norm for piecewise linear functions
#' 
#' Intended to compute the L1norm of the difference between two distribution
#' functions.
#' 
#' Both F and G should be of class \code{stepfun}, and they should be
#' non-defective distribution functions.  There are some tolerance issues in
#' checking whether both functions are proper distribution functions at the
#' extremes of their support.  For simulations it may be prudent to wrap
#' \code{L1norm} in \code{try}.
#' 
#' @param F A stepfunction
#' @param G Another stepfunction
#' @param eps A tolerance parameter
#' @return A real number.
#' @author R. Koenker
#' @keywords utility
#' @export
#' @importFrom graphics plot.default
#' @importFrom stats approxfun is.stepfun knots rnorm runif var
#' @examples
#' 
#' # Make a random step (distribution) function with Gaussian knots
#' rstep <- function(n){
#'         x <- sort(rnorm(n))
#'         y <- runif(n)
#'         y <- c(0,cumsum(y/sum(y)))
#'         stepfun(x,y)
#'         }
#' F <- rstep(20)
#' G <- rstep(10)
#' S <- L1norm(F,G)
#' plot(F,main = paste("||F - G|| = ", round(S,4)))
#' lines(G,col = 2)
#' 
L1norm <- function(F,G, eps = 1e-6) {
        if(!is.stepfun(F) || !is.stepfun(G)) stop("Both F and G must be stepfun")
        xk <- sort(c(knots(F), knots(G)))
        n <- length(xk)
        if(!all.equal(F(xk[1] - eps), G(xk[1] - eps))) stop("F(x[1]-) != G(x[1]-)")
        if(!all.equal(F(xk[n]), G(xk[n]))) stop("F(x[n]) != G(x[n])")
        dy <- (abs(F(xk) - G(xk)))[-n]
        dx <- diff(xk)
        sum(dy * dx)
}
#' Huber density function
#'
#' Huber (1964) least favorable density for the Gaussian contamination model 
#' @param x points to evaluate the density
#' @param mu center of symmetry of the density
#' @param sigma  standard deviation of the nominal Gaussian model
#' @param k Huber k value
#' @param heps Huber epsilson value
#' @return A vector of density values 
#' @importFrom stats dnorm
#' @author R. Koenker
#' @keywords utility
#' @export
#' 
dhuber <- function(x, mu = 0, sigma = 1, k = 1.642, heps = hubereps(k)){
    (x > k)*(1-heps)*dnorm(k)*exp(-k*(x-k)) + 
    (abs(x) <= k)*(1-heps)*dnorm(x) + 
    (x < -k)*(1-heps)*dnorm(k)*exp(k*(x+k)) 
}
#' Huber epsilon 
#'
#' Find the epsilon corresponding to a Huber k value
#' @param k Huber k value
#' @return Huber epsilon value
#' @author R. Koenker
#' @importFrom stats dnorm pnorm
#' @export
#' 
hubereps <- function(k){
    v <- function(z, k)
	2*dnorm(k)/k - 2*pnorm(-k) - z/(1-z)
    uniroot(v, c(0.01, 0.3), k = k)$root
}


