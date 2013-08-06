Gompertzmix <- function(x, v, alpha, theta, m=300, weight = rep(1, length(x)),
	eps = 1e-06, hist = FALSE, rtol = 1e-06, verb = 0){
    # Log-scale for frailty Gompertz model
    # f(t|alpha,theta,v) = theta * exp(v) * exp(alpha * t) * 
    #      exp(-(theta/alpha) * exp(v) * (exp(alpha*t)-1))
    # Support is bounded by max and min of MLE for v
    require(reliaR)
    n <- length(x)
    if (missing(v)) # Support is bounded by max and min of MLE for v
        v <- seq(min(-log((theta/alpha)*(exp(alpha*x)-1)))-eps,
	   max(-log((theta/alpha)*(exp(alpha*x)-1)))+eps,length=m) 
    if (hist) {
        u <- seq(min(x) - eps, max(x) + eps, length = m)
        w <- tabulate(findInterval(x, u))
        x <- (u[-1] + u[-m])/2
        wnz <- (w > 0)
        w <- w[wnz]/sum(w[wnz])
        x <- x[wnz]
    }
    else 
        w <- weight
    d <- diff(v)
    d <- c(d[1], d)
    A <- matrix(0, nrow=n,ncol=length(v))
    for (i in 1:n){
        for (j in 1:length(v)){
            A[i,j] <- dgompertz(x[i],alpha=alpha, theta = theta * exp(v[j]))
        }
    }
    A <- Matrix(A, sparse = TRUE)
    f = KWDual(x, w, d, A, rtol = rtol, verb = verb)
    z <- list(x = v, y = f$f, g = f$g, logLik = n * f$logLik, 
    flag = f$status)
    class(z) <- "density"
    return(z)
}
