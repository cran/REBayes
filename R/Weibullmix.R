Weibullmix <- function(x, v, alpha, lambda, weight = rep(1, length(x)), 
	m=300, eps = 1e-06, hist = FALSE, rtol = 1e-06, verb = 0, control = NULL){
    # log-scale for frailty Weibull model 
    # f(t|alpha, lambda) = alpha * exp(v) * (lambda * t )^(alpha-1) * 
	# exp(-(lambda * t)^alpha * exp(v))
    # shape = alpha, scale = lambda^(-1) * (exp(v))^(-1/alpha)
    # support bounded by MLE for v = -alpha * log (lambda * t)
    n <- length(x)
    if (missing(v)) 
        v <- seq(min(-alpha*log(lambda * x))-eps,max(-alpha*log(lambda * x))+eps,length=m)
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
    for(j in 1:length(v))
       A[,j] <- dweibull(x,shape=alpha, scale = (1/lambda) * (exp(v[j]))^(-1/alpha))
    A <- Matrix(A, sparse = TRUE)
    f = KWDual(x, w, d, A, rtol = rtol, verb = verb, control = control)
    z <- list(x = v, y = f$f, g = f$g, logLik = n * f$logLik, 
    flag = f$status)
    class(z) <- "density"
    return(z)
}
