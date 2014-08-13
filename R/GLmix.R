GLmix <- function (x, v, m = 300, sigma = 1, eps = 1e-06, hist = FALSE, 
    rtol = 1e-06, verb = 0, control = NULL) 
{
    n <- length(x)
    if (missing(v)) 
        v <- seq(min(x) - eps, max(x) + eps, length = m)
    if (hist) {
      histbin <- function(x, m = 300, eps = 1e-06) {
        u <- seq(min(x) - eps, max(x) + eps, length = m)
        w <- tabulate(findInterval(x, u))
        x <- (u[-1] + u[-m])/2
        wnz <- (w > 0)
        w <- w[wnz]/sum(w[wnz])
        list(x = x[wnz], w = w)
        }
      if(length(sigma) == 1){
	  h <- histbin(x, m)
	  x <- h$x
	  w <- h$w
	}
      else { # create sigma bins for histogram
	    sus <- sort(unique(sigma))
	    us <- match(sigma, sus)
	    nus <- table(us)
	    if(min(nus) < 100) stop("too few obs in some sigma bin")
	    h <- as.list(1:length(nus))
	    for(i in 1:length(sus))
		h[[i]] <- histbin(x[us == i],m)
	    x <- unlist(lapply(h, function(f) f$x))
	    w <- unlist(lapply(h, function(f) f$w))
	    w <- w/sum(w)
	    sigma <- rep(sus,unlist(lapply(h, function(f) length(f$x))))
      }
    }
    else w <- rep(1, length(x))/length(x)
    d <- diff(v)
    d <- c(d[1], d)
    A <- dnorm(outer(x, v, "-"), sd = sigma)
    A <- Matrix(A, sparse = TRUE)
    f <- KWDual(x, w, d, A, rtol = rtol, verb = verb, control = control)
    y <- f$f
    dy <- as.vector((A %*% (y * v))/(A %*% y))
    o <- order(x)
    g <- approxfun(x[o], w[o]/(sum(f$f) * f$g[o]), rule = 2)
    z <- list(x = v, y = f$f, g = g, sigma = sigma, dy = dy, 
        logLik = n * f$logLik, flag = f$status)
    class(z) <- c("GLmix", "density")
    return(z)
}
