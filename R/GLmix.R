GLmix <- function(x, v, m = 300, eps = 1e-6, hist = FALSE, rtol = 1.0e-06, verb=0){

   # Kiefer-Wolfowitz Estimation Gaussian Location mixtures
   # Input:
   #   x is an n vector of observed values
   #   v is a grid of points on which we evaluate (optional)
   #   u is a grid of points on which we bin the x's
   n <- length(x)
   if(missing(v)) v <- seq(min(x)-eps,max(x)+eps,length = m)
   if(hist){
      u <- seq(min(x)-eps,max(x)+eps,length = m)
      w <- tabulate(findInterval(x,u))
      x <- (u[-1] + u[-m])/2
      wnz <- (w > 0)
      w <- w[wnz]/sum(w[wnz])
      x <- x[wnz]
      }
   else 
      w <- rep(1,length(x))/length(x)
   d <- diff(v)
   d <- c(d[1],d)
   A <- dnorm(outer(x,v,"-")) 
   A <- Matrix(A, sparse = TRUE)
   f <- KWDual(x,w,d,A, rtol = rtol, verb = verb)
   y <- f$f
   dy <- as.vector((A %*% (y * v))/(A %*% y))
   o <- order(x)
   g <- approxfun(x[o], w[o]/(sum(f$f)*f$g[o]), rule = 2) 
z <- list(x = v, y = f$f, g = g, dy = dy, logLik = n * f$logLik, flag = f$status)
class(z) <- "density"
return(z)
}
