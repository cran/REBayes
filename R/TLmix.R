TLmix <- function(x, v, u, df = 1, m = 300, eps = 1e-6, 
	hist = FALSE, rtol = 1.0e-06, verb=0){

   # Kiefer-Wolfowitz Estimation for Student Location mixtures
   # Input:
   #   x is an n vector of observed values
   #   v is a grid of points on which we evaluate (optional)
   #   u is a grid of points on which we bin the x's
   #   df is the degrees of freedom parameter of the Student base density
   n <- length(x)
   if(missing(v)) v <- seq(min(x)-eps, max(x)+eps, length = m)
   if(hist){
      if(missing(u)) u <- seq(min(x)-eps,max(x)+eps,length = m)
      w <- tabulate(findInterval(x,u))
      x <- (u[-1] + u[-m])/2
      wnz <- (w > 0)
      w <- w[wnz]/sum(w[wnz])
      x <- x[wnz]
      }
   else 
   w <- rep(1,n)/n
   d <- diff(v)
   d <- c(d[1],d)
   A <- dt(outer(x,v,"-"),df = df) 
   A <- Matrix(A, sparse = TRUE)
   f = KWDual(x,w,d,A, rtol = rtol, verb = verb)
z <- list(x = v, y = f$f, g = f$g, logLik = n * f$logLik, flag = f$status)
class(z) <- "density"
return(z)
}
