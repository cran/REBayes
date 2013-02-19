Bmix <- function(x, k, v, m = 300, eps = 1e-6, rtol = 1.0e-9, collapse = TRUE, verb=0){

   # Kiefer-Wolfowitz Estimation of Binomial Mixtures
   # Input:
   #   x is an n vector of observed values from B(k,p)
   #   v is a grid of p in [0,1] on which we evaluate (optional)
   #   y on return is a vector of function values for the mixing density
   #   g on return is a vector of function values for the Bayes Rule
   if(collapse){ #collapse observations into cell counts
      T <- table(x,k)
      x <- rep(as.numeric(dimnames(T)[[1]]), NCOL(T))
      k <- rep(as.numeric(dimnames(T)[[2]]), each = NROW(T))
      y <- c(T)
      s <- y > 0
      y <- y[s]
      x <- x[s]
      k <- k[s]
      w <- y/sum(y)
      }
   else
      w <- rep(1,length(x))/length(x)
   if(missing(v)) v <- seq(eps, 1 - eps, length = m)
   d <- diff(v)
   d <- c(d[1],d)
   A <- outer(x,v,function(x, v, k) dbinom(x,size = k, prob = v), k = k)
   A <- Matrix(A, sparse = TRUE)
   z <- KWDual(x,w,d,A, rtol = rtol, verb = verb)
   f <- z$f
   z <- list(x = v, y = f, logLik = z$logLik, flag = z$status)
class(z) <- "density"
return(z)
}
