Pmix <- function(X, v, m = 300, eps = 1e-6, rtol = 1e-06, verb=0, control = NULL){

   # Kiefer-Wolfowitz Estimation of Poisson Mixtures
   # Input:
   #   X is an n vector of observed values
   #   v is a grid of rates on which we evaluate (optional)
   #   g on return is a vector of function values for the Bayes Rule
   # Note that X is reduced to cell counts and rehydrogenated for likelihood computation
   if(missing(v)) v <- seq(max(2*eps, min(X))-eps,max(X)+eps,length = m)
   y <- table(X)
   w <- y/sum(y)
   x <- as.integer(unlist(dimnames(y)))
   d <- diff(v)
   d <- c(d[1],d)
   A <- outer(x,v,"dpois")
   A <- Matrix(A, sparse = TRUE)
   z <- KWDual(x,w,d,A,rtol,verb, control = control)
   f <- z$f
   y <- 0:(max(x) + 1)
   g <- outer(y,v,"dpois") %*% f
   nx <- tabulate(X+1, length(y))
   logLik <- sum(nx * log(g/sum(g)))
z <- list(x = v, y = f, g = g, logLik = logLik, flag = z$status)
class(z) <- "density"
return(z)
}
