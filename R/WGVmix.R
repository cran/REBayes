WGVmix <- function(y, id, w, v, pv = 300, eps = 1e-6, rtol = 1.0e-6, 
		   verb=0, control = NULL){

   # Kiefer-Wolfowitz Estimation of Gaussian Variance Mixtures for repeated measures 
   # Input:
   #   y is an N vector of observed values
   #   w is an N vector of weights
   #   id is an N vector of indices for the n "individuals"
   #   v is a grid of points on which we evaluate individual variances
   # Output:
   #   v as above
   #   fv mixing density for the variances
   #   g mixture density at the observed s's
   #   flag indicating (non)convergence code (0 is OK)

m <- tapply(y,id,"length")
n <- length(m)
wsum <- tapply(w, id, "sum")
s <- (tapply(w * y^2, id, "sum") - tapply(w * y, id, sum)^2/wsum)/(m-1)
if(missing(v)) v <- seq(min(s) - eps, max(s) + eps, length = pv)
pv <- length(v)
dv <- diff(v)
dv <- c(dv[1],dv)
wv <- rep(1,n)/n

# Note that 2*r*s/theta ~ chisq_2r  so A needs to be an n by p matrix with entries
# f(s,theta) = (r*s/theta)^r exp(-r*s/theta)/(s Gamma(r))

r <- (m-1)/2
R <- outer(r*s,v,"/")  
sgamma <- outer(s * gamma(r),rep(1,pv))
r <- outer((m - 1)/2, rep(1,pv))
A <- (exp(-R) * R^r)/sgamma
A <- Matrix(A, sparse = TRUE)
f <- KWDual(s,wv,dv,A,rtol = rtol, verb = verb, control = control)
y <- f$f/sum(f$f * dv)
g <- A %*% y
z <- list(x = v, y = y, g = g, logLik = f$logLik, flag = f$status)
class(z) <- "density"
return(z)
}
