GVmix <- function(y, id, v, pv = 300, eps = 1e-6, rtol = 1.0e-6, verb=0){

   # Kiefer-Wolfowitz Estimation of Gaussian Variance Mixtures for repeated measures 
   # Input:
   #   y is an N vector of observed values
   #   id is an N vector of indices for the n "individuals"
   #   v is a grid of points on which we evaluate individual variances
   # Output:
   #   v as above
   #   fv mixing density for the variances
   #   flag indicating (non)convergence code (0 is OK)

s <- tapply(y,id,"var")
m <- tapply(y,id,"length")
r <- (m - 1)/2
n <- length(s)
if(missing(v)) v <- seq(min(s) - eps, max(s) + eps, length = pv)
pv <- length(v)
dv <- diff(v)
dv <- c(dv[1],dv)
wv <- rep(1,n)/n

# Note that 2*r*s/theta ~ chisq_2r  so A needs to be an n by p matrix with entries
# f(s,theta) = (r*s/theta)^(r-1) exp(-r*s/theta)/(Gamma(r)*theta/r)

R <- outer(r*s,v,"/")  
vgamma <- outer(gamma(r)/r, v)
r <- outer((m - 1)/2, rep(1,pv))
A <- (exp(-R) * R^(r-1))/vgamma
A <- Matrix(A, sparse = TRUE)
f <- KWDual(s,wv,dv,A, rtol = rtol, verb = verb)
y <- f$f/sum(f$f * dv)
z <- list(x = v, y = y, logLik = f$logLik, flag = f$status)
class(z) <- "density"
return(z)
}
