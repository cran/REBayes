WGLVmix <- function(y, id, w, u, v, pu = 30, pv = 30, eps = 1e-6, rtol = 1.0e-6,
		    verb=0, control = NULL){

# Weighted Kiefer-Wolfowitz Estimation of Gaussian Location and Variance Mixtures
# Produces a bivariate (joint) distribution for the Location and Variance Mixing Distribution
# Input:
#   y is an N vector of observed values
#   id is an N vector of indices for the n "individuals"
#   w is the known variance component for individual i period t
#   u is a grid of points on which we evaluate individual means
#   v is a grid of points on which we evaluate individual variances
#      if v is scalar then it is treated as a fixed variance parameter
# Output:
#   u as above
#   v as above
#   fuv mixing density for the means and variances
#   logLik value
#   flag indicating (non)convergence code

wsum <- tapply(w, id, "sum")
t <- tapply(w * y, id, "sum")/wsum
m <- tapply(y, id, "length")
r <- (m - 1)/2
s <- (tapply(w * y^2, id, "sum") - t^2 * wsum)/(m - 1)
n <- length(s)
if(missing(u)) u <- seq(min(t) - eps, max(t) + eps, length = pu)
if(missing(v)) v <- seq(min(s) - eps, max(s) + eps, length = pv)
if(length(v) == 1) 
    v <- seq(min(s) - eps, max(s) + eps, length = pv)
du <- diff(u)
du <- c(du[1],du)
wu <- rep(1,n)/n
dv <- diff(v)
dv <- c(dv[1],dv)
wv <- rep(1,n)/n

# Note that 2*r*s/theta ~ chisq_2r  so A needs to be an n by p matrix with entries
# f(s,theta) = (r*s/theta)^r exp(-r*s/theta)/(s Gamma(r))

R <- outer(r*s,v,"/")  
sgamma <- outer(s * gamma(r),rep(1,pv))
r <- outer((m - 1)/2, rep(1,pv))
Av <- outer((exp(-R) * R^r)/sgamma, rep(1,pu))
Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum),rep(1,pu)), sqrt(v), "/"))
Au <- Au/outer(outer(1/sqrt(wsum),rep(1,pu)),sqrt(v))
Au <- aperm(Au,c(1,3,2)) # permute Au indices so that they are aligned with those of Av
A <- Av * Au

B <- NULL
for (j in 1:pu)
    B <- cbind(B,A[,,j])
B <- Matrix(B, sparse = TRUE)
duv=kronecker(du,dv)
f <- KWDual(y,wu,duv,B,rtol = rtol, verb = verb)
fuv <- f$f
flag <- f$status
list(u = u, v = v, fuv = fuv, logLik = f$logLik, flag = flag)
}
