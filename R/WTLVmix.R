WTLVmix <- function(y, id, w, u, v, pu = 300, pv = 300, eps = 1e-6, 
	rtol = 1.0e-6, verb=0, control = NULL)
{

   # Kiefer-Wolfowitz Estimation of Gaussian Location and Variance Mixtures Student t Version
   # Input:
   #   y is an N vector of observed values
   #   id is an N vector of indices for the n "individuals"
   #   w is an N vector of weights for the y observations
   #   u is a grid of points on which we evaluate individual means
   #   v is a grid of points on which we evaluate individual variances
   #       if v is scalar then it is treated as a fixed variance parameter
   # Output:
   #   u as above
   #   v as above
   #   fu mixing density for the means
   #   fv mixing density for the variances
   #   logLik value
   #   flag indicating (non)convergence code

wsum <- tapply(w,id,"sum")
t <- tapply(w*y,id,"sum")/wsum
m <- tapply(y,id,"length")
r <- (m-1)/2
s <- (tapply(w*y^2,id,"sum") - t^2*wsum)/(m-1)
n <- length(s)
if(missing(u)) u <- seq(min(t) - eps, max(t) + eps, length = pu)
if(missing(v)) v <- seq(min(s) - eps, max(s) + eps, length = pv)
if(length(v) == 1) {
	v0 <- v
	v <- seq(min(s) - eps, max(s) + eps, length = pv)
	}
du <- diff(u)
du <- c(du[1],du)
wu <- rep(1,n)/n
wv <- rep(1,n)/n

Au <- dt(outer(t, u, "-") * outer(sqrt(wsum/s),rep(1,pu)),df = m-1)
Au <- Au/outer(sqrt(s/wsum),rep(1,pu))

# Initialize the variances 

if(exists("v0")){ 
	fv <- rep(0,length(v))
	v0 <- findInterval(v0,v)
	fv[v0] <- 1/diff(v)[v0]
	}
else{
	f<-WGVmix(y,id,w, pv = pv, rtol = rtol, verb=verb)
	fv<-f$y
	}

#Now do the mean part. 

Au <- Matrix(Au, sparse = TRUE)
f <- KWDual(t,wu,du,Au, rtol = rtol, verb = verb, control = control)
fu <- f$f
flag <- f$status
list(u = u, fu = fu, v = v, fv = fv, logLik = f$logLik, flag = flag)
}
