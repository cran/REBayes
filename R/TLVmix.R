TLVmix <- function(y, id, u, v, pu = 300, pv = 300, eps = 1e-6, rtol = 1.0e-6,verb=0){

   # Kiefer-Wolfowitz Estimation of Gaussian Location and Variance Mixtures
   # (This differs from GLVmix in that the likelihood formulation is Student t)
   # Input:
   #   y is an N vector of observed values
   #   id is an N vector of indices for the n "individuals"
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

s <- tapply(y,id,"var")
t <- tapply(y,id,"mean")
m <- tapply(y,id,"length")
r <- (m-1)/2
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

Au <- dt(outer(t, u, "-") * outer(sqrt(m/s),rep(1,pu)),df = m-1)
Au <- Au/outer(sqrt(s/m),rep(1,pu))

# Initialize the variances 

if(exists("v0")){ 
	fv <- rep(0,length(v))
	v0 <- findInterval(v0,v)
	fv[v0] <- 1/diff(v)[v0]
	}
else{
	f<-GVmix(y,id,rtol = rtol, pv = pv, verb=verb)
	fv<-f$y
	}

#Now do the mean part. 

Au <- Matrix(Au, sparse = TRUE)
f <- KWDual(t,wu,du,Au,rtol = rtol, verb = verb)
fu <- f$f
flag <- f$status
list(u = u, fu = fu, v = v, fv = fv, logLik = f$logLik, flag = flag)
}
