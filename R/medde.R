medde <- function(x, v = 300, lambda = 0.5, alpha = 1, Dorder = 1, 
	rtol = 1e-06, verb = 0){
############################################################################
# fit <- medde(x, v, lambda, alpha, Dorder, rtol, verb)
#
# Maximum Entropy Deregularized Density Estimation (The RMosek Version)
#
#       This version implements a class of penalized density estimators solving:
#
#       min_x phi(x_1)
#
#               s.t. A_1 x_1 - A_2 x_2 <- b, 
#                             0 <<- x_1,  
#                       -lambda <<- x_2 <= lambda
#
#       where x is a vector with two component subvectors:
#
#               x1 is a vector of function values of the density
#               x2 is a vector of dual values
#		lambda is typically positive, if negative then the x_2 constraint
#			is replaced by x_2 <-> 0, which for alpha = Dorder = 1, 
#			constrains the fitted density to be log-concave, for
#			negative lambda fitting for other alphas and Dorders 
#			proceed at your own risk. See also below...
#		phi is an additive convex function in the coordinates of x_1,
#			interpretable as a negative Renyi entropy, see below.
#-------------------------------------------------------------------------------------
#
# This formulation purports to solve a class of dual penalized maximum (Renyi) entropy
# problems in which regularization is achieved by constraining the sup norm of the
# dual vector.  In the primal representation of the problems this corresponds to
# a roughness penalty on the total variation of the Dorder derivative of some
# transformation of the fitted density.  See Renyi.m for examples. 
#
# Input:
#     x      - n observations (aka data)
#     v      - m virtual observations (aka undata), or length m of equispaced grid
#     lambda - total variation penalty parameter, if lambda is in [-1,0], a
#               concavity constraint is imposed on the log density. If lambda 
#               is in (-oo, -1) a convexity constraint on .5 x^2 + log f 
#               is imposed. See Koenker and Mizera (2012) for further 
#               details on this last option, and Koenker and Mizera (2010) 
#               for further details on the concavity constrained options.
#     alpha  - Renyi  entropy parameter
#     Dorder - order of the TV penalty, e.g. 0 for TV(phi(f)), 1 for TV(phi(f)'), ... 
#     rtol   - Mosek convergence tolerance
#     verb   - Mosek verbosity parameter 0 for silent, 5 for verbose
#
# Output:  an object called fit consisting of three components:
#
#	x -- grid values v on which density estimate is evaluated.
#	y -- estimated density values at these points 
#	phi -- n by 2 matrix with original data and log-density at these values
#	T -- Delone triangulation of the evaluation  points (not yet implemented)
#
# Roger Koenker and Ivan Mizera
#
# First version: 22 Dec 2006  (only does Dorder <- 1, and a few alphas)
#                26 Dec 2006  (does Dorder <- 1,2,3, and a few more alphas)
#                28 Dec 2006  (added merge flag) 
#                31 Dec 2006  (log-concave case) 
#                 5 Jan 2007  (bivariate case) 
#                 9 Apr 2008  (well-tempered case) 
#                12 Nov 2010  (RMosek Version:  only 1d) 
############################################################################


if(length(v) == 1){
	eps <- ifelse(lambda < 0, 0.0001,1)
	v <- seq(min(x) - eps, max(x) + eps, length = v)
	}
dimx <- NCOL(x)
dimv <- NCOL(v)
if(dimx != dimv) 
	stop('x and v of different dimensions')
if(dimx == 1)
	mesh <- mesh1(x,v,Dorder)
else if(dimx == 2){
   if(Dorder != 1)
	stop('Dorder must be 1 for bivariate data')
   else {
	stop("Bivariate smoothing not (yet) implemented for REBayes")
	#mesh <- mesh2(x,v,merge)
	#T <- mesh$T
        }
   }
else
   stop('x and v must be either 1d or 2d')


XL <- mesh$XL
XC <- mesh$XC
XJ <- mesh$XJ
XU <- mesh$XU

n <- NROW(x)
p <- NCOL(XJ)
q <- NROW(XJ)
H <- as(XC,"matrix.diag.csr")
C <- rep(0,p+q)
if(lambda > -1)
   A <- cbind(H , t(XJ))
else {
   A <- cbind(H , -t(XJ))
   C[(p+1):(p+q)] <- rep(1,q)
   }
A <- as(as.matrix.csc(A),"dgCMatrix")
L <- c(XL %*% rep(1,n)/n)

if(lambda > 0) { # TV Constraint
   LX <- c(rep(0,p),  -rep(lambda,q))
   UX <- c(rep(Inf,p), rep(lambda,q))
   }
else {  # Concavity/Convexity Constraint
   LX <- rep(0,p+q)
   UX <- rep(Inf,p+q)
   }

P <- list(sense = "min")
P$c <- C
P$A <- A
P$bx <- rbind(LX,UX)
P$bc <- rbind(L,L)

opro <- matrix(list(),nrow = 5, ncol = p)
rownames(opro) <- c(" type ", "j", "f", "g", "h")

#switch solver case{'mosek'} # MOSEK solution  (For another day, only mosek for now)
if(alpha == 1){ # Shannon
    opro[1, ] <- as.list(rep("ent", p))
    opro[2, ] <- as.list(1:p)
    opro[3, ] <- as.list(XC)
    opro[4, ] <- as.list(rep(0, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 0.5){ # Hellinger
    opro[1, ] <- as.list(rep("pow", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(-XC)
    opro[4, ] <- as.list(rep(0.5, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 0){ # Berg
    opro[1, ] <- as.list(rep("log", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(-XC)
    opro[4, ] <- as.list(rep(0, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 2){ # Pearson
    opro[1, ] <- as.list(rep("pow", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(XC)
    opro[4, ] <- as.list(rep(2, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 3){ # Silverman for Good
    opro[1, ] <- as.list(rep("pow", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(XC)
    opro[4, ] <- as.list(rep(3, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else
   stop('specified alpha not (yet) implemented')


P$scopt <- list(opro = opro)
P$dparam$intpnt_nl_tol_rel_gap <- rtol
z <- mosek(P, opts = list(verbose = verb))
status = z$sol$itr$solsta
if(status != "OPTIMAL") warning(paste("Solution status = ", status))
f <- z$sol$itr$xx[1:p]
g <- t(XL) %*% f
o <- order(x)
phi <- cbind(x[o],log(g[o]))
logLik <- sum(phi[,2])
z <- list(x = v, y = f, phi = phi, logLik = logLik, status = status)
class(z) <- "medde"
z
}
plot.medde <- function(x, xlab = "x", ylab = "f(x)", ...){
	plot.default(x, type = "l", xlab = xlab, ylab = ylab, ...)
	}

mesh1 <- function(x, v, Dorder = 1) {
############################################################################
# mesh <- mesh1(x, v, Dorder, merge)
#
# Meshing for univariate MEDDE, old option to merge x and v removed in R version
#
# Input:
#     x      - vector of n observations (aka data)
#     v      - vector of m virtual observations (aka undata)
#     Dorder - order of the TV penalty, e.g. 0 for TV(phi(f)), 1 for TV(phi(f)'), ... 
#
# Output:  
#
#       XL  - the evaluation operator 
#       XC  - the Riemann weights used to integrate 
#       XJ  - the penalty block of the constraint matrix
#       XU  - the unique values at which the density is to be estimated
#
# Roger Koenker, last MeddeR revision 4 Jan 2007
#                     RMosek version 12 Nov 2012
############################################################################

n <- length(x)
p <- length(v)
k <- findInterval(x,v)
h <- diff(v)
ia <- c(1:n,1:n)
ja <- c(k+1,k)
ra <- c((x - v[k])/h[k], (v[k+1] - x)/h[k])
R <- t(new("matrix.coo",ra = ra, ja = as.integer(ja), ia = as.integer(ia), 
	dimension = as.integer(c(n,p))))
s <- 1/h
d <- 0.5*(c(h, 0) + c(0,h))
H <- as(d,"matrix.diag.csr")

if(Dorder %in% 0:2){
  switch(Dorder + 1,
  {#case 0
     q <-  p-1
     IA <- c(1:(p-1), 1:(p-1))
     JA <- c(1:(p-1), 2:p)
     XA <- c(-s[1:(p-1)], s[1:(p-1)])
     },
  {#case 1
     q <-  p-2
     IA <- c(1:(p-2), 1:(p-2), 1:(p-2))
     JA <- c(1:(p-2), 2:(p-1), 3:p)
     XA <- c(-s[1:(p-2)], s[1:(p-2)]+s[2:(p-1)], -s[2:(p-1)])
     D  <- .5 * (h[1:(p-2)]+h[2:(p-1)])
     XA <- XA / c(D, D, D)
     },
  {#case 2
     q <-  p-3
     IA <- c(1:(p-3), 1:(p-3), 1:(p-3), 1:(p-3))
     JA <- c(1:(p-3), 2:(p-2), 3:(p-1), 4:p)
     D1 <-  .5 * (h[1:(p-3)] + h[2:(p-2)])
     D2 <-  .5 * (h[2:(p-2)] + h[3:(p-1)])
     DA <- .5 * (D1 + D2)
     XA <- c(s[1:(p-3)]/D1, - (s[1:(p-3)] +  s[2:(p-2)])/D1 - s[2:(p-2)]/D2,
             s[2:(p-2)]/D1 + s[2:(p-2)]/D2 + s[3:(p-1)]/D2, -s[3:(p-1)]/D2)
     XA <- XA/c(DA, DA, DA, DA)
     }
  )
}
else
   stop("Dorder must be in {0,1,2}")

   XJ <- new("matrix.coo",ra = -XA, ja = as.integer(JA), ia = as.integer(IA), 
	dimension = as.integer(c(q,p)))
   XJ <- as.matrix.csr(XJ)
   XL <- as.matrix.csr(R)
   list(XJ = XJ, XC = d, XL = XL, XU = v)
}

# Functions borrowed from logcondens for computing [rq]medde, etc (slightly modified)
qmedde <- function (p, medde) { 
    if (any(p < 0 | p > 1)) stop("All p's must be in [0, 1]!\n")
    x <- medde$phi[,1]
    phi <- medde$phi[,2]

    Fhat <- function (x, phi) {
       Jexp <- function (x, y, v = 1) {
          m <- length(x)
          z <- exp(x)
          d <- y - x
          k <- (1:m)[abs(d) > 0.005]
          z[k] <- z[k] * (exp(v * d[k]) - 1)/d[k]
          k <- (1:m)[abs(d) <= 0.005]
          z[k] <- z[k] * 
              (v + d[k] * (v/2 + d[k] * (v/6 + d[k] * (v/24 + d[k] * v/120))))
          return(z)
          }
    n <- length(x)
    Fhat <- 1:n * 0
    dx <- diff(x)
    Fhat[2:n] <- cumsum(dx * Jexp(phi[1:(n - 1)], phi[2:n]))
    Fhat[n] <- max(Fhat[n], 1)
    Fhat
    }
    qloglin <- function (u, v) 
        ifelse(abs(v) > 1e-06, log(1 + ((exp(v) - 1) * u))/v, u + v * u * (1 - u)/2)

    Fhat <- Fhat(x,phi)
    n <- length(x)
    m <- length(p)
    qs <- rep(0,m)
    for (i in 1:m) {
        p0 <- p[i]
        if (p0 == 0) {
            q <- -Inf
        }
        if (p0 == 1) {
            q <- x[n]
        }
        if ((p0 > 0) && (p0 < 1)) {
            xj <- max(x[Fhat <= p0])
            j <- length(x[x <= xj])
            u <- (p0 - Fhat[j])/(Fhat[j + 1] - Fhat[j]) 
            v <- (x[j + 1] - x[j]) * (phi[j + 1] - phi[j])
            q <- xj + (x[j + 1] - x[j]) * qloglin(u,v)
        }
        qs[i] <- as.numeric(q)
    }
    qs
}
rmedde <- function(n, medde, smooth = TRUE) {
	z <- qmedde(sort(runif(n)), medde) 
	if(smooth){
	   x <- medde$x
	   fx <- medde$y
	   X <- medde$phi[,1]
	   dx <- diff(medde$x)[1]
	   v <- var(X) - sum(((x - mean(X))^2)*fx)*dx
	   z <- z + rnorm(n,sd = sqrt(v))
	   }
	z
}
