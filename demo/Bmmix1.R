# A test problem to explore Binomial Mixture model with no parametric 
# model for number of trials only marginal data on m

# Model and Data Structure
set.seed(14)
n = 5000
p = sample(c(.1,.5,.9), n, prob = c(.2, .6, .2), replace = TRUE)
m = (p < .5) * rpois(n,10) +  (p == .5) * rpois(n,15) +  (p > .5) * rpois(n,10)
y = rbinom(n,m,p)


Bmmix = function(x, m, v = 50, weights = NULL, ...){
# Binomial mixture with possibly dependent sample sizes
    n <- length(x)
    w <- weights
    if(!length(w))
	w <- rep(1,n)/n
    M <- table(m)
    u <- as.numeric(names(M))
    h <- as.vector(M)/n
    H <- length(h)
    eps <- 1e-04
    if(length(v) == 1)
	v <- seq(eps,1-eps, length = v)
    J <- length(v)
    K <- length(u)
    d <- rep(1,J*K)
    A <- array(NA,c(n, J, K))
    for(k in 1:K)
	A[,,k] <- outer(x,v,function(x,v,uk) dbinom(x, size = uk, prob = v), uk = u[k])  
    m <- J*K
    A <- matrix(A,n,m) 
    B <- kronecker(diag(H), matrix(1,1,J))
    A <- rbind(A, B)

    # Setup Mosek problem a hackish version of KWPrimal
    A <- Matrix::Matrix(A, sparse = TRUE)
    dots <- list(...)
    rtol <- ifelse(length(dots$rtol), dots$rtol, 1e-6)
    verb <- ifelse(length(dots$verb), dots$verb, 0)
    if(length(dots$control)) control <- dots$control
    else control <- NULL
    P <- list(sense = "min")
    A0 <- Matrix::Diagonal(n)
    B0 <- Matrix::Matrix(0, H, n)
    A1 <- Matrix::Matrix(0, n, n)
    B1 <- Matrix::Matrix(0, H, n)
    A0 <- rbind(A0, B0)
    A1 <- rbind(A1, B1)
    P$c <- c(rep(0, m + n), -w)
    P$A <- rbind(cbind(A, -A0, A1),c(d,rep(0,2*n)))
    P$bc <- rbind(c(rep(0, n), h, 1), c(rep(0,n), h, 1))
    P$bx <- rbind(c(rep(0, m + n), rep(-Inf,n)), rep(Inf, m + 2*n))
    P$F <- Matrix::sparseMatrix(c(seq(1,3*n, by = 3), seq(3, 3*n, by = 3)),
		 m + c(1:n, (n+1):(2*n)), x = rep(1,2*n))
    P$g <- rep(c(0,1,0), n)
    P$cones <- matrix(list("PEXP", 3, NULL), nrow = 3, ncol = n)
    rownames(P$cones) <- c("type", "dim", "conepar")
    P$dparam$intpnt_co_tol_rel_gap <- rtol
    z <- Rmosek::mosek(P, opts = list(verbose = verb))
    status <- z$sol$itr$solsta
    if (status != "OPTIMAL")
	warning(paste("Solution status = ", status))
    f <- z$sol$itr$xx[1:m]
    if(min(f) < -rtol)
	warning("mixing distribution has some negative values: consider reducing rtol")
    else f[f < 0] <- 0
    g <- as.vector(A %*% (f * d))[1:n]
    logLik <- n * sum(w * log(g))
    z <- list(v = v, u = u, y = f, g = g, logLik = logLik, status = z$status)
    class(z) <- c("Bmmix", "density")
    return(z)
}
# Now try fitting and visualization
z <- Bmmix(y, m, verb = 5)
g <- expand.grid(v = z$v, u = z$u)
g$y <- z$y
require(lattice)
pl <- cloud(y ~ u * v, data = g, type = "h", lwd = 2, 
        zlim = c(0, max(z$y)), scales = list(arrows = FALSE,
        xlab = expression(u), ylab = expression(v), zlab = "density",
        screen = list(z = 10, x = -70)))
print(pl)





