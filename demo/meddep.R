require(Matrix)
require(REBayes)
require(Rmosek)
meddep <- function(x, v = 300, mass=1, alpha=1, verb = 0) {
    # This is a primal version of the 1d Renyi estimator
    # Structured after Ivan's ndpccm code
    n <- length(x)
    if(length(v) == 1) v <- seq(min(x), max(x), length = v)
    p <- length(v)
    A <- diff(diff(Diagonal(p)))
    q <- dim(A)[1] # q = p-2
    w <- (c(diff(v),0) + c(0,diff(v)))/2
    E <- apply(Diagonal(p), 2, function(e,v,x) approx(v,e,x)$y, v=v, x=x)
    e <- apply(E,2,mean)
    beta <- alpha/(alpha-1)
    opro <- matrix(list(), nrow=5, ncol=p)
    opro[5,] <- rep(0,p)
    if (alpha == 1) {
	opro[1,] <- "exp"
	opro[2,] <- (1:p)
	opro[3,] <- w
	opro[4,] <- rep(-1,p)
    } 
    else if (alpha == 0) {
	opro[1,] <- "log"
	opro[2,] <- (1:p)
	opro[3,] <- -w
	opro[4,] <- rep(1,p)
    } 
    else {
	opro[1,] <- "pow"
	opro[4,] <- rep(beta,p)
	if (alpha > 1) { 
	    opro[2,] <- ((p+1):(2*p))
	    opro[3,] <- w*rep(1/beta,p)
	    } 
	else { 
	    opro[2,] <- (1:p)
	    opro[3,] <- w*rep(-1/beta,p)
	    }
	}

    P <- list()
    if(alpha > 1) stop("Not implemented")
    P$sense <- "min"
    P$c <- e
    P$A <- A
    P$bc <- rbind(rep(0,q),rep(Inf,q))   
    P$bx <- rbind(rep(0,p),rep(Inf,p))
    P$scopt <- list(opro=opro)

    R <- mosek(P, opts = list(verbose = verb))


    status <- R$sol$itr$solsta
    if(status != "OPTIMAL")
	warning(paste("Solution status = ", status))
    g <- R$sol$itr$xx[1:p]

    if (alpha == 1) f <- exp(-g)
    else if (alpha > 1) f <- pmax(0,g)^(beta-1)
    else f <- abs(g)^(beta-1)
  
    z <- list(x = v, y = f, g = g)
    class(z) <- c("density", "meddep")
    z
}
par(mfrow = c(2,2))
x <- rt(100, 3)
alphas <- c(1,0.5,0,-1)
for(a in alphas){
    zd <- medde(x, lambda = -1, alpha = a)
    zp <- meddep(x, alpha = a)
    plot(zp, main = paste("a = ", a), xlab = "")
    lines(zd, col = 2)
}
