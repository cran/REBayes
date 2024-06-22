Brown <- function(x, v = 300, rtol = 1e-06, verb = 0, control = NULL){
    # Monotone Bayes Rule Density Estimation for Gaussian Mixture Model
    if(length(v) == 1){
	eps <- ifelse(lambda < 0, 0.0001,1)
	v <- seq(min(x) - eps, max(x) + eps, length = v)
	}
    n <- length(x)
    p <- length(v)
    k <- findInterval(x,v)
    h <- diff(v)
    ia <- c(1:n,1:n)
    ja <- c(k+1,k)
    ra <- c((x - v[k])/h[k], (v[k+1] - x)/h[k])
    XL <- t(sparseMatrix(ia, ja, x = ra, dims = c(n,p)))
    s <- 1/h
    XC <- 0.5*(c(h, 0) + c(0,h))
    q <-  p-2
    IA <- c(1:(p-2), 1:(p-2), 1:(p-2))
    JA <- c(1:(p-2), 2:(p-1), 3:p)
    XA <- c(-s[1:(p-2)], s[1:(p-2)]+s[2:(p-1)], -s[2:(p-1)])
    D  <- .5 * (h[1:(p-2)]+h[2:(p-1)])
    XA <- XA / c(D, D, D)
    XJ <- sparseMatrix(IA, JA, x = -XA, dims = c(q,p))
    H <- Diagonal(x = XC)
    C <- rep(0,p+q)
    A <- cbind(H , -t(XJ))
    C[(p+1):(p+q)] <- rep(1,q)
    w <- rep(1,n)/n
    L <- as.vector(XL %*% w)
    LX <- rep(0,p+q)
    UX <- rep(Inf,p+q)

    P <- list(sense = "min")
    if(utils::packageVersion("Rmosek") < "9") {
	P$c <- C
	P$A <- A
	P$bx <- rbind(LX,UX)
	P$bc <- rbind(L,L)
	if(length(control)){
	    P$iparam <- control$iparam
	    P$dparam <- control$dparam
	    P$sparam <- control$sparam
	}
	opro <- matrix(list(),nrow = 5, ncol = p)
	rownames(opro) <- c(" type ", "j", "f", "g", "h")
	opro[1, ] <- as.list(rep("ent", p))
	opro[2, ] <- as.list(1:p)
	opro[3, ] <- as.list(XC)
	opro[4, ] <- as.list(rep(0, p))
	opro[5, ] <- as.list(rep(0, p))
	P$scopt <- list(opro = opro)
	P$dparam$intpnt_nl_tol_rel_gap <- rtol
    }
    else{
        P$c <- c(C, rep(-1,p))
        P$A <- cbind(A, Matrix(0,p,p))
	P$bx <- rbind(c(LX, rep(-Inf, p)),c(UX, rep(Inf, p)))
	P$bc <- rbind(L,L)
	P$F <- sparseMatrix(c(seq(3, 3 * p, by = 3), seq(1, 3 * p, by = 3)), 
	    c(1:p, (p + q + 1):(2 * p + q)), x = rep(-1, 2 * p))
        P$g <- rep(c(0, 1, 0), p)
        P$cones <- matrix(list("PEXP", 3, NULL), nrow = 3, ncol = p)
        rownames(P$cones) <- c("type", "dim", "conepar")
        P$dparam$intpnt_co_tol_rel_gap <- rtol
    }
    if (length(control)) {
        P$iparam <- control$iparam
        P$dparam <- control$dparam
        P$sparam <- control$sparam
    }
z <- Rmosek::mosek(P, opts = list(verbose = verb))
status = z$sol$itr$solsta
if(status != "OPTIMAL") warning(paste("Solution status = ", status))
f <- z$sol$itr$xx[1:p]
z <- list(x = v, y = f, status = status)
class(z) <- "medde"
z
}


# Example 1 for monotonized Bayes rule estimator 

require(Rmosek)
 par(mfrow = c(1,2))
 set.seed(1984)
 n <- 100
 m <- runif(n,5,15)
 y <- rnorm(n,m)
 v <- 1:500/25
 f <- Brown(y, v, verb = 5)
 plot(f$x, f$y, type = "l",xlab = "y", ylab = "f(y)")
 x <- 1:200/10
 h <- function(x) (pnorm(15-x)-pnorm(5-x))/10 #da truth
 lines(x,h(x),col="blue")
 K <- .5 * f$x^2 + log(f$y)
 plot(f$x[-1], diff(K)/diff(f$x), type = "l",xlim = c(min(x), max(x)),
              xlab = "y", ylab = expression(delta (y) ))
 g <- function(v) (pnorm(15-v)-pnorm(5-v))/10 #da truth
 lines(v[-1],v[-1] + diff(log(g(v)))/diff(v),col = "blue")
