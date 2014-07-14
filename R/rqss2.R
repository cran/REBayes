rqss2 <- function(x, y, tau = .5, lambda = 1, dual = TRUE, rtol = 1e-06, verb = 0){
# Smoothing spline estimation for QR with L2 penalty  Dual Version
n <- length(x)
if(length(y) != n) stop("lengths of x and y don't match")

MakeQR <- function(x) {
	u <- sort(unique(x))
	h <- diff(u)
	n <- length(h)
	R <- bandSparse(n-1, n-1, -1:1, list(h[-1]/3, c(2*(h[-1] + h[-n])/3, 0), h[-1]/3))
	Q <- bandSparse(n+1,n-1, -2:0, list(1/h[-1], c(-(1/h[-1] + 1/h[-n]), 0), 1/h[-n]))
	list(Q = Q, R = R)
	}
G <- MakeQR(x)
T <- backsolve(chol(G$R), Diagonal(n-2)) %*% t(G$Q)

#  Now formulate a Mosek problem to finish the job...

QP <- list()
QP$sense <- "min"
QP$dparam$intpnt_nl_tol_rel_gap <- rtol
if(dual){
   QP$c <- c(rep(0, n-2), -y)
   QP$A <- cBind(sqrt(lambda) * t(T), -Diagonal(n))
   QP$bc <- rbind(blc = rep(0,n), buc = rep(0,n)) 
   QP$bx <- rbind(blx = c(rep(-Inf,n-2), rep(tau - 1, n)), 
	       bux = c(rep(Inf,n-2), rep(tau,n)))
   QP$qobj <- list(i = 1:(n-2), j = 1:(n-2), v = rep(1,n-2))
   r <- mosek(QP, opts = list(verbose = verb))
   fit <- r$sol$itr$slc - r$sol$itr$suc
   }
else{
   D <- Diagonal(n)
   D2 <- Diagonal(n-2)
   M1 <- Matrix(0,n,n-2)
   M2 <- Matrix(0,n-2,n)
   QP$c <- c(rep(0, 2*n-2), rep(tau,n), rep(1-tau,n))
   QP$A <- rBind(cBind(M1,D,D,-D), cBind(D2,-T,M2,M2))
   QP$bc <- rbind(blc = c(y, rep(0,n-2)), buc = c(y, rep(0,n-2))) 
   QP$bx <- rbind(blx = c(rep(-Inf,2*n-2), rep(0, 2*n)), 
	          bux = c(rep(Inf,2*n-2), rep(Inf, 2*n)))
   QP$qobj <- list(i = 1:(n-2), j = 1:(n-2), v = rep(lambda,n-2))
   r <- mosek(QP, opts = list(verbose = verb))
   fit <- r$sol$itr$xx[(n-1):(2*n-2)]
}
status <- r$sol$itr$solsta
list(fit = fit, status = status)
}
