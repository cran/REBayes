KWDual <- function(x, w, d, A, rtol = 1.0e-6, verb = 0){
############################################################################
# Useage:  KWDual(x, w, d, A, rtol = 1.0e-6, verb = 0)
#
# Dual Kiefer-Wolfowitz MLE for Mixture Problems
#
#       This version implements a class of density estimators solving:
#
#       min_x  {F(x) := sum -log (x_i)}  s.t. A' x <= d, 0 <= x,  
#
#
#	where e.g.  A = phi(outer(Y,g,"-")), with Y data and g a grid on the support of Y.
#
#	This is an Rmosek version of an earlier MeddeR function
#-------------------------------------------------------------------------------------
#
# Roger Koenker 
#
# First   version: 24 Feb 2012  
############################################################################


n <- nrow(A)
m <- ncol(A)
A <- t(A) 
C <- rep(0,n)
P <- list(sense = "min")
P$c <- C
P$A <- A
P$bc <- rbind(rep(0,m),d)
P$bx <- rbind(rep(0,n),rep(Inf,n))


opro <- matrix ( list (), nrow =5, ncol = n)
rownames ( opro ) <- c(" type ","j","f","g","h")

opro[1,] <-  as.list(rep('log',n))
opro[2,] <-  as.list(1:n)
opro[3,] <-  as.list(-w)
opro[4,] <-  as.list(rep(1,n))
opro[5,] <-  as.list(rep(0,n))
P$scopt<- list(opro = opro)
P$dparam$intpnt_nl_tol_rel_gap <- rtol
z <- mosek(P, opts = list(verbose = verb))
status <- z$sol$itr$solsta
g <- z$sol$itr$xx
f <- z$sol$itr$suc
logLik <- sum(w * log(w/(g * sum(f))))
list(f = f, g = g, logLik = logLik, status = status)
}
