traprule <- function(x,y) sum(diff(x) * (y[-1] + y[-length(y)]))/2
L1norm <- function(F,G, eps = 1e-6) {
        if(!is.stepfun(F) || !is.stepfun(G)) stop("Both F and G must be stepfun")
        xk <- sort(c(knots(F), knots(G)))
        n <- length(xk)
        if(!all.equal(F(xk[1] - eps), G(xk[1] - eps))) stop("F(x[1]-) != G(x[1]-)")
        if(!all.equal(F(xk[n]), G(xk[n]))) stop("F(x[n]) != G(x[n])")
        dy <- (abs(F(xk) - G(xk)))[-length(xk)]
        dx <- diff(xk)
        sum(dy * dx)
        }

