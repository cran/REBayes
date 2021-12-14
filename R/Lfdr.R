#' Local False Discovery Rate Estimation
#'
#' A Generic function for estimation of Local FDR
#'
#' Given an estimated mixing distribution, G, Lfdr computes
#' an estimated local false discovery rate at a specified set
#' of points and threshold value cnull.  The argument G can be
#' specified as the fitted object from one of several possible
#' fitting routines for nonparametric mixing distributions.
#' 
#' 
#' @param G A fitted object from some G-modeling function.
#' @param newdata data frame to in which to evaluate Lfdr
#' @param cnull threshold for evaluation of Lfdr
#' @param tail either "R" or "L" to specify tail focus
#' @param ... other arguments 
#' @export
Lfdr <- function(G, ...) UseMethod("Lfdr")
#' @rdname Lfdr
#' @export
Lfdr.GLVmix <- function(G, newdata, cnull, tail = "R", ...){
    if(missing(newdata)) A = G$A
    else{
	with(newdata,{
	    n <- length(y)
	    wsum <- tapply(w, id, "sum")
	    t <- tapply(w * y, id, "sum")/wsum
	    m <- tapply(y, id, "length")
	    r <- (m - 1)/2
	    s <- (tapply(w * y^2, id, "sum") - t^2 * wsum)/(m - 1)
	    n <- length(s)
	    u = G$u
	    v = G$v
	    pu = length(u)
	    pv = length(v)
	    fuv = G$fuv
	    uv <- expand.grid(alpha = u, theta = v)
	    R <- outer(r * s, v, "/")
	    G <- outer(s * gamma(r), rep(1, pv))
	    r <- outer((m - 1)/2, rep(1, pv))
	    Av <- outer((exp(-R) * R^r)/G, rep(1, pu))
	    Av <- aperm(Av, c(1, 3, 2))
	    Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum), rep(1, 
		pu)), sqrt(v), "/"))
	    Au <- Au/outer(outer(1/sqrt(wsum), rep(1, pu)), sqrt(v))
	    Auv <- Av * Au
	    A <- NULL
	    for (j in 1:pv) A <- cbind(A, Auv[, , j])})
    }
    uv <- expand.grid(G$u, G$v)
    if(tail == "R")  z = 1 - c((A %*% (G$fuv * (uv[,1] < cnull)))/(A %*% G$fuv))
    else z = 1 - c((A %*% (G$fuv * (uv[,1] >= cnull)))/(A %*% G$fuv))
    z
}
#' @rdname Lfdr
#' @export
Lfdr.WGLVmix <- function(G, newdata, cnull, tail = "R", ...){
    if(missing(newdata)) A = G$A
    else{
	with(newdata,{
	    n <- length(y)
	    wsum <- tapply(w, id, "sum")
	    t <- tapply(w * y, id, "sum")/wsum
	    m <- tapply(y, id, "length")
	    r <- (m - 1)/2
	    s <- (tapply(w * y^2, id, "sum") - t^2 * wsum)/(m - 1)
	    n <- length(s)
	    u = G$u
	    v = G$v
	    pu = length(u)
	    pv = length(v)
	    fuv = G$fuv
	    uv <- expand.grid(alpha = u, theta = v)
	    R <- outer(r * s, v, "/")
	    G <- outer(s * gamma(r), rep(1, pv))
	    r <- outer((m - 1)/2, rep(1, pv))
	    Av <- outer((exp(-R) * R^r)/G, rep(1, pu))
	    Av <- aperm(Av, c(1, 3, 2))
	    Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum), rep(1, 
		pu)), sqrt(v), "/"))
	    Au <- Au/outer(outer(1/sqrt(wsum), rep(1, pu)), sqrt(v))
	    Auv <- Av * Au
	    A <- NULL
	    for (j in 1:pv) A <- cbind(A, Auv[, , j])})
    }
    uv = expand.grid(G$u, G$v)
    if(tail == "R")  z = 1 - c((A %*% (G$fuv * (uv[,1] < cnull)))/(A %*% G$fuv))
    else z = 1 - c((A %*% (G$fuv * (uv[,1] >= cnull)))/(A %*% G$fuv))
    z
}


#' @rdname Lfdr
#' @export
Lfdr.GLmix <- function(G, newdata, cnull, tail = "R", ...){ # changed for left tail selection
    v = G$x
    fv = G$y
    if(missing(newdata)) A = G$A
    else{
	with(newdata,{
	    x = newdata$x
	    s = newdata$sigma
	    A = dnorm(outer(x, v, "-"), sd = s)
	    })
    }
    if(tail == "R")  v = 1 - c((A %*% (fv * (v < cnull)))/(A %*% fv))
    else v = 1 - c((A %*% (fv * (v >= cnull)))/(A %*% fv))
    v
}


