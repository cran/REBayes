#' Predict Method for Bmix
#' 
#' Predict Method for Binomial Mixtures
#' 
#' The predict method for B2mix objects will compute posterior means.
#' 
#' @param object fitted object of class "B2mix"
#' @param newdata Values at which prediction is desired an n by 2 matrix
#' @param Loss Loss function used to generate prediction:  Currently supported values:
#' 2 to get mean predictions, 1 to get median predictions, 0 to get modal predictions
#' or any tau in (0,1) to get tau-th quantile predictions.
#' @param newk k values (number of trials) for the predictions an n by 2 matrix
#' @param ... optional arguments to predict
#' @return A vector of predictions
#' @author Jiaying Gu and Roger Koenker
#' @keywords nonparametric
#' @export
predict.B2mix <- function(object, newdata, Loss = 2, newk, ...) {
    x <- newdata
    u <- object$u
    v <- object$v
    y <- object$y
    k = newk
    uv <- c(outer(u,v))
    if(NCOL(x) != 2) stop("newdata must be n by 2 matrix")
    if(NCOL(k) != 2) stop("newk must be n by 2 matrix")
    if (NROW(newk)!=NROW(x)) stop("nrow(newk) must equal to nrow(newdata)")
    makeA = function(x,k,u,v)
	outer(outer(x, u, function(x,u,k) dbinom(x,k,u), k = k), rep(1,length(v)))
    A1 = makeA(x[,1], k[,1], u, v)
    A2 = makeA(x[,2], k[,2], v, u)
    A = A1 * aperm(A2, c(1,3,2))
    B = NULL
    for(j in 1:length(v)) B = cbind(B, A[,,j])
    if(Loss == 2) { # mean case equivalent to object$dy when x == original data
	xhat <- as.vector((B %*% (y * uv))/(B %*% y))
    }
    else 
	stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}




