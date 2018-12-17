#' Predict Method for Pmix
#' 
#' Predict Method for Poisson Mixtures
#' 
#' The predict method for \code{Pmix} objects will compute means, quantiles or
#' modes of the posterior according to the \code{Loss} argument.  Typically,
#' \code{newdata} would be passed to \code{predict}
#' 
#' @param object fitted object of class "Pmix"
#' @param newdata Values at which prediction is desired
#' @param Loss Loss function used to generate prediction.  Currently supported values:
#' 2 to get mean predictions, 1 to get median predictions, 0 to get modal predictions
#' or any tau in (0,1) to get tau-th quantile predictions.
#' @param newexposure exposure values for the predictions 
#' @param ... optional arguments to predict
#' @return A vector of predictions
#' @author Jiaying Gu
#' @keywords nonparametric
#' @importFrom stats dpois
#' @export
predict.Pmix <- function(object, newdata, Loss = 2, newexposure = NULL, ...) {
    x <- newdata
    n = length(x)
    v <- object$x
    fv <- object$y
    if (length(newexposure)){
    	A0 <- outer(x,v,"dpois")}
    else if (length(newexposure)==n){
    	A0 = matrix(0, n, length(v))
    	for (i in 1:n) A0[i,] = sapply(x[i], v*newexposure[i], FUN = dpois)
    	}
    else 	stop("length(newexposure) must equal to length(newdata)")
    if(Loss == 2) { # mean case equivalent to object$dy when x == original data
	A <- A0
	xhat <- as.vector((A %*% (fv * v))/(A %*% fv))
    }
    else if(Loss > 0 && Loss <= 1){ #median case
	if(Loss == 1) Loss <- 1/2
       A <- t(t(A0) * fv)
       B <- apply(A/apply(A,1,sum),1,cumsum) < Loss
       j <- apply(B,2,sum)
       if(any(j == 0)) { # Should only happen when v grid is very restricted
	   j <- j + 1
	   warning("zeros in posterior median indices")
       }
       xhat <- v[j]
    }
    else if(Loss == 0) { # mode case
       A <- t(t(A0) * fv)
       xhat <- v[apply(A/apply(A,1,sum),1,which.max)]
    }   
    else 
	stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}
