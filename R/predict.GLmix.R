#' Predict Method for GLmix
#' 
#' Predict Method for Gaussian Location Mixtures
#' 
#' The predict method for \code{GLmix} objects will compute means, quantiles or
#' modes of the posterior according to the \code{Loss} argument.  Typically,
#' \code{newdata} would be passed to \code{predict}
#' 
#' @param object fitted object of class "GLmix"
#' @param newdata Values at which prediction is desired
#' @param Loss Loss function used to generate prediction:  Currently supported values: 
#' 2 to get mean predictions, 1 to get median predictions, 0 to get modal predictions
#' or any tau in (0,1) to get tau-th quantile predictions.
#' @param newsigma sigma values for the predictions 
#' @param ... optional arguments to predict
#' @return A vector of predictions
#' @author Roger Koenker
#' @keywords nonparametric
#' @export
predict.GLmix <- function(object, newdata, Loss = 2, newsigma = NULL, ...) {
    x <- newdata
    n <- length(x)
    v <- object$x
    fv <- object$y
    if (length(newsigma)) object$sigma = newsigma
    if(Loss == 2) { # mean case equivalent to object$dy when x == original data
	A <- dnorm(outer(x, v, "-"), sd = object$sigma)
	xhat <- as.vector((A %*% (fv * v))/(A %*% fv))
    }
    else if(Loss > 0 && Loss <= 1){ #quantile case
	if(Loss == 1) Loss <- 1/2
       A <- dnorm(outer(x, v, "-"), sd = object$sigma) * outer(rep(1,n),fv)
       B <- apply(A/apply(A,1,sum),1,cumsum) < Loss
       j <- apply(B,2,sum)
       if(any(j == 0)) { # Should only happen when v grid is very restricted
	   j <- j + 1
	   warning("zeros in posterior median indices")
       }
       xhat <- v[j]
    }
    else if(Loss == 0) { # mode case
       A <- dnorm(outer(x, v, "-"), sd = object$sigma) * outer(rep(1,n), fv)
       xhat <- v[apply(A/apply(A,1,sum),1,which.max)]
    }
    else 
	stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}
