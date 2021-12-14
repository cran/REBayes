#' Predict Method for GLVmix
#' 
#' Predict Method for Gaussian Location-scale Mixtures
#' 
#' The predict method for \code{GLmix} objects will compute means, quantiles or
#' modes of the posterior according to the \code{Loss} argument.  Typically,
#' \code{newdata} would be passed to \code{predict}.  Note that these predictions
#' are for the location parameter only.
#' 
#' @param object Fitted object of class "GLVmix"
#' @param newdata data.frame with components(t,s,m) at which prediction is desired
#' @param Loss Loss function used to generate prediction:  Currently supported values: 
#' 2 to get mean predictions, 1 to get median predictions, 0 to get modal predictions
#' or any tau in (0,1) to get tau-th quantile predictions.
#' @param ... optional arguments to predict
#' @return A vector of predictions
#' @author Roger Koenker
#' @keywords nonparametric
#' @export
predict.GLVmix <- function(object, newdata, Loss = 2, ...) {
    if(missing(newdata) && Loss == 2) return(object$du)
    u <- object$u
    v <- object$v
    pu <- length(u)
    pv <- length(v)
    du <- rep(1,pu)
    dv <- rep(1,pv)
    fuv <- object$fuv
    uv <- expand.grid(u, v)
    r <- (newdata$m - 1)/2
    R <- outer(r * newdata$s, v, "/")
    G <- outer(newdata$s * gamma(r), rep(1,pv))
    r <- outer(r, rep(1, pv))
    Av <- outer((exp(-R) * R^r)/G, rep(1,pu))
    Av <- aperm(Av, c(1,3,2))
    Au <- dnorm(outer(outer(newdata$t, u, "-") * 
	outer(sqrt(newdata$m), rep(1, pu)), sqrt(v), "/"))
    Au <- Au/outer(outer(1/sqrt(newdata$m), rep(1, pu)), sqrt(v))
    A <- Av * Au
    B <- NULL
    for(j in 1:pv) B <- cbind(B, A[,,j])
    g <- as.vector(B %*% fuv)
    if(Loss == 2) { # mean case equivalent to object$dy when x == original data
	xhat <- as.vector(B %*% (uv[, 1] * fuv))/g
    }
    else if(Loss > 0 && Loss <= 1){ #quantile case
	if(Loss == 1) Loss <- 1/2
        A <- B * outer(rep(1,nrow(B)), fuv)
        B <- apply(A/apply(A,1,sum),1,cumsum) < Loss
        j <- apply(B,2,sum)
        if(any(j == 0)) { # Should only happen when v grid is very restricted
	   j <- j + 1
	   warning("zeros in posterior quantile indices")
       }
       xhat <- uv[j, 1]
    }
    else if(Loss == 0) { # mode case
	A <- B * outer(rep(1,nrow(B)), fuv)
        xhat <- uv[apply(A/apply(A,1,sum),1,which.max),1]
    }
    else 
	stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}
