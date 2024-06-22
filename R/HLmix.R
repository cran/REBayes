#' Kiefer-Wolfowitz NPMLE for Huber Location Mixtures
#' 
#' Kiefer Wolfowitz Nonparametric MLE for Huber Location Mixtures
#' 
#' Kiefer Wolfowitz NPMLE for location mixtures with Huber (1964) base density
#' The Huber \code{k} specifies the point at which the influence function of
#' the Huber M-estimator kinks.  
#' The predict method for \code{HLmix} objects compute means, medians or
#' modes of the posterior according to whether the \code{Loss} argument is 2, 1
#' or 0, or posterior quantiles if \code{Loss} is in (0,1).
#' 
#' @param x  Data: Sample Observations 
#' @param v  Undata: Grid Values defaults equal spacing of with v bins, when v is
#' a scalar 
#' @param sigma  scale parameter of the Gaussian noise, may take vector values
#' of length(x) 
#' @param k  Huber k value 
#' @param heps  Huber epsilon contamination value, should match k, by default
#' this is automatically enforced. 
#' @param ...  other parameters to pass to KWDual to control optimization 
#' @return An object of class density with components: 
#'      \item{x}{points of  evaluation on the domain of the density} 
#'      \item{y}{estimated function values at the points v, the mixing density} 
#'      \item{g}{marginal density values}
#'      \item{logLik}{log likelihood}
#'      \item{sigma}{sigma}
#'      \item{dy}{posterior means at the observed \code{x} values}
#'      \item{k}{Huber k}
#'      \item{heps}{Huber epsilon}
#' @author Roger Koenker
#' @export
                                                                                            

HLmix <- function(x, v = 300, sigma = 1, k = 1.345, heps = hubereps(k), ...){
    n = length(x)
    eps <- 1e-04
    if (length(v) == 1) 
        v <- seq(min(x) - eps, max(x) + eps, length = v)
    m <- length(v)
    w <- rep(1,n)/n
    d <- rep(1, length(v))
    A <- dhuber(outer(x, v, "-"), sigma = sigma, k = k, heps = heps)
    f <- KWDual(A, d, w, ...)
    y <- f$f
    g <- f$g
    logLik <- n * sum(w * log(g))
    dy <- as.vector((A %*% (y * d * v))/g)
    z <- list(x = v, y = y, g = g, logLik = logLik, sigma = sigma, 
        dy = dy, k = k, heps = heps, status = f$status)
    class(z) <- c("HLmix", "density")
    return(z)
}
#' Predict Method for HLmix
#' 
#' Predict Method for Huber Location Mixtures
#' 
#' The predict method for \code{HLmix} objects computes means, quantiles or
#' modes of the posterior according to the \code{Loss} argument.  Typically,
#' \code{newdata} would be passed to \code{predict}.  Note that if newdata
#' is simply equal to the original observations (denoising case) then the
#; dy component of object contains the posterior means.
#' 
#' @param object  fitted object of class "HLmix" 
#' @param newdata  Values at which prediction is desired 
#' @param Loss  Loss function used to generate prediction:  Currently supported values: 
#' 2 to get mean predictions, 1 to get median predictions, 0 to get modal predictions
#' or any tau in (0,1) to get tau-th quantile predictions. 
#' @param newsigma  sigma values for the predictions 
#' @param ...  optional arguments to predict 
#' @return A vector of predictions
#' @author Roger Koenker
#' @export
#'
predict.HLmix <- function (object, newdata, Loss = 2, newsigma = NULL, ...) 
{
    x <- newdata
    n <- length(x)
    v <- object$x
    k <- object$k
    heps <- object$heps
    fv <- object$y
    if (length(newsigma)) 
        object$sigma = newsigma
    A <- dhuber(outer(x, v, "-"), sigma = object$sigma, k = k, heps = heps)
    if (Loss == 2) {
        xhat <- as.vector((A %*% (fv * v))/(A %*% fv))
    }
    else if (Loss > 0 && Loss <= 1) {
        if (Loss == 1) 
            Loss <- 1/2
        B <- apply(A/apply(A, 1, sum), 1, cumsum) < Loss
        j <- apply(B, 2, sum)
        if (any(j == 0)) {
            j <- j + 1
            warning("zeros in posterior median indices")
        }
        xhat <- v[j]
    }
    else if (Loss == 0) {
        xhat <- v[apply(A/apply(A, 1, sum), 1, which.max)]
    }
    else stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}

