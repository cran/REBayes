
predict.GLmix <- function(object, newdata, Loss = 2, ...) {
    #Given a fitted Gaussian Location Mixture predict at the vector x:
    #	object is a fitted GLmix object
    #	newdata  the points at which predictions are desired
    #   Loss is the p of the Lp loss function
    x <- newdata
    v <- object$x
    fv <- object$y

    if(Loss == 2) { # mean case equivalent to object$dy when x == original data
	A <- dnorm(outer(x, v, "-"), sd = object$sigma)
	xhat <- as.vector((A %*% (fv * v))/(A %*% fv))
    }
    else if(Loss == 1){ #median case
       A <- t(t(dnorm(outer(x, v, "-"), sd = object$sigma)) * fv)
       B <- apply(A/apply(A,1,sum),1,cumsum) < 1/2
       j <- apply(B,2,sum)
       if(any(j == 0)) { # Should only happen when v grid is very restricted
	   j <- j + 1
	   warning("zeros in posterior median indices")
       }
       xhat <- v[j]
    }
    else if(Loss == 0) { # mode case
       A <- t(t(dnorm(outer(x, v, "-"), sd = object$sigma)) * fv)
       xhat <- v[apply(A/apply(A,1,sum),1,which.max)]
    }
    else 
	stop(paste("Loss", Loss, "not (yet) implemented"))
    xhat
}
