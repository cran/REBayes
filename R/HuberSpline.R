#' Huber (1974) Minimal Fisher (location) information spline via conic optimization
#'
#'   \deqn{\min \sum \frac{(f_{i+1}-f_i )^2}{( f_{i+1}+f_i )/2} = \sum \frac{u_i^2}{v_i/2} \approx I(F)}
#'
#'   \deqn{\Leftrightarrow  \quad  \min \sum w_i \;  \mbox{s.t.} \;  u_i^2  \leq 2 v_i w_i}
#'
#'   subject to interpolation of constraints \eqn{F(x_j) = p_j, \;  j=1,...,n.}  
#'   When \eqn{\kappa > 0}, \eqn{I(F)} is minimized within a Kolmogorov neighborhood 
#'   of the constraint points, rather than interpolating them.
#'   The generalization to Kolmogorov neighborhoods is due to Donoho
#'   and Reeves (2013).
#'
#'   N.B.  When the grid is not equispaced, one would have to include grid spacings.
#'
#' @param x quantiles to be interpolated
#' @param p probabilities associated with x
#' @param grid grid values for fitted object
#' @param kappa width of Kolmogorov neighborhood
#' @return An object of class density with solution \eqn{f*}
#' @author R. Koenker and J. Gu

#' @references P. J. Huber. (1974) "Fisher Information and Spline Interpolation." 
#' Ann. Statist. 2 (5) 1029 - 1033,  
#'
#' D. L. Donoho and G. Reeves, (2013) Achieving Bayes MMSE Performance in the Sparse Signal 
#' Gaussian White Noise Model when the Noise Level is Unknown, Proc. IEEE Symposium Istanbul, Turkey.
#'
#' @keywords nonparametric
#' @importFrom Rmosek Matrix
#' @export 

HuberSpline <- function(x, p, grid, kappa = 0){

    m <- length(grid)
    s <- c(x, max(grid))
    p <- c(p, 1)
    h <- findInterval(s,grid)
    A <- spMatrix(length(x) + 1, m)
    for(i in 1:(length(x) + 1)) A[i,1:h[i]] <- 1
    pr <- list()
    pr$sense <- "min"
    pr$c <- c(rep(0,3*m-2), rep(1,m-1))
    pr$A <- rbind(
	 cbind(diff(Diagonal(m)), -Diagonal(m-1), spMatrix(m-1,2*m-2)),  
         cbind(abs(diff(Diagonal(m))), spMatrix(m-1,m-1), -Diagonal(m-1),
               spMatrix(m-1,m-1)),     
         cbind(1, spMatrix(1,4*m-4)), 
         cbind(A, spMatrix(length(x) + 1,3*m-3)),
         c(rep(1,m), rep(0,3*m-3)))
    lhs <- c(rep(0,2*m-2),0,p - kappa,1)
    rhs <- c(rep(0,2*m-2),0,p + kappa,1)
    pr$bc <- rbind(lhs,rhs)   
    pr$bx <- rbind(c(rep(0,m),rep(-Inf,m-1),rep(0,2*m-2)),
               c(rep(Inf,4*m-3)))
    js = cbind((2*m-1)+(1:(m-1)), (3*m-2)+(1:(m-1)), m+(1:(m-1)))               
    pr$F = sparseMatrix(i = 1:(3*(m-1)), j = c(t(js)), 
		    x = rep(1, 3*(m-1)), dims = c(3*(m-1), 4*m-3))
    pr$g = c(rep(0, 3*(m-1)))
    pr$cones <-  matrix(list(),nrow=3,ncol=m-1)
    rownames(pr$cones) = c("type", "dim", "conepar")
    for (k in 1:(m-1))
	pr$cones[,k] <- list("RQUAD", 3, NULL)
    res <- Rmosek::mosek(pr)
    z <- list(x = grid, y = res$sol$itr$xx[1:m])
    class(z) <- "density"
    return(z)
    }
