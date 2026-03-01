#'  Hodges-Lehmann (1952) Modification of Bayes Priors for Compound Decisions
#'  
#'  Given a prior \eqn{G} find a modified prior \eqn{f} that bounds  minimax risk.
#'
#'  There are two variants both  minimize Fisher information for location 
#'  via conic optimization:
#'   \deqn{\min \sum \frac{(f_{i+1}-f_i )^2}{( f_{i+1}+f_i )/2} = \sum \frac{u_i^2}{v_i/2} \approx I(F)}
#'   \deqn{\Leftrightarrow  \quad  \min \sum w_i \;  \mbox{s.t.} \;  u_i^2  \leq 2 v_i w_i}
#'  Huber Variant as proposed in Efron and Morris (1971) imposing constraint
#'  	\deqn{f(x) = \alpha \Phi * G + (1-\alpha) h(x)} 
#'  Mallows Variant as proposed in Bickel (1983) imposing constraints
#'   	\deqn{f(x) = \alpha \Phi * G + (1-\alpha) h(x), \; h(x) = \Phi * H}
#'  N.B. When the grid is not equispaced, one would have to include grid spacings.
#'
#' @param grid grid on which to interpolate Hodges-Lehmann solution
#' @param G initial prior (should integrate to 1)
#' @param alpha contamination proportion
#' @param type either "Huber" or "Mallows" 
#' @param sd standard deviation of the Gaussian noise
#' @param ... other arguments to be passed to Mosek.
#' @return A list containing:
#' \itemize{
#'	\code{x}: grid for domain of marginal density
#'
#'	\code{y}: function values for modified marginal density at \code{x}
#'
#'	\code{h}: function values for contamination portion at \code{x}
#'
#'	\code{d}: Bayes rule for modified prior at \code{x}
#'
#'	\code{H}: function values for contamination prior distribution, only for the "Mallows" option
#' }
#' @author R. Koenker and J. Gu
#'
#' @references 
#' Bickel, P. (1983), Minimax estimation of the mean of a normal distribution subject
#' to doing well at a point, in M. H. Rizvi, J. S. Rustagi & D. Siegmund, eds,
#' ‘Recent Advances in Statistics: Papers in Honor of Herman Chernoff on his
#' Sixtieth Birthday’, Academic Press, pp. 511–528
#'
#' Efron, B. & Morris, C. (1971), ‘Limiting the risk of Bayes and empirical Bayes
#' estimators part I: the Bayes case’, Journal of the American Statistical Association
#' 66, 807–815.
#'
#' Hodges, J. L. & Lehmann, E. L. (1952), ‘The use of previous experience in reaching
#' statistical decisions’, The Annals of Mathematical Statistics pp. 396–407.
#'
#' Huber, P. (1964), ‘Robust estimation of a location parameter’, The Annals of
#' Mathematical Statistics pp. 73–101.
#'
#' Huber, P. (1974) "Fisher Information and Spline Interpolation." 
#' Ann. Statist. 2 (5) 1029 - 1033,  
#'
#' Mallows, C. (1978), ‘Problem 78-4, minimizing an integral’, SIAM Review 20, 183– 183.
#'
#' @keywords nonparametric
#' @seealso HuberSpline
#' @importFrom Rmosek Matrix
#' @export 


HodgesLehmann <- function(grid, G, alpha, type = "Huber", sd = 1, ...){
  require(REBayes)
  if(!all.equal(sum(G$y),1)) warning("Initial G mass should sum to 1")
  gridg <- G$x
  m <- length(grid)
  p <- length(gridg)
  fG <- dnorm(outer(grid, gridg, "-"))%*%G$y
  fG <- fG/sum(fG)
  dots <- list(...)
  rtol <- ifelse(length(dots$rtol), dots$rtol, 1e-06)
  verb <- ifelse(length(dots$verb), dots$verb, 0)
  P <- list()
  P$sense <- "min"
  if(type == "Huber"){
    A = cbind(Diagonal(m), spMatrix(m, 3*m-3), -(1-alpha)*Diagonal(m))
    P$c <- c(rep(0,3*m-2), rep(1,m-1), rep(0, m))
    P$A <- rbind(A, 
               cbind(diff(Diagonal(m)), -Diagonal(m-1), spMatrix(m-1,2*m-2+m)),  
               cbind(abs(diff(Diagonal(m))), spMatrix(m-1,m-1), -Diagonal(m-1),
                     spMatrix(m-1,m-1+m)),
               cbind(matrix(1, 1,m), spMatrix(1, 4*m-3)))
    lhs <- c(alpha * fG, rep(0, 2*m-2),1)
    P$bc <- rbind(lhs,lhs)   
    P$bx <- rbind(c(rep(0,m),rep(-Inf,m-1),rep(0,2*m-2),  rep(0, m)),
                c(rep(Inf,4*m-3+m)))
    js = cbind((2*m-1)+(1:(m-1)), (3*m-2)+(1:(m-1)), m+(1:(m-1)))               
    P$F = sparseMatrix(i = 1:(3*(m-1)), j = c(t(js)), 
                     x = rep(1, 3*(m-1)), dims = c(3*(m-1), 4*m-3+m))
    P$g = c(rep(0, 3*(m-1)))
  }
  else if(type == "Mallows"){
    A = cbind(Diagonal(m), spMatrix(m, 3*m-3), -dnorm(outer(grid, gridg, "-")), spMatrix(m, p))
    P$c <- c(rep(0,3*m-2), rep(1,m-1), rep(0, 2*p))
    P$A <- rbind(A, 
        cbind(spMatrix(p, 4*m-3), Diagonal(p), -(1-alpha) * Diagonal(p)), 
        c(rep(0, 4*m-3+p), rep(1, p)),
        cbind(diff(Diagonal(m)), -Diagonal(m-1), spMatrix(m-1,2*m-2+2*p)),  
        cbind(abs(diff(Diagonal(m))), spMatrix(m-1,m-1), -Diagonal(m-1),
             spMatrix(m-1,m-1+2*p)))
    lhs <- c(rep(0, m), alpha * G$y, 1, rep(0, 2*m-2))
    P$bc <- rbind(lhs,lhs)   
    P$bx <- rbind(c(rep(0,m),rep(-Inf,m-1),rep(0,2*m-2),  rep(0, 2*p)),
	c(rep(Inf,4*m-3+2*p)))
    js = cbind((2*m-1)+(1:(m-1)), (3*m-2)+(1:(m-1)), m+(1:(m-1)))               
    P$F = sparseMatrix(i = 1:(3*(m-1)), j = c(t(js)), 
	x = rep(1, 3*(m-1)), dims = c(3*(m-1), 4*m-3+2*p))
    P$g = c(rep(0, 3*(m-1)))
  }
  else
      stop("type must be either 'Huber' or 'Mallows'")
  P$cones <-  matrix(list(),nrow=3,ncol=m-1)
  rownames(P$cones) = c("type", "dim", "conepar")
  for (k in 1:(m-1))
    P$cones[,k] <- list("RQUAD", 3, NULL)
  P$dparam$intpnt_co_tol_rel_gap <- rtol
  res <- Rmosek::mosek(P, opts = list(verbose = verb))
  x <- grid
  y = res$sol$itr$xx[1:m]
  if(type == "Huber"){
      h <-  res$sol$itr$xx[-c(1:(4*m-3))] 
      d <- x[-1] + diff(log(y))/diff(x)
      d <- c(d[1],d)
      z <- list(x = grid, y = y, h = h, d = d)
  }
  else{
    h <- dnorm(outer(grid, gridg, "-"))%*%res$sol$itr$xx[-c(1:(4*m-3+p))]
    h <- h/sum(h)
    H <- res$sol$itr$xx[-c(1:(4 * m - 3 + p))]
    d = x[-1] + diff(log(y))/diff(x)
    d = c(d[1],d)
    z <- list(x = grid, y = y, h = h, d = d, H = H)
  }
  class(z) <- "density"
  return(z)
}
