#' Plot a GLVmix object
#'
#' Given a fitted mixture model by GLVmix plot the estimated mass points
#'
#' @param  x is the fitted object
#' @param  ... other arguments to pass to \code{symbols}, notably e.g. \code{add = TRUE}
#' @return nothing (invisibly)
#' @importFrom lattice cloud
#' @export
plot.GLVmix <- function(x, ...){
    g <- expand.grid(alpha = x$u, theta = x$v)
    g$fuv <- x$fuv
    pl <- lattice::cloud(fuv ~ alpha * theta, data = g, type = "h", lwd = 2,
       zlim = c(0, max(g$fuv)), scales = list(arrows = FALSE,
       xlab = "\u03B1", ylab = "\u03B8", zlab = "density",
       screen = list(z = 10, x = -70)), ...)
    print(pl)
}
