

# Example 1 for monotonized Bayes rule estimator 

 require(REBayes)
 par(mfrow = c(1,2))
 set.seed(1984)
 n <- 100
 m <- runif(n,5,15)
 x <- rnorm(n,m)
 v <- 1:500/25
 f <- medde(x, v, lambda = -2, verb = 5)
 plot(f$x, f$y, type = "l",xlab = "x", ylab = "g(x)")
 x <- 1:200/10
 h <- function(x) (pnorm(15-x)-pnorm(5-x))/10 #da truth
 lines(x,h(x),col="blue")
 K <- .5 * f$x^2 + log(f$y)
 plot(f$x[-1], diff(K)/diff(f$x), type = "l",xlim = c(min(x), max(x)),
              xlab = "y", ylab = expression(delta (y) ))
 g <- function(v) (pnorm(15-v)-pnorm(5-v))/10 #da truth
 lines(v[-1],v[-1] + diff(log(g(v)))/diff(v),col = "blue")
