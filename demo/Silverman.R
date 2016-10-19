# Test of Silverman log-spline estimator in medde
# Silverman BW (1982). On the Estimation of a Probability Density Function
# by the Maximum Penalized Likelihood Method. Annals of Statistics, 10, 795-810
# Moral:  lambda selection can be tricky!!

require(REBayes)
n <- 500
x <- rnorm(n)
main = "Histogram and two Silverman Estimates"
hist(x, 70, freq = FALSE, main = main, col = grey(.9))
f <- medde(x, Dorder = 2, lambda = 0.005, verb = 5)
lines(f, col = "red")
f <- medde(x, Dorder = 2, lambda = 0.0000005, verb = 5)
lines(f, col = "blue")
