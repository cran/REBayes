# Test of Silverman log-spline estimator in medde
# Silverman BW (1982). On the Estimation of a Probability Density Function
# by the Maximum Penalized Likelihood Method. Annals of Statistics, 10, 795-810

require(REBayes)
n <- 500
x <- rnorm(n)
main = "Histogram and two Silverman Estimates"
hist(x, 70, freq = FALSE, main = main, col = grey(.9))
f <- medde(x, Dorder = 2, lambda = 1, verb = 5, mass = 1)
lines(f, col = "red")
f <- medde(x, Dorder = 2, lambda = 0.05, verb = 5, mass = 1)
lines(f, col = "blue")
