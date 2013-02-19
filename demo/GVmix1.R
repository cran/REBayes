# Compound decisions for Gaussian Variances  

# Initialize Rmosek
require(Rmosek)

n <- 100
r <- 10
m <- 2*r + 1
theta <- rep(c(.2,.3),n/2)
y <- rnorm(n*m, mean = 0, sd = rep(theta, each = m))
id <- rep(1:n,each = m)
f <- GVmix(y,id, verb = 0)
plot(f,xlab = expression(sigma^2), main = "Estimated Mixing Density")
abline(v = .04, col = "red")
abline(v = .09, col = "red")
