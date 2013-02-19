# Compound decisions for Gaussian Means and Variances  

require(Rmosek)

n <- 100
r <- 10
m <- 2*r + 1
theta <- rep(c(1,2),n/2)
mu <- rep(c(2,3),each = n/2)
y <- rnorm(n*m, mean = rep(mu, each = m), sd = rep(theta, each = m))
id <- rep(1:n,each = m)
f <- GLVmix(y,id, verb = 5)
X11(width = 8, height = 5)
par(mfrow = c(1,2))
plot(f$u,f$fu,type = "l")
plot(f$v,f$fv,type = "l")
