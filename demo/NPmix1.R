# A simple test problem to explore bivariate Normal-Poisson model

# %Model and Data Structure
set.seed(14)
n = 5000
a = sample(c(6,7.5,9), n, prob = rep(1/3,3), replace = TRUE)
m = (a < 7.5) * rpois(n,30) +  (a == 7.5) * rpois(n,60) +  (a > 7.5) * rpois(n,30)
y = rnorm(n,a, 1/sqrt(m))

# Now try fitting and visualization
z <- NPmix(y, m, verb = 5)
g <- expand.grid(v = z$v, u = z$u)
g$y <- z$y
require(lattice)
pl <- cloud(y ~ u * v, data = g, type = "h", lwd = 2, 
        zlim = c(0, max(z$y)), scales = list(arrows = FALSE,
        xlab = expression(u), ylab = expression(v), zlab = "density",
        screen = list(z = 10, x = -70)))
print(pl)





