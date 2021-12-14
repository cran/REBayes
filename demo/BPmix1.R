# A simple test problem to explore bivariate Binomial-Poisson mixture model

# Model and Data Structure

set.seed(14)
n = 5000
p = sample(c(.1,.5,.9), n, prob = c(.2, .6, .2), replace = TRUE)
m = (p < .5) * rpois(n,10) +  (p == .5) * rpois(n,15) +  (p > .5) * rpois(n,10)
y = rbinom(n,m,p)

# Now try fitting and visualization
z <- BPmix(y, m, verb = 5)
g <- expand.grid(v = z$v, u = z$u)
g$y <- z$y
require(lattice)
pl <- cloud(y ~ u * v, data = g, type = "h", lwd = 2, 
        zlim = c(0, max(z$y)), scales = list(arrows = FALSE,
        xlab = expression(u), ylab = expression(v), zlab = "density",
        screen = list(z = 10, x = -70)))
print(pl)





