# Test DGP for B2mix
n = 1000
p = rbind(c(0.5, 0.75),c(0.33, 0.66))
s = sample(1:2, n, replace = TRUE)
k = matrix(sample(8:12,2*n, replace = TRUE), n,2)
x = matrix(NA, n, 2)
for(i in 1:n){
    for(j in 1:2){
	x[i,j] = rbinom(1,k[i,j],p[s[i],j])
    }
}
OraclelogLik = sum(log(0.5 * dbinom(x[,1],k[,1],p[1,1]) * dbinom(x[,2],k[,2],p[1,2]) + 
	0.5 * dbinom(x[,1],k[,1],p[2,1]) * dbinom(x[,2],k[,2],p[2,2])))
require(Rmosek)
f = B2mix(x, k, verb = 5)
uv = expand.grid(f$u,f$v)
puv = cbind(uv[f$y > 0.001,],f$y[f$y > 0.001])

# Now try to plot as in demo(WGLVmix1)
g <- expand.grid(u = f$u, v = f$v)
g$y <- f$y
require(lattice)
pl <- cloud(y ~ u * v, data = g, type = "h", lwd = 2, 
       zlim = c(0, max(f$y)), scales = list(arrows = FALSE,
       xlab = expression(u), ylab = expression(v), zlab = "density",
       screen = list(z = 10, x = -70)))
print(pl)
