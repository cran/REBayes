par(mfrow=c(2,2))
set.seed(9)
a = 30
n = 500
grid = seq(-a,a,length = 500)
mu <- sample(c(-2,2), n, replace=TRUE, prob = c(0.5,0.5))
y <- mu + rnorm(n)
v <- seq(-5,8, by = 0.02)
G <- GLmix(y, v, rtol = 1e-10)
f0 <- HodgesLehmann(grid, G, alpha = 0.9, rtol = 1e-10)
f1 <- HodgesLehmann(grid, G, alpha = 0.9, type = "Mallows", rtol = 1e-10)

plot(f1$x, f1$h, type="l", xlab = "x", ylab = "h(x)", main = "n = 500", ylim = c(0, 0.05))
lines(f0$x, f0$h,col=2)
legend("topleft", c("Mallows", "Huber"), lty = 1, col = 1:2, cex = .7)

d1f <- approxfun(f1$x, f1$d)
d0f <- approxfun(f0$x,f0$d)
d2f <- approxfun(G$x,predict(G,G$x))

newy <- seq(-6,6,by = 0.01)

plot(newy, d1f(newy), type="l", xlab = "x", ylab = "d(x)")
lines(newy, d0f(newy), col = 2)
lines(newy, d2f(newy), col = 4)
legend("topleft", c("Mallows", "Huber", "Gauss"), col = c(1,2,4), lty = 1, cex = .7)
abline(c(0,1), lty = 3)

set.seed(9)
a = 30
n = 5000
grid = seq(-a,a,length = 500)
mu <- sample(c(-2,2), n, replace=TRUE, prob = c(0.5,0.5))
#mu <- rnorm(n, 0, 1)
y <- mu + rnorm(n)
v <- seq(-5,8, by = 0.02)
G <- GLmix(y, v, rtol = 1e-10)
f0 <- HodgesLehmann(grid, G, alpha = 0.9, rtol = 1e-10)
f1 <- HodgesLehmann(grid, G, alpha = 0.9, type = "Mallows", rtol = 1e-10)

plot(f1$x, f1$h, type="l", main = "n= 5000", xlab = "x", ylab = "h(x)")
lines(f0$x, f0$h,col=2)
legend("topleft", c("Mallows", "Huber"), lty = 1, col = 1:2, cex = .7)

d1f <- approxfun(f1$x, f1$d)
d0f <- approxfun(f0$x, f0$d)
d2f <- approxfun(G$x,predict(G,G$x))

newy <- seq(-6,6,by = 0.01)

plot(newy, d1f(newy), type="l", xlab = "x", ylab = "d(x)")
lines(newy, d0f(newy), col = 2)
lines(newy, d2f(newy), col = 4)
legend("topleft", c("Mallows", "Huber", "Gauss"), col = c(1,2,4), lty = 1, cex = .7)
abline(c(0,1), lty = 3)
