#Demo of Huber Spline function
x <- sort(rchisq(5,3))
p <- pchisq(x,3)
f0 <- HuberSpline(x, p, seq(0,16,length = 500), kappa = 0)
f1 <- HuberSpline(x, p, seq(0,16,length = 500), kappa = 0.02)
par(mfrow = c(1,2))
plot(f0, xlab = "x", main = "")
lines(f1, col = 2)
plot(f0$x,cumsum(f0$y),type = "l", xlab = "x", ylab = "F(x)")
lines(f1$x, cumsum(f1$y), col = 2)
points(x,p)

