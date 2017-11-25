
# Toy problem with two mass points

require(REBayes)
uv0 <- rbind(rep(pi/2,2),rep(pi/4,2))
B0 <- cbind(sin(uv0[,1]) * cos(uv0[,2]),
       sin(uv0[,1]) * sin(uv0[,2]),
       cos(uv0[,1]))
n <- 500
m = 50
x <- matrix(rnorm(2 * n), n, 2)
XB0 <- cbind(1,x) %*% t(B0)
s <- sample(0:1, n, replace = TRUE)
utility <- s * XB0[,1] + (1-s) * XB0[,2]
y <- (utility > 0) - 0
f <- Cosslett2(x, y, m = m, rtol = 1e-10, verb = 5)
u <- v <- seq(0, pi, length = m)
contour(u,v, matrix(f$y,m,m))
points(uv0, col = 2)

