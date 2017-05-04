# Quasi concave density demo of medde()
require(REBayes)
data(velo)
x <- velo[!(velo == 0)] #discard zeros
hist(x, 100, freq = FALSE, main = "Histogram of Star Velocity and Three Quasi Fits")
f <- medde(x, v = 1000, lambda = -0.5, alpha = 0.5, rtol = 1e-8, mass = 1)
lines(f$x, f$y, col = "red")
f <- medde(x, v = 1000, lambda = -0.5, alpha = 0, rtol = 1e-8, mass = 1)
lines(f$x, f$y, col = "blue")
f <- medde(x, v = 1000, lambda = -0.5, alpha = -2, rtol = 1e-8, mass = 1)
lines(f$x, f$y, col = "green")
leg <- as.expression(lapply(c(1/2,0,-2),function(x) bquote(alpha ==.(x))))
legend(300,0.025,leg, lty = 1, col = c("red","blue","green"))



