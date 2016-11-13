# Quasi concave density demo of medde()
require(REBayes)
data(velo)
x <- velo[!(velo == 0)] #discard zeros
v <- 0:460 + 0.05
w <- tabulate(x)
w <- w/sum(w)
y <- 1:length(w)
hist(x, 100, freq = FALSE, main = "Histogram of Star Velocity and Some Quasi Fits")
f <- medde(y, v = 1000, lambda = -0.5, alpha = 0.5, w = w, rtol = 1e-8)
lines(f$x, f$y, col = "red")
f <- medde(y, v = 1000, lambda = -0.5, alpha = 0, w = w, rtol = 1e-8)
lines(f$x, f$y, col = "blue")
f <- medde(y, v = 1000, lambda = -0.5, alpha = -1, w = w, rtol = 1e-8)
lines(f$x, f$y, col = "green")
leg <- as.expression(lapply(c(1/2,0,-1),function(x) bquote(alpha ==.(x))))
legend(300,0.025,leg, lty = 1, col = c("red","blue","green"))



