# A Model for Shakespeare's vocabulary
require(deconvolveR) # For data
data(bardWordCount)
w <- bardWordCount
x <- 1:100
x <- rep(1:100, times = w)
v <- exp(seq(-8, 4.5, 0.005)) 
f <- Pmix(x, v = v, support = c(0,100), rtol = 1e-14)
plot(f$x, cumsum(f$y), type = "l", xlab = expression(lambda), ylab = expression(G[n] (lambda)))
title("Mixing Distribution for Word Frequency")



