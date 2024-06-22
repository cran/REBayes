require(REBayes)
par(mfrow = c(2,2))
x <- rt(100, 3)
alphas <- c(1, 0.5,0,-1)
main = c(expression(alpha ==  1.0),
	 expression(alpha ==  0.5),
	 expression(alpha ==  0),
	 expression(alpha ==  -1))
for(i in 1:length(alphas)){
    z <- medde(x, lambda = -1, alpha = alphas[i])
    plot(z, type = "l", xlab = "x", ylab = "f(x)", main = main[i])
}
par(mfrow = c(1,1))
