require(REBayes)
par(mfrow = c(2,2))
x <- rt(100, 3)
alphas <- c(1, 0.5,0,-1)
for(a in alphas){
    z <- medde(x, lambda = -1, alpha = a)
    plot(z, type = "l", main = paste("a = ", a))
}
par(mfrow = c(1,1))
