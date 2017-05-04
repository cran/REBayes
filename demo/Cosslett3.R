# Profile likelihood for Cosslett Estimator.
n <- 2000
x1 <- rnorm(n)
x2 <- rnorm(n)
util = x1 * 0.6 + x2 * 0.8

eta <- sample(c(-1,1),n,replace=TRUE, prob = c(1/2,1/2))
y = (eta < util)-0
bs <- 1:38/40
ll <- bs
for (r in 1:length(bs)){
    x <- bs[r] * x1 + x2 * sqrt(1- bs[r]^2)
    fit <- Cosslett(x = x, y = y)
    ll[r] <- fit$logL
}
plot(bs, ll, type="o", cex = 0.3, main = "Profile Likelihood for beta1 (n = 2000)", 
     xlab = "b1", ylab = "logL")



