# Poisson mixture prediction
# Illustrates that the weighted loss function of Clevenson and Zidek
# which implies the posterior harmonic mean as Bayes rule is closely
# linked to quadratic loss and the posterior mean Bayes rule.
# NB. You can try to trick the default predict command with Loss = 2
# into computing the posterior harmonic mean by evaluating with
# newdata = X - 1, but if there are any X = 0 observations, they produce NaNs. 
n = 500
lam = runif(n, 0.5, 15)
X = rpois(n, lam)
f = Pmix(X)
d2 = predict(f, newdata = X, newexposure = rep(1,n))
d1 = predict(f, newdata = X+1, Loss = 1, newexposure = rep(1,n))
plot(d1,d2, cex = 0.5, xlab = expression(delta[1](X)), ylab = expression(delta[2](X)))
abline(c(0,1))
title(expression(delta[1] (X + 1) == delta[2] (X)))
