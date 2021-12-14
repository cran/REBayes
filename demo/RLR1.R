# A Simple Lasso Plot for Logistic Regression
set.seed(1729)
n = 100
p = 10
X <- matrix(rnorm(n*p),n,p)
y <- sample(0:1, n, replace = TRUE)
lambdas <- 1:25/3
B <- matrix(0, length(lambdas),10)
for(i in 1:length(lambdas))
    B[i,] <- RLR(X,y,diag(p),lambdas[i])$coef
matplot(lambdas, B, type = "l", xlab = expression(lambda),ylab = "Coefficients")
title("A Simple Logistic Regression Lasso Plot")

