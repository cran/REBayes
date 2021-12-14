# Bradley Terry Estimation with Ranking Lasso Regularization
dgp <- function(a,m){
    # Roundrobin tournament: n teams, m rounds
    n <- length(a)
    A <- combn(n, 2)
    i <- factor(rep(A[1,], each = m), levels = 1:n)
    j <- factor(rep(A[2,], each = m), levels = 1:n)
    p <- a[i]/(a[i] + a[j])
    y <- (runif(length(i)) < p) * 1
    data.frame(i = i, j = j, y = y)
}
set.seed(123)
n = 20
m = 10
a <- exp(rnorm(n))
D <- dgp(a, m)
y <- D$y
X <- with(D, model.matrix(~ i - 1) - model.matrix(~ j - 1))
X <- X[,-1] # 1 is the reference team

# Construct Ranking Penalty Matrix
teams <- factor(levels(D$i)[-1])
pairs <- t(combn(teams,2)) 
P <- model.matrix(~ pairs[,1] - 1) - model.matrix(~ pairs[,2] - 1)
P <- rbind(diag(length(teams)),P)

lambdas <- 1:10/2
groups <- lambdas
for(i in 1:length(lambdas)){
    f = RLR(X, y, P, lambdas[i])
    groups[i] <- length(unique(round(c(0,f$coef),4)))
}
plot(lambdas, groups, type = "b")

