# Comparison of classical least concave majorant implementation of the
# Grenander estimator of a decreasing density with two direct optimization
# methods: one based on the function medde from the REBayes package, the 
# other based on the NPMLE of Kiefer and Wolfowitz for uniform scale mixtures.
# The functions LCM and Grenander are slightly edited versions of code from 
# https://github.com/pietg/book/blob/master/Simple_Rscripts/Fig2-4.R

LCM <- function(x){ # Least Concave Majorant 
    n = length(x)
    x = sort(x)
    y = rep(0, n+1)
    s = seq(1/n, 1, by = 1/n)
    y[1] = s[1]/x[1]
    for(i in 2:n){
	y[i] = (s[i] - s[i-1])/(x[i] - x[i-1])
	if(y[i-1] < y[i]){
	    j = i
	    while(y[j-1] < y[j] && j > 1){
		j = j-1
		if(j > 1) 
		    y[i] = (s[i] - s[j-1])/(x[i] - x[j-1])
		else
		    y[i] = s[i]/x[i]
		for(k in j:(i-1))
		    y[k] = y[i]
	    }
	}
    }
    y
}

Grenander = function(x){
    n = length(x)
    x = sort(x)
    y = LCM(x)
    y[n+1] = 0
    j = 0
    for(i in 1:n)
	if (y[i+1] < y[i]) j = j+1
    m = j
    t = u = rep(0,m)
    j = 0
    for (i in 1:n) {
	if (y[i+1] < y[i]) {
		j = j+1
		t[j] = x[i]
		u[j] = y[i]
	}
    }
    z = list(x = x, y = y, t = t, u = u)
    class(z) = "Grenander"
    return(z)
}


plot.Grenander = function(x, ...){
    plot(c(x$x,max(x$x)),x$y,xlab="x",ylab="f(x)",type="n")
    points(x$t,x$u,pch=19, cex = .5)
    lines(c(0,x$x,max(x$x)),c(x$y[1],x$y),lwd=2,type = "S",col="blue")
}

n = 50
x = rexp(n, 1) 
x = runif(n, 0, x)
f = Grenander(x)
plot(f)
title("Grenander 3 Ways")
g = medde(x, 1000, lambda = -1, alpha = 1, Dorder = 0)
lines(g$x, g$y, type = "s", col = 2)
h = Umix(sort(x))
lines(c(0,h$x), c(h$g[1],h$g), type = "S", col = 3)
legend("topright", c("LCM", "medde", "NPMLE"), lty = 1, col = 1:3)
# Further comparison with fdrtool version
#e = ecdf(x)
#h = grenander(e)
#lines(h$x.knots, h$f.knots, type = "s")
#s = sum( h$f.knots[-length(h$f.knots)]*diff(h$x.knots) )
