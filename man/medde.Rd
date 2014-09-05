\name{medde}
\alias{medde}
\alias{qmedde}
\alias{rmedde}
\alias{mesh1}
\alias{plot.medde}
\title{Maximum Entropy [De]Regularized Density Estimation}
\description{Density estimation based on maximum entropy methods}
\usage{ medde(x, v = 300, lambda, alpha = 1, Dorder = 1, rtol = 1e-06, 
	verb = 0, control = NULL)
	qmedde(p,medde)
	rmedde(n,medde,smooth = TRUE)
}
\arguments{
  \item{x}{Data:  either univariate or bivariate (not yet implemented in Rmosek)} 
  \item{v}{Undata:  either univariate or bivariate, by default there is an
	equally spaced grid of 300 values} 
  \item{lambda}{ total variation penalty parameter, if lambda is in [-1,0], 
	a concavity constraint is imposed. If lambda is in (-oo, -1)
	a convexity constraint on .5 x^2 + log f is imposed. See Koenker and Mizera 
	(2013) for further details on this last option, and Koenker and Mizera (2010)
	for further details on the concavity constrained options.}
  \item{alpha}{ Renyi entropy parameter characterizing fidelity criterion }
  \item{Dorder}{ Order of the derivative operator for the penalty} 
  \item{rtol}{ Convergence tolerance for Mosek algorithm}
  \item{verb}{ Parameter controlling verbosity of solution, 0 for silent, 5 
	gives rather detailed iteration log.} 
  \item{control}{ Mosek control list see KWDual documentation}
  \item{n}{ number of observations desired in calls to rmedde}
  \item{p}{ probabilities to evaluate in calls to qmedde}
  \item{medde}{ fitted medde object for calls in qmedde and rmedde}
  \item{smooth}{ option to draw random meddes from the smoothed density
	estimate, as described in Chen and Samworth (2013)}
}
\details{
See the references for further details. And also Mosek "Manuals".
The acronym has a nice connection to 

        http://www.urbandictionary.com/define.php?term=medder

        Term used in Bahamian dialect, mostly on the Family Islands like 
        Eleuthera and Cat Island meaning "mess with" "get involved" 
        "get entangled" "fool around" "bother"

        "I don't like to medder up with all kinda people"

        "Don't medder with people (chirren)"

        "Why you think she medderin up in their business."

}
\value{An object of class "Medde" with components
\item{x}{points of evaluation on the domain of the density}
\item{y}{estimated function values at the evaluation points x}
\item{phi}{n by 2 matrix of (sorted) original data and estimated log 
	density values these data points}
\item{logLik}{log likelihood value}
\item{status}{exit status from Mosek}
}
\references{a 
Chen, Y. and R.J. Samworth, (2013) "Smoothed log-concave maximum 
likelihood estimation with applications", \emph{Statistica Sinica}, 23, 1373--1398.

Koenker, R and I. Mizera, (2007) 
``Density Estimation by Total Variation Regularization,'' 
\emph{Advances in Statistical Modeling and Inference:  Essays in Honor of
Kjell Doksum}, V.N. Nair (ed.), 613-634.

Koenker, R and I. Mizera, (2006)
``The alter egos of the regularized maximum likelihood density estimators:
deregularized maximum-entropy, Shannon, Renyi, Simpson, Gini, and stretched
strings,'' \emph{ Proceedings of the  7th Prague Symposium on
Asymptotic Statistics}.

Koenker, R and I. Mizera, (2010) ``Quasi-Concave Density Estimation'' 
\emph{Annals of Statistics}, 38, 2998-3027.

Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
Compound Decisions, and Empirical Bayes Rules,'' JASA, 109, 674--685. 

Koenker, R and I. Mizera, (2014) ``Convex Optimization in R.'', 
\emph{Journal of Statistical Software}, 60, 1-23.
}

\author{Roger Koenker}

\seealso{ This function is based on an earlier function of the same name
in the deprecated package MeddeR that was based on an R-Matlab interface.
}
\examples{
#Maximum Likelihood Estimation of a Log-Concave Density
set.seed(1968)
x <- rgamma(50,10)
m <- medde(x, v = 50, lambda = -.5, verb = 5)
plot(m)
lines(m$x,dgamma(m$x,10),col = 2)
points(x,m$g,cex = 0.5)
rug(x)
title(paste("log likelihood = ", round(m$logLik,2)))
legend(14,.12,c("ghat","true"),lty = 1, col = 1:2)

#Maximum Likelihood Estimation of a Gamma Density with TV constraint
set.seed(1968)
x <- rgamma(50,5)
f <- medde(x, v = 50, lambda = 0.005, verb = 5)
plot(f)
lines(f$x,dgamma(f$x,5),col = 2)
legend(10,.15,c("ghat","true"),lty = 1, col = 1:2)

}
\keyword{nonparametric}
