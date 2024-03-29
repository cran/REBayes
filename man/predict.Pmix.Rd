% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.Pmix.R
\name{predict.Pmix}
\alias{predict.Pmix}
\title{Predict Method for Pmix}
\usage{
\method{predict}{Pmix}(object, newdata, Loss = 2, newexposure = NULL, ...)
}
\arguments{
\item{object}{fitted object of class "Pmix"}

\item{newdata}{Values at which prediction is desired}

\item{Loss}{Loss function used to generate prediction.  Currently supported values:
2 to get mean predictions, 1 to get harmonic mean predictions, 0 to get modal predictions
or any tau in (0,1) to get tau-th quantile predictions.  The posterior harmonic mean is
the Bayes rule for quadratic loss weighted by variances as in Clevenson and Zidek (1975).}

\item{newexposure}{exposure values for the predictions}

\item{...}{optional arguments to predict}
}
\value{
A vector of predictions
}
\description{
Predict Method for Poisson Mixtures
}
\details{
The predict method for \code{Pmix} objects will compute means, quantiles or
modes of the posterior according to the \code{Loss} argument.  Typically,
\code{newdata} would be passed to \code{predict}
}
\references{
Clevenson, M. L. and Zidek, J. V. 1975. Simultaneous Estimation of the
Means of Independent Poisson Laws, Journal of the American Statistical Association,
70, 698-705.
}
\author{
Jiaying Gu and Roger Koenker
}
\keyword{nonparametric}
