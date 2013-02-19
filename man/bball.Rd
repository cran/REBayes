\name{bball}
\alias{bball}
\docType{data}
\title{
Baseball Batting Averages
}
\description{
U.S. Major League Batting Average Data:  2002-2011
}
\usage{data(bball)}
\format{
  A data frame with 10575 longitudinal observations on the following 8 variables.
  \describe{
    \item{\code{name}}{name of player}
    \item{\code{id}}{numeric (unique) identifier of the player}
    \item{\code{year}}{year played, see details}
    \item{\code{season}}{season played, see details}
    \item{\code{pitcher}}{indicator: 1 if pitcher}
    \item{\code{HA}}{transformed batting averages, see details}
    \item{\code{AB}}{at bats}
    \item{\code{H}}{hits}
  }
}
\details{
Data is aggregated into half seasons:  so season indicates whether the observation
is in the first or second half of the season of a given year.
Only players who have more than 10 at bats in any half season are included,
and only players who have more than three half seasons are represented.
The transformed batting average is arcsin(sqrt((H + 1/4)/(AB + 1/2))).
R programs to extract the data from the original sources are available on request.
}
\source{
ESPN Website: {http://espn.go.com/mlb/statistics}
}
\examples{
data(bball)
}
\keyword{datasets}
