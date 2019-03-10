\name{varwci}
\alias{varwci}
\title{
  varwci
}
\description{
  Surround the univariate variance estimator of the function var with a
  confidence interval, not assuming normality
}
\usage{
varwci(x, conf.level=0.95)
}
\arguments{
  \item{x}{A one-dimensional numeric vector}
  \item{conf.level}{The confidence level for the confidence
    interval. Defaults to 0.95}
}


\value{
  Returns a vector with two entries: the lower and the upper bound of the
  confidence interval, and the following attributes:
  \describe{
    \item{point.estimator}{The usual sample variance at the center of the interval}
    \item{conf.level}{The confidence level used}
    \item{var.SampleVariance}{The estimated variance of the sample variance}
  }
}

\references{
  http://dx.doi.org/10.1080/15598608.2016.1158675 and https://mathiasfuchs.de/b3.html
}

\author{
Mathias Fuchs
}

\note{
The underlying theory is that of U-statistics. See Hoeffding 1948.
}

\section{Warning }{
On very small sample sizes, the result is NA because there is insufficient information on the variance estimation
}

\examples{
##
## Example: throwing a dice
## 

# throw a dice 100 times
s <- sample(6, 100, replace=TRUE)

# the standard point estimator for the variance
print(var(s))

# contains the true value 2.9166 with a probability of 95 percent.
print(varwci(s))

##
## Check the coverage probability of the confidence interval
##

                                        # True quantities that do not depend on n
trueMeanOfDice <- mean(1:6)
trueVarianceOfDice <- mean((1:6)^2) - trueMeanOfDice^2

## see package description for more details
                                        # number of times we draw a
                                        # sample and compute a confidence interval
N <- 1e4
trueValueCovered <- rep(NA, N)
for (i in 1:N) {
    if (i \%\% 1e3 == 0) print(i)
                                        # throw a dice 100 times
    x <- sample(6, 100, replace=TRUE)
                                        # compute our confidence interval
    ci <- varwci(x)
                                        # We know that the true variance
                                        # of the dice is 91/6 - 49/4 = 2.916666...
                                        # did the confidence interval contain the correct value?
    trueValueCovered[i] <- (trueVarianceOfDice > ci[1] && trueVarianceOfDice < ci[2])
}

                                        # Result of simulation study: should be close to 0.95
print(mean(trueValueCovered))
}

\keyword{variance estimation}
\keyword{confidence interval}
