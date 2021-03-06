% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alpha2prob.R
\name{ci2prob}
\alias{ci2prob}
\title{Confidence Intervals to Probabilities}
\usage{
ci2prob(ci = c(0.999, 0.99, 0.95))
}
\arguments{
\item{ci}{Numeric vector.
Confidence interval.
By default,
\code{ci} is set to conventional
confidence intervals
\code{ci = c(0.999, 0.99, 0.95)}.}
}
\value{
Returns
probabilities associated with the specified confidence interval/s.
The results are sorted from smallest to largest.
}
\description{
Calculates the cumulative probabilities of confidence limits
associated with the specified confidence interval/s
for a two-tailed test.
}
\examples{
# vector
ci2prob(ci = c(0.999, 0.99, 0.95))
# single numeric value
ci2prob(ci = 0.95)
}
\references{
\href{https://en.wikipedia.org/wiki/Statistical_significance}{Wikipedia: Statistical significance}

\href{https://en.wikipedia.org/wiki/Confidence_interval}{Wikipedia: Confidence interval}
}
\seealso{
Other alpha functions: 
\code{\link{alpha2crit}()},
\code{\link{alpha2prob}()},
\code{\link{ci2crit}()},
\code{\link{nhstplot}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{alpha functions}
\keyword{alpha}
