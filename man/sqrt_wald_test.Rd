% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci.R
\name{sqrt_wald_test}
\alias{sqrt_wald_test}
\title{Square Root of the Wald Test Statistic for a Single Parameter}
\usage{
sqrt_wald_test(thetahat, se_thetahat, theta_null = 0, distribution = "z", df)
}
\arguments{
\item{thetahat}{Numeric.
Parameter estimate
(\eqn{\hat{\theta}}).}

\item{se_thetahat}{Numeric.
Standard error of \code{thetahat}
(\eqn{\mathrm{se} \left( \hat{\theta} \right)}).}

\item{theta_null}{Numeric.
Hypothesized value of \code{theta}
(\eqn{\theta_{0}}).
Set to zero by default.}

\item{distribution}{Character string.
\code{distribution = "z"} for the standard normal distribution.
\code{distribution = "t"} for the t distribution.}

\item{df}{Numeric.
Degrees of freedom (df) if \code{distribution = "t"}.
Ignored if \code{distribution = "z"}.}
}
\value{
Returns a vector with the following elements:
\describe{
\item{statistic}{Square root of Wald test statistic.}
\item{p}{p-value.}
}
}
\description{
Calculates the square root of the Wald test statistic for a single parameter.
}
\details{
The square root of the Wald test statistic for a single parameter is calculated
as follows:
\deqn{
    \sqrt{W}
    =
    \frac{
      \hat{\theta} - \theta_0
    }
    {
      \mathrm{se}
      \left(
        \hat{\theta}
      \right)
    }
  }
The associated \code{p}-value
from the \code{z} or \code{t} distribution is calculated
using \code{\link[=pnorm]{pnorm()}} or \code{\link[=pt]{pt()}}.
}
\references{
\href{https://en.wikipedia.org/wiki/Wald_test}{Wikipedia: Wald test}
}
\author{
Ivan Jacob Agaloos Pesigan
}