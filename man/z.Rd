% This is a boilerplate example.
% Delete this file if you are building documentation files with roxygen.
\name{z}
\alias{z}
\title{Calculate z-scores}
\usage{
z(x, mu, sigma)
}
\arguments{
\item{x}{Numeric vector.}
\item{mu}{Single value.
Population mean.}
\item{sigma}{Single value.
Population standard deviation.}
}
\value{
Returns a numeric vector of z-scores.
}
\description{
Standardizes a numeric vector
by making its mean 0 and variance 1.
}
\examples{
x <- rnorm(n = 10, mean = 100, sd = 15)
z(x = x, mu = 100, sigma = 15)
}
\author{
Ivan Jacob Agaloos Pesigan
}

