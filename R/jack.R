#' Jackknife
#'
#' Generates jackknife samples.
#'
#' \deqn{
#'   \mathbf{x}_i
#'   =
#'   \left\{
#'     x_1,
#'     x_2,
#'     x_3,
#'     \dots,
#'     x_n
#'   \right\}
#' }
#' for \eqn{i = \left\{ 1, 2, 3, \dots, n \right\}}
#' the \eqn{i}th jackknife sample
#' consists of the original sample data
#' with the \eqn{i}th observation removed.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param data Vector, matrix, or data frame.
#'   Sample data.
#' @inheritParams nb
#' @return Returns a list of jackknife samples
#'   of length `n`,
#'   where `n` is the sample size of `data`.
#'   Each item in the list,
#'   will have a sample size of
#'   `n - 1`.
#' @examples
#' n <- 5
#' # vector
#' data_vector <- rnorm(n = n)
#' jack(
#'   data = data_vector
#' )
#' # matrix
#' X <- rnorm(n = n)
#' Y <- rnorm(n = n)
#' Z <- rnorm(n = n)
#' data_matrix <- cbind(X, Y, Z)
#' jack(
#'   data = data_matrix
#' )
#' # data frame
#' data_dataframe <- data.frame(X, Y, Z)
#' jack(
#'   data = data_dataframe
#' )
#' @family jackknife functions
#' @keywords jackknife
#' @references
#'   [Wikipedia: Jackknife resampling](https://en.wikipedia.org/wiki/Jackknife_resampling)
#' @export
jack <- function(data,
                 par = FALSE,
                 ncores = NULL) {
  if (is.vector(data)) {
    n <- length(data)
  }
  if (is.matrix(data) | is.data.frame(data)) {
    n <- nrow(data)
  }
  foo <- function(i, data) {
    if (is.data.frame(data) | is.matrix(data)) {
      return(data[-i, ])
    }
    if (is.vector(data)) {
      return(data[-i])
    }
  }
  util_lapply(
    FUN = foo,
    args = list(
      i = 1:n,
      data = data
    ),
    par = par,
    ncores = ncores
  )
}

#' Jackknife Estimates
#'
#' Calculates jackknife estimates.
#'
#' The jackknife estimate of bias is given by
#' \deqn{
#'   \hat{\mathrm{bias}}_{\mathrm{jack}}
#'   \left(
#'     \theta
#'   \right)
#'   =
#'   \left(
#'     n - 1
#'   \right)
#'   \left(
#'     \hat{\theta}_{
#'       \left(
#'         \cdot
#'       \right)
#'     }
#'     -
#'     \hat{\theta}
#'   \right)
#' }
#'
#' where
#'
#' \deqn{
#'   \hat{\theta}_{
#'     \left(
#'       \cdot
#'     \right)
#'   }
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \hat{\theta}_{
#'     \left(
#'       i
#'     \right)
#'   } .
#' }
#'
#' The jackknife estimate of standard error
#' is given by
#'
#' \deqn{
#'   \hat{\mathrm{se}}_{\mathrm{jack}}
#'   \left(
#'     \hat{\theta}
#'   \right)
#'   =
#'   \sqrt{
#'     \frac{n - 1}{n}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       \hat{\theta}_{\left( i \right)}
#'       -
#'       \hat{\theta}_{\left( \cdot \right)}
#'     \right)^2
#'   } .
#' }
#'
#' The bias-corrected jackknife estimate
#' is given by
#'
#' \deqn{
#'   \hat{\theta}_{\mathrm{jack}}
#'   =
#'   \hat{\theta}
#'   -
#'   \hat{\mathrm{se}}_{\mathrm{jack}}
#'   \left(
#'     \hat{\theta}
#'   \right)
#'   =
#'   n
#'   \hat{\theta}
#'   -
#'   \left(
#'     n - 1
#'   \right)
#'   \hat{\theta}_{
#'     \left(
#'       \cdot
#'     \right)
#'   } .
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_star Numeric vector.
#'   Jackknife sampling distribution,
#'   that is,
#'   the sampling distribution of `thetahat`
#'   estimated for each `i` jackknife sample.
#'   \eqn{
#'     \hat{\theta}_{\left( 1 \right)},
#'     \hat{\theta}_{\left( 2 \right)},
#'     \hat{\theta}_{\left( 3 \right)},
#'     \dots,
#'     \hat{\theta}_{\left( n \right)}
#'   } .
#' @param thetahat Numeric.
#'   Parameter estimate
#'   (\eqn{\hat{\theta}})
#'   from the original sample data.
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{mean}{Mean of `thetahat_star` (\eqn{\hat{\theta}_{\left( \cdot \right) }}).}
#'     \item{bias}{Jackknife estimate of bias (\eqn{\hat{\mathrm{bias}}_{\mathrm{jack}} \left( \theta \right)}).}
#'     \item{se}{Jackknife estimate of standard error (\eqn{\hat{\mathrm{se}}_{\mathrm{jack}} \left( \hat{\theta} \right)}).}
#'     \item{thetahat_jack}{Bias-corrected jackknife estimate (\eqn{\hat{\theta}_{\mathrm{jack}}}).}
#'   }
#' @family jackknife functions
#' @keywords jackknife
#' @inherit jack references
#' @examples
#' n <- 100
#' # vector
#' x <- rnorm(n = n)
#' thetahat <- mean(x)
#' jack_samples <- jack(
#'   data = x
#' )
#' thetahat_star <- sapply(
#'   X = jack_samples,
#'   FUN = mean
#' )
#' hist(thetahat_star)
#' jack_hat(
#'   thetahat_star = thetahat_star,
#'   thetahat = thetahat
#' )
#' @export
jack_hat <- function(thetahat_star,
                     thetahat) {
  n <- length(thetahat_star)
  mean_thetahat <- mean(thetahat_star)
  bias <- (n - 1) * (mean_thetahat - thetahat)
  se <- sqrt(((n - 1) / n) * sum(thetahat_star - mean_thetahat)^2)
  thetahat_jack <- mean_thetahat - bias
  c(
    mean = mean_thetahat,
    bias = bias,
    se = se,
    thetahat_jack = thetahat_jack
  )
}
