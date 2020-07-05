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
#' #################################
#' # vector
#' #################################
#' x <- rnorm(n = n)
#' x_star <- jack(
#'   data = x
#' )
#' str(x_star)
#' #################################
#' # matrix
#' #################################
#' x1 <- rnorm(n = n)
#' x2 <- rnorm(n = n)
#' x3 <- rnorm(n = n)
#' X <- cbind(x1, x2, x3)
#' X_star <- jack(
#'   data = X
#' )
#' str(X_star)
#' #################################
#' # data frame
#' #################################
#' X <- as.data.frame(X)
#' X_star <- jack(
#'   data = X
#' )
#' str(X_star)
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
#'   \widehat{\mathrm{bias}}_{\mathrm{jack}}
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
#'   \widehat{\mathrm{se}}_{\mathrm{jack}}
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
#'   \hat{\mathrm{bias}}_{\mathrm{jack}}
#'   \left(
#'     \theta
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
#' Pseudo-values can be computed using
#'
#' \deqn{
#'   \tilde{\theta}_{i}
#'   =
#'   n \hat{\theta}
#'   -
#'   \left(
#'     n - 1
#'   \right)
#'   \hat{\theta}_{\left( i \right)} .
#' }
#'
#' The standard error can be estimated using the pseudo-values
#'
#' \deqn{
#'   \widehat{\mathrm{se}}_{\mathrm{jack}}
#'   \left(
#'     \tilde{\theta}
#'   \right)
#'   =
#'   \sqrt{
#'     \sum_{i = 1}^{n}
#'     \frac{
#'       \left(
#'         \tilde{\theta}_{i}
#'         -
#'         \tilde{\theta}
#'       \right)^2
#'     }
#'     {
#'       \left(
#'         n
#'         -
#'         1
#'       \right)
#'       n
#'     }
#'   }
#' }
#'
#' where
#'
#' \deqn{
#'   \tilde{\theta}
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \tilde{\theta}_{i} .
#' }
#'
#' An interval can be generated using
#'
#' \deqn{
#'   \tilde{\theta}
#'   \pm
#'   t_{\frac{\alpha}{2}}
#'   \times
#'   \widehat{\mathrm{se}}_{\mathrm{jack}}
#'     \left(
#'       \tilde{\theta}
#'     \right)
#' }
#' with degrees of freedom \eqn{\nu = n - 1}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_star_jack Numeric vector.
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
#'   \eqn{\left( \hat{\theta} \right)}
#'   from the original sample data.
#' @inheritParams wald
#' @return Returns a list with the following elements:
#'   \describe{
#'     \item{hat}{Jackknife estimates.}
#'     \item{ps}{Pseudo-values.}
#'     \item{ci}{Confidence intervals using pseudo-values.}
#'   }
#' The first list element `hat` contains the following:
#'   \describe{
#'     \item{mean}{Mean of `thetahat_star_jack` \eqn{\left( \hat{\theta}_{\left( \cdot \right) } \right)}.}
#'     \item{bias}{Jackknife estimate of bias \eqn{\left( \widehat{\mathrm{bias}}_{\mathrm{jack}} \left( \theta \right) \right)}.}
#'     \item{se}{Jackknife estimate of standard error \eqn{\left( \widehat{\mathrm{se}}_{\mathrm{jack}} \left( \hat{\theta} \right) \right)}.}
#'     \item{thetahat_jack}{Bias-corrected jackknife estimate \eqn{\left( \hat{\theta}_{\mathrm{jack}} \right)}.}
#'   }
#' @family jackknife functions
#' @keywords jackknife
#' @inherit jack references
#' @examples
#' n <- 100
#' x <- rnorm(n = n)
#' thetahat <- mean(x)
#' x_star <- jack(
#'   data = x
#' )
#' thetahat_star_jack <- sapply(
#'   X = x_star,
#'   FUN = mean
#' )
#' str(x_star)
#' hist(
#'   thetahat_star_jack,
#'   main = expression(
#'     paste(
#'       "Histogram of ",
#'       hat(theta),
#'       "*"
#'     )
#'   ),
#'   xlab = expression(
#'     paste(
#'       hat(theta),
#'       "*"
#'     )
#'   )
#' )
#' jack_hat(
#'   thetahat_star_jack = thetahat_star_jack,
#'   thetahat = thetahat
#' )
#' @export
jack_hat <- function(thetahat_star_jack,
                     thetahat,
                     alpha = c(
                       0.001,
                       0.01,
                       0.05
                     ),
                     eval = FALSE,
                     theta = 0) {
  n <- length(thetahat_star_jack)
  mean_thetahat <- mean(thetahat_star_jack)
  bias <- (n - 1) * (mean_thetahat - thetahat)
  se <- sqrt(((n - 1) / n) * sum(thetahat_star_jack - mean_thetahat)^2)
  thetahat_jack <- mean_thetahat - bias
  # pseudo-values
  pseudo_values <- n * thetahat - (n - 1) * thetahat_star_jack
  mean_pseudo_values <- mean(pseudo_values)
  # research more on formula for se_pseudo_values
  se_pseudo_values <- sqrt(sum(((pseudo_values - mean_pseudo_values)^2) / ((n - 1) * n)))
  ci <- wald(
    thetahat = mean_pseudo_values,
    sehat_thetahat = se_pseudo_values,
    theta_null = 0,
    alpha = alpha,
    distribution = "t",
    df = n - 1,
    eval = eval,
    theta = theta
  )
  ci[1] <- NA
  ci[2] <- NA
  hat <- c(
    mean = mean_thetahat,
    bias = bias,
    se = se,
    thetahat_jack = thetahat_jack,
    mean_ps = mean_pseudo_values,
    se_ps = se_pseudo_values
  )
  list(
    hat = hat,
    ps = pseudo_values,
    ci = ci
  )
}
