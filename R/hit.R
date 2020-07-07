#' Confidence Interval - Theta Hit
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family confidence interval hit functions
#' @keywords confidence interval
#' @param lo Numeric.
#'   Lower limit of the estimated confidence interval
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}.
#' @param theta Numeric.
#'   Population parameter
#'   \eqn{\left( \theta \right)}.
#' @param up Numeric.
#'   Upper limit of the estimated confidence interval
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#' @return Returns
#'   `TRUE` if `theta` \eqn{\left( \theta \right)} is between the interval
#'   `lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#'   Returns
#'   `FALSE` if `theta` \eqn{\left( \theta \right)} is outside the interval
#'   `lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#' @examples
#' # FALSE
#' theta_hit(lo = 1, theta = 0, up = 2)
#' # TRUE
#' theta_hit(lo = -1, theta = 0, up = 1)
#' @export
theta_hit <- function(lo,
                      theta,
                      up) {
  lo < theta & theta < up
}

#' Confidence Interval - Zero Hit
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family confidence interval hit functions
#' @keywords confidence interval
#' @inheritParams theta_hit
#' @return Returns
#'   `TRUE` if zero is between the interval
#'   `lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#'   Returns
#'   `FALSE` if zero is outside the interval
#'   `lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#' @examples
#' # FALSE
#' zero_hit(lo = 1, up = 2)
#' # TRUE
#' zero_hit(lo = -1, up = 1)
#' @export
zero_hit <- function(lo,
                     up) {
  theta_hit(
    lo,
    theta = 0,
    up
  )
}
