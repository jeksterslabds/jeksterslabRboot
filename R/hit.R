#' Confidence Interval - Theta Hit
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_lo Numeric.
#'   Lower limit of the estimated confidence interval
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}.
#' @param theta Numeric.
#'   Population parameter
#'   \eqn{\left( \theta \right)}.
#' @param thetahat_up Numeric.
#'   Upper limit of the estimated confidence interval
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#' @return Returns
#'   `TRUE` if `theta` \eqn{\left( \theta \right)} is between the interval
#'   `thetahat_lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `thetahat_up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#'   Returns
#'   `FALSE` if `theta` \eqn{\left( \theta \right)} is outside the interval
#'   `thetahat_lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `thetahat_up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#' @examples
#' # FALSE
#' theta_hit(thetahat_lo = 1, theta = 0, thetahat_up = 2)
#' # TRUE
#' theta_hit(thetahat_lo = -1, theta = 0, thetahat_up = 1)
#' @family confidence interval hit functions
#' @keywords confidence interval
#' @export
theta_hit <- function(thetahat_lo,
                      theta,
                      thetahat_up) {
  thetahat_lo < theta & theta < thetahat_up
}

#' Confidence Interval - Zero Hit
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams theta_hit
#' @return Returns
#'   `TRUE` if zero is between the interval
#'   `thetahat_lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `thetahat_up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#'   Returns
#'   `FALSE` if zero is outside the interval
#'   `thetahat_lo`
#'   \eqn{\left( \hat{\theta}_{\mathrm{lo}} \right)}
#'   to
#'   `thetahat_up`
#'   \eqn{\left( \hat{\theta}_{\mathrm{up}} \right)}.
#' @examples
#' # FALSE
#' zero_hit(thetahat_lo = 1, thetahat_up = 2)
#' # TRUE
#' zero_hit(thetahat_lo = -1, thetahat_up = 1)
#' @family confidence interval hit functions
#' @keywords confidence interval
#' @export
zero_hit <- function(thetahat_lo,
                     thetahat_up) {
  theta_hit(
    thetahat_lo,
    theta = 0,
    thetahat_up
  )
}
