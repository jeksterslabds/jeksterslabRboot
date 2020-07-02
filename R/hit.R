#' Confidence Interval - Theta Hit
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_lo Numeric.
#'   Lower limit of the estimated confidence interval
#'   (\eqn{\hat{\theta}_{\mathrm{lo}}}).
#' @param theta Numeric.
#'   Population parameter
#'   (\eqn{\theta}).
#' @param thetahat_up Numeric.
#'   Upper limit of the estimated confidence interval
#'   (\eqn{\hat{\theta}_{\mathrm{up}}}).
#' @return Returns
#'   `TRUE` if `theta` (\eqn{\theta}) is between the interval
#'   `thetahat_lo`
#'   (\eqn{\hat{\theta}_{\mathrm{lo}}})
#'   to
#'   `thetahat_up`
#'   (\eqn{\hat{\theta}_{\mathrm{up}}}).
#'   Returns
#'   `FALSE` if `theta` (\eqn{\theta}) is outside the interval
#'   `thetahat_lo`
#'   (\eqn{\hat{\theta}_{\mathrm{lo}}})
#'   to
#'   `thetahat_up`
#'   (\eqn{\hat{\theta}_{\mathrm{up}}}).
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
#'   (\eqn{\hat{\theta}_{\mathrm{lo}}})
#'   to
#'   `thetahat_up`
#'   (\eqn{\hat{\theta}_{\mathrm{up}}}).
#'   Returns
#'   `FALSE` if zero is outside the interval
#'   `thetahat_lo`
#'   (\eqn{\hat{\theta}_{\mathrm{lo}}})
#'   to
#'   `thetahat_up`
#'   (\eqn{\hat{\theta}_{\mathrm{up}}}).
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
