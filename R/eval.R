#' Confidence Interval - Length
#'
#' Calculates confidence interval length.
#'
#' The confidence interval length is given by
#' \deqn{
#'   \mathrm{confidence \ interval \ length}
#'   =
#'   \hat{\theta}_{\mathrm{up}}
#'   -
#'   \hat{\theta}_{\mathrm{lo}}
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams theta_hit
#' @family confidence interval evaluation functions
#' @keywords confidence interval
#' @export
len <- function(thetahat_lo,
                thetahat_up) {
  thetahat_up - thetahat_lo
}

#' Confidence Interval - Shape
#'
#' Calculates confidence interval shape.
#'
#' The confidence interval shape is given by
#' \deqn{
#'   \mathrm{confidence \ interval \ shape}
#'   =
#'   \frac{
#'     \hat{\theta}_{\mathrm{up}} - \hat{\theta}
#'   }
#'   {
#'     \hat{\theta} - \hat{\theta}_{\mathrm{lo}}
#'   }
#' }
#'
#' The shape measures the asymmetry of the confidence interval
#' around the point estimate \eqn{\hat{\theta}}.
#' Shape \eqn{> 1.00} is indicative of greater distance between
#' \eqn{\hat{\theta}_{\mathrm{up}}} and \eqn{\hat{\theta}}
#' than
#' \eqn{\hat{\theta}} and \eqn{\hat{\theta}_{\mathrm{lo}}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat Numeric.
#'   Parameter estimate
#'   (\eqn{\hat{\theta}}).
#' @inheritParams theta_hit
#' @family confidence interval evaluation functions
#' @keywords confidence interval
#' @export
shape <- function(thetahat_lo,
                  thetahat,
                  thetahat_up) {
  (thetahat_up - thetahat) / (thetahat - thetahat_lo)
}
