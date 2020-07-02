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

#' Confidence Interval Evaluation
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param ci Vector.
#'   Confidence intervals sorted from smallest to largest.
#'   Length should be even.
#'   The first and the last element correspond to the widest confidence interval.
#'   The second and the second to the last element correspond to the second widest confidence interval.
#'   And so on and so forth.
#' @param thetahat Numeric.
#'   Parameter estimate
#'   (\eqn{\hat{\theta}}).
#' @param theta Numeric.
#'   Population parameter
#'   (\eqn{\theta}).
#' @param label Vector.
#'   Vector used to label results.
#'   If not provided defaults to `label = 1:(length(ci)/2)`.
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{zero_hit_}{Logical. Tests if confidence interval contains zero.}
#'     \item{theta_hit_}{Logical. Tests if confidence interval contains theta.}
#'     \item{length_}{Length of confidence interval.}
#'     \item{shape_}{Shape of confidence interval.}
#'   }
#' @export
ci_eval <- function(ci,
                    thetahat,
                    theta = 0,
                    label = NULL) {
  if (is.null(label)) {
    label <- 1:half_ci
  }
  half_ci <- length(ci) / 2
  shape_vector <- len_vector <- zero_hit_vector <- theta_hit_vector <- rep(x = NA, times = half_ci)
  for (i in 1:half_ci) {
    thetahat_lo <- ci[i]
    thetahat_up <- ci[1 + length(ci) - i]
    theta_hit_vector[i] <- theta_hit(
      thetahat_lo = thetahat_lo,
      theta = theta,
      thetahat_up = thetahat_up
    )
    zero_hit_vector[i] <- zero_hit(
      thetahat_lo = thetahat_lo,
      thetahat_up = thetahat_up
    )
    len_vector[i] <- len(
      thetahat_lo = thetahat_lo,
      thetahat_up = thetahat_up
    )
    shape_vector[i] <- shape(
      thetahat_lo = thetahat_lo,
      thetahat = thetahat,
      thetahat_up = thetahat_up
    )
  }
  names(theta_hit_vector) <- paste0(
    "theta_hit_",
    label
  )
  names(zero_hit_vector) <- paste0(
    "zero_hit_",
    label
  )
  names(len_vector) <- paste0(
    "length_",
    label
  )
  names(shape_vector) <- paste0(
    "shape_",
    label
  )
  c(
    zero_hit_vector,
    theta_hit_vector,
    len_vector,
    shape_vector
  )
}
