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
len <- function(lo,
                up) {
  up - lo
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
#'   \eqn{\left( \hat{\theta} \right)}.
#' @inheritParams theta_hit
#' @family confidence interval evaluation functions
#' @keywords confidence interval
#' @export
shape <- function(lo,
                  thetahat,
                  up) {
  (up - thetahat) / (thetahat - lo)
}

#' Confidence Interval Evaluation
#'
#' Evaluates confidence interval using
#' [`zero_hit()`],
#' [`theta_hit()`],
#' [`len()`],
#' and
#' [`shape()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param ci Vector.
#'   Confidence intervals sorted from smallest to largest.
#'   The length should be even.
#'   The first and the last element correspond to the widest confidence interval.
#'   The second and the second to the last element correspond to the second widest confidence interval.
#'   And so on and so forth.
#' @param thetahat Numeric.
#'   Parameter estimate
#'   \eqn{\left( \hat{\theta} \right)}.
#' @param theta Numeric.
#'   Population parameter
#'   \eqn{\left( \theta \right)}.
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
#' @examples
#' ci <- c(
#'   98.04786,
#'   98.38773,
#'   98.68060,
#'   100.5447,
#'   100.8375,
#'   101.1774
#' )
#' thetahat <- 99.6126336
#' theta <- 100
#' label <- c(
#'   0.001,
#'   0.01,
#'   0.05
#' )
#' ci_eval(
#'   ci = ci,
#'   thetahat = thetahat,
#'   theta = theta,
#'   label = label
#' )
#' @export
ci_eval <- function(ci,
                    thetahat,
                    theta = 0,
                    label = NULL) {
  ci <- sort(ci)
  len_ci <- length(ci)
  half_ci <- len_ci / 2
  if ((len_ci %% 2) != 0) {
    stop(
      "Length of ci should be even."
    )
  }
  if (length(label) != half_ci) {
    message(
      "Length of label is incompatible. Defaulting to label = 1:(length(ci)/2)."
    )
    label <- NULL
  }
  if (is.null(label)) {
    label <- 1:half_ci
  }
  shape_vector <- len_vector <- zero_hit_vector <- theta_hit_vector <- rep(x = NA, times = half_ci)
  for (i in 1:half_ci) {
    lo <- ci[i]
    up <- ci[1 + len_ci - i]
    theta_hit_vector[i] <- theta_hit(
      lo = lo,
      theta = theta,
      up = up
    )
    zero_hit_vector[i] <- zero_hit(
      lo = lo,
      up = up
    )
    len_vector[i] <- len(
      lo = lo,
      up = up
    )
    shape_vector[i] <- shape(
      lo = lo,
      thetahat = thetahat,
      up = up
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
