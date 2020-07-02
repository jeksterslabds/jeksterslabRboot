#' Alpha to Probabilities
#'
#' Calculates the probabilities of confidence limits
#' associated with specified alpha levels
#' for a two-tailed test.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param alpha Numeric vector.
#'   Alpha level.
#'   By default,
#'   `alpha` is set to conventional
#'   alpha levels
#'   `alpha = c(0.001, 0.01, 0.05)`.
#' @examples
#' # vector
#' alpha2prob(alpha = c(0.001, 0.01, 0.05))
#' # single numeric value
#' alpha2prob(alpha = 0.05)
#' @return Returns
#'   probabilities associated with specified alpha levels.
#'   The results are sorted from smallest to biggest.
#' @family alpha functions
#' @keywords alpha
#' @export
alpha2prob <- function(alpha = c(
                         0.001,
                         0.01,
                         0.05
                       )) {
  alpha <- sort(alpha)
  prob_ll <- alpha / 2
  prob_ul <- rev(1 - prob_ll)
  c(prob_ll, prob_ul)
}

#' Confidence Intervals to Probabilities
#'
#' Calculates the probabilities of confidence limits
#' associated with specified confidence intervals
#' for a two-tailed test.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param ci Numeric vector.
#'   Confidence interval.
#'   By default,
#'   `ci` is set to conventional
#'   confidence intervals
#'   `ci = c(0.999, 0.99, 0.95)`.
#' @examples
#' # vector
#' ci2prob(ci = c(0.999, 0.99, 0.95))
#' # single numeric value
#' ci2prob(ci = 0.95)
#' @return Returns
#'   probabilities associated with specified confidence intervals.
#'   The results are sorted from smallest to biggest.
#' @family alpha functions
#' @keywords alpha
#' @export
ci2prob <- function(ci = c(
                      0.999,
                      0.99,
                      0.95
                    )) {
  alpha2prob(alpha = 1 - ci)
}
