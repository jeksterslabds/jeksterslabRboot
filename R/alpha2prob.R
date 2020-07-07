#' Alpha to Probabilities
#'
#' Calculates the cumulative probabilities of confidence limits
#' associated with the specified significance level/s
#' \eqn{\left( \alpha \right)}
#' for a two-tailed test.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family alpha functions
#' @keywords alpha
#' @param alpha Numeric vector.
#'   Significance level
#'   \eqn{\left( \alpha \right)} .
#'   By default,
#'   `alpha` is set to conventional
#'   significance levels
#'   `alpha = c(0.001, 0.01, 0.05)`.
#' @return Returns
#'   probabilities associated with
#'   the specified significance level/s
#'   \eqn{\left( \alpha \right)}
#'   for a two-tailed test.
#'   The results are sorted from smallest to largest.
#' @examples
#' # vector
#' alpha2prob(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   )
#' )
#' # single numeric value
#' alpha2prob(alpha = 0.05)
#' @references
#' [Wikipedia: Statistical significance](https://en.wikipedia.org/wiki/Statistical_significance)
#'
#' [Wikipedia: Confidence interval](https://en.wikipedia.org/wiki/Confidence_interval)
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
#' Calculates the cumulative probabilities of confidence limits
#' associated with the specified confidence interval/s
#' for a two-tailed test.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family alpha functions
#' @keywords alpha
#' @inherit alpha2prob references
#' @param ci Numeric vector.
#'   Confidence interval.
#'   By default,
#'   `ci` is set to conventional
#'   confidence intervals
#'   `ci = c(0.999, 0.99, 0.95)`.
#' @return Returns
#'   probabilities associated with the specified confidence interval/s.
#'   The results are sorted from smallest to largest.
#' @examples
#' # vector
#' ci2prob(ci = c(0.999, 0.99, 0.95))
#' # single numeric value
#' ci2prob(ci = 0.95)
#' @export
ci2prob <- function(ci = c(
                      0.999,
                      0.99,
                      0.95
                    )) {
  alpha2prob(alpha = 1 - ci)
}
