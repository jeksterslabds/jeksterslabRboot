#' Alpha to Critical Values
#'
#' Calculates the \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#' critical value/s
#' of confidence limits
#' associated with the specified significance level/s
#' \eqn{\left( \alpha \right)} .
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family alpha functions
#' @keywords alpha
#' @inheritParams alpha2prob
#' @inherit alpha2prob references
#' @param dist Character string.
#'   `dist = "z"` for the standard normal distribution.
#'   `dist = "t"` for the t distribution.
#'   `dist = "F"` for the F distribution.
#'   `dist = "chisq"` for the chi-square distribution.
#' @param two.tailed Logical.
#'   If `TRUE`, two-tailed alpha.
#'   If `FALSE`, one-tailed alpha.
#'   Ignored if `dist = "F"` or `dist = "chisq"`
#'   as both tests are one-tailed using the right tail.
#' @param right.tail Logical.
#'   If `TRUE`, right tail (positive critical value).
#'   If `FALSE`, left tail (negative critical value).
#'   Ignored if `two.tailed = TRUE`.
#'   Ignored if `dist = "F"` or `dist = "chisq"`
#'   as both tests are one-tailed using the right tail.
#' @param ... Degrees of freedom.
#'   `df` for `dist = "t"` and `dist = "chisq"`.
#'   `df1` and `df2` for `dist = "F"`.
#' @return Returns
#'   \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#'   critical value/s
#'   associated with
#'   the specified significance level/s
#'   \eqn{\left( \alpha \right)}.
#'   The results are sorted from smallest to largest.
#' @examples
#' # z two-tailed
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   )
#' )
#' ## single numeric value
#' alpha2crit(alpha = 0.05)
#' # t two-tailed
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "t",
#'   df = 1000
#' )
#' ## single numeric value
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "t",
#'   df = 1000
#' )
#' # z one-tailed right
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "z",
#'   two.tailed = FALSE
#' )
#' # t one-tailed right
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "t",
#'   two.tailed = FALSE,
#'   df = 5
#' )
#' # F one-tailed
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "F",
#'   two.tailed = FALSE,
#'   df1 = 2,
#'   df2 = 2
#' )
#' # chisq one-tailed
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "chisq",
#'   two.tailed = FALSE,
#'   df = 1
#' )
#' # z one-tailed left
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "z",
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' # t one-tailed left
#' ## vector
#' alpha2crit(
#'   alpha = c(
#'     0.001,
#'     0.01,
#'     0.05
#'   ),
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @importFrom stats qf
#' @importFrom stats qchisq
#' @export
alpha2crit <- function(alpha = c(
                         0.001,
                         0.01,
                         0.05
                       ),
                       dist = "z",
                       two.tailed = TRUE,
                       right.tail = TRUE,
                       ...) {
  alpha <- sort(alpha)
  ci <- 1 - alpha
  tail <- 1 - ci
  if (dist == "F" | dist == "chisq") {
    two.tailed <- FALSE
    right.tail <- TRUE
  }
  if (right.tail) {
    lower.tail <- FALSE
  } else {
    lower.tail <- TRUE
  }
  if (two.tailed) {
    if (dist == "z") {
      return(
        qnorm(
          p = alpha2prob(
            alpha = alpha
          ),
          lower.tail = TRUE
        )
      )
    }
    if (dist == "t") {
      return(
        qt(
          p = alpha2prob(
            alpha = alpha
          ),
          lower.tail = TRUE,
          ...
        )
      )
    }
  } else {
    if (dist == "z") {
      return(
        qnorm(
          p = tail,
          lower.tail = lower.tail,
        )
      )
    }
    if (dist == "t") {
      return(
        qt(
          p = tail,
          lower.tail = lower.tail,
          ...
        )
      )
    }
    if (dist == "F") {
      return(
        qf(
          p = tail,
          lower.tail = lower.tail,
          ...
        )
      )
    }
    if (dist == "chisq") {
      return(
        qchisq(
          p = tail,
          lower.tail = lower.tail,
          ...
        )
      )
    }
  }
}

#' Confidence Intervals to Critical Values
#'
#' Calculates the \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#' critical value/s
#' of confidence limits
#' associated with the specified confidence interval/s.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family alpha functions
#' @keywords alpha
#' @inheritParams ci2prob
#' @inheritParams alpha2crit
#' @inherit alpha2prob references
#' @return Returns
#'   \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#'   critical value/s
#'   associated with
#'   the specified confidence interval/s.
#'   The results are sorted from smallest to largest.
#' @examples
#' # z two-tailed
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#' )
#' ## single numeric value
#' ci2crit(ci = 0.95)
#' # t two-tailed
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "t",
#'   df = 1000
#' )
#' ## single numeric value
#' ci2crit(
#'   ci = 0.95,
#'   dist = "t",
#'   df = 1000
#' )
#' # z one-tailed right
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "z",
#'   two.tailed = FALSE
#' )
#' # t one-tailed right
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "t",
#'   two.tailed = FALSE,
#'   df = 5
#' )
#' # F one-tailed
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "F",
#'   two.tailed = FALSE,
#'   df1 = 2,
#'   df2 = 2
#' )
#' # chisq one-tailed
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "chisq",
#'   two.tailed = FALSE,
#'   df = 1
#' )
#' # z one-tailed left
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "z",
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' # t one-tailed left
#' ## vector
#' ci2crit(
#'   ci = c(
#'     0.999,
#'     0.99,
#'     0.95
#'   ),
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' @export
ci2crit <- function(ci = c(
                      0.999,
                      0.99,
                      0.95
                    ),
                    dist = "z",
                    two.tailed = TRUE,
                    right.tail = TRUE,
                    ...) {
  alpha2crit(
    dist = dist,
    alpha = 1 - ci,
    two.tailed = two.tailed,
    right.tail = right.tail,
    ...
  )
}
