#' Alpha to Probabilities
#'
#' Calculates the cumulative probabilities of confidence limits
#' associated with the specified significance level/s
#' \eqn{\left( \alpha \right)}
#' for a two-tailed test.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param alpha Numeric vector.
#'   Significance level
#'   \eqn{\left( \alpha \right)} .
#'   By default,
#'   `alpha` is set to conventional
#'   significance levels
#'   `alpha = c(0.001, 0.01, 0.05)`.
#' @examples
#' # vector
#' alpha2prob(alpha = c(0.001, 0.01, 0.05))
#' # single numeric value
#' alpha2prob(alpha = 0.05)
#' @return Returns
#'   probabilities associated with
#'   the specified significance level/s
#'   \eqn{\left( \alpha \right)}
#'   for a two-tailed test.
#'   The results are sorted from smallest to largest.
#' @family alpha functions
#' @keywords alpha
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

#' Alpha to Critical Values
#'
#' Calculates the \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#' critical value/s
#' of confidence limits
#' associated with the specified significance level/s
#' \eqn{\left( \alpha \right)} .
#'
#' @author Ivan Jacob Agaloos Pesigan
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
#' @inheritParams alpha2prob
#' @return Returns
#'   \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#'   critical value/s
#'   associated with
#'   the specified significance level/s
#'   \eqn{\left( \alpha \right)}.
#'   The results are sorted from smallest to largest.
#' @inherit alpha2prob references
#' @examples
#' # z two-tailed
#' ## vector
#' alpha2crit(alpha = c(0.001, 0.01, 0.05))
#' ## single numeric value
#' alpha2crit(alpha = 0.05)
#' # t two-tailed
#' ## vector
#' alpha2crit(alpha = c(0.001, 0.01, 0.05), dist = "t", df = 1000)
#' ## single numeric value
#' alpha2crit(alpha = 0.05, , dist = "t", df = 1000)
#' @family alpha functions
#' @keywords alpha
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @importFrom stats qf
#' @importFrom stats qchisq
#' @examples
#' # vector
#' ## two-tailed
#' ### z
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "z"
#' )
#' ### t
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "t",
#'   df = 5
#' )
#' ## one-tailed
#' ### right
#' #### z
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "z",
#'   two.tailed = FALSE
#' )
#' #### t
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   df = 5
#' )
#' #### F
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "F",
#'   two.tailed = FALSE,
#'   df1 = 5,
#'   df2 = 5
#' )
#' #### chi-square
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "chisq",
#'   two.tailed = FALSE,
#'   df = 5
#' )
#' ### left
#' #### z
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "z",
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' #### t
#' alpha2crit(
#'   alpha = 0.05,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
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

#' Confidence Intervals to Probabilities
#'
#' Calculates the cumulative probabilities of confidence limits
#' associated with the specified confidence interval/s
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
#'   probabilities associated with the specified confidence interval/s.
#'   The results are sorted from smallest to largest.
#' @family alpha functions
#' @keywords alpha
#' @inherit alpha2prob references
#' @export
ci2prob <- function(ci = c(
                      0.999,
                      0.99,
                      0.95
                    )) {
  alpha2prob(alpha = 1 - ci)
}

#' Confidence Intervals to Critical Values
#'
#' Calculates the \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#' critical value/s
#' of confidence limits
#' associated with the specified confidence interval/s.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ci2prob
#' @inheritParams alpha2crit
#' @examples
#' # vector
#' ## two-tailed
#' ### z
#' ci2crit(
#'   ci = 0.95,
#'   dist = "z"
#' )
#' ### t
#' ci2crit(
#'   ci = 0.95,
#'   dist = "t",
#'   df = 5
#' )
#' ## one-tailed
#' ### right
#' #### z
#' ci2crit(
#'   ci = 0.95,
#'   dist = "z",
#'   two.tailed = FALSE
#' )
#' #### t
#' ci2crit(
#'   ci = 0.95,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   df = 5
#' )
#' #### F
#' ci2crit(
#'   ci = 0.95,
#'   dist = "F",
#'   two.tailed = FALSE,
#'   df1 = 5,
#'   df2 = 5
#' )
#' #### chi-square
#' ci2crit(
#'   ci = 0.95,
#'   dist = "chisq",
#'   two.tailed = FALSE,
#'   df = 5
#' )
#' ### left
#' #### z
#' ci2crit(
#'   ci = 0.95,
#'   dist = "z",
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' #### t
#' ci2crit(
#'   ci = 0.95,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' @return Returns
#'   \eqn{z}, \eqn{t}, \eqn{\chi^2}, or \eqn{F}
#'   critical value/s
#'   associated with
#'   the specified confidence interval/s.
#'   The results are sorted from smallest to largest.
#' @inherit alpha2prob references
#' @family alpha functions
#' @keywords alpha
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

#' Alpha to Plot
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams alpha2crit
#' @param alpha Numeric.
#'   Significance level
#'   \eqn{\left( \alpha \right)} .
#' @param ci Logical.
#'   `ci = FALSE` by default.
#'   If `TRUE`, plot confidence interval.
#'   If `FALSE`, plot alpha level.
#' @param statistic Numeric.
#'   Test statistic.
#'   This is an optional argument.
#' @examples
#' # ci = FALSE
#' # two-tailed
#' alpha2plot(
#'   alpha = 0.05,
#'   dist = "z",
#'   statistic = 3
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   dist = "t",
#'   df = 5
#' )
#' # one-tailed
#' alpha2plot(
#'   alpha = 0.05,
#'   two.tailed = FALSE,
#'   right.tail = TRUE
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = TRUE,
#'   df = 5
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   dist = "F",
#'   df1 = 10,
#'   df2 = 3
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   dist = "chisq",
#'   df = 3
#' )
#' # ci = TRUE
#' # two-tailed
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "z",
#'   statistic = 3
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "t",
#'   df = 5
#' )
#' # one-tailed
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   two.tailed = FALSE,
#'   right.tail = TRUE
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = TRUE,
#'   df = 5
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "F",
#'   df1 = 10,
#'   df2 = 3
#' )
#' alpha2plot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "chisq",
#'   df = 3
#' )
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics polygon
#' @importFrom graphics abline
#' @importFrom stats dnorm
#' @importFrom stats dt
#' @importFrom stats dchisq
#' @export
alpha2plot <- function(alpha = 0.05,
                       ci = FALSE,
                       dist = "z",
                       two.tailed = TRUE,
                       right.tail = TRUE,
                       statistic = NULL,
                       ...) {
  if (length(alpha) > 1) {
    stop(
      "alpha should be a single numeric value."
    )
  }
  if (dist == "F" | dist == "chisq") {
    two.tailed <- FALSE
    right.tail <- TRUE
  }
  critical <- alpha2crit(
    alpha = alpha,
    dist = dist,
    two.tailed = two.tailed,
    right.tail = right.tail,
    ...
  )
  # two-tailed
  if (two.tailed) {
    # common to z and t
    lo <- critical[1]
    up <- critical[2]
    x_axis <- c(
      lo,
      0,
      up
    )
    x <- seq(
      from = lo - 2,
      to = up + 2,
      length = 1000
    )
    # z two-tailed
    if (dist == "z") {
      lo_density <- dnorm(
        x[x <= lo],
        ...
      )
      up_density <- dnorm(
        x[x >= up]
      )
      if (ci) {
        ci_center <- dnorm(
          x[x >= lo & x <= up]
        )
        subtitle <- bquote(
          P(.(lo) ~ "to" ~ .(up)) == .(1 - alpha)
        )
      } else {
        subtitle <- bquote(
          P(z <= .(lo)) + P(z >= .(up)) == .(alpha)
        )
      }
    }
    # t two-tailed
    if (dist == "t") {
      lo_density <- dt(
        x[x <= lo],
        ...
      )
      up_density <- dt(
        x[x >= up],
        ...
      )
      if (ci) {
        ci_center <- dt(
          x[x >= lo & x <= up],
          ...
        )
        subtitle <- bquote(
          P(.(lo) ~ "to" ~ .(up)) == .(1 - alpha)
        )
      } else {
        subtitle <- bquote(
          P(t <= .(lo)) + P(t >= .(up)) == .(alpha)
        )
      }
    }
  } else { # one-tailed
    # right
    if (right.tail) {
      up <- critical[1]
      x_axis <- c(
        0,
        up
      )
      x <- seq(
        from = -1 * up - 2,
        to = up + 2,
        length = 1000
      )
      if (dist == "F" | dist == "chisq") {
        x <- seq(
          from = 0,
          to = up + 2,
          length = 1000
        )
      }
      if (dist == "z") {
        up_density <- dnorm(
          x[x >= up]
        )
        if (ci) {
          ci_right_tail <- dnorm(
            x[x <= up]
          )
          subtitle <- bquote(
            P(z <= .(up)) == .(alpha)
          )
        } else {
          subtitle <- bquote(
            P(z >= .(up)) == .(alpha)
          )
        }
      }
      if (dist == "t") {
        up_density <- dt(
          x[x >= up],
          ...
        )
        if (ci) {
          ci_right_tail <- dt(
            x[x <= up],
            ...
          )
          subtitle <- bquote(
            P(z <= .(up)) == .(1 - alpha)
          )
        } else {
          subtitle <- bquote(
            P(t >= .(up)) == .(alpha)
          )
        }
      }
      if (dist == "F") {
        up_density <- df(
          x[x >= up],
          ...
        )
        if (ci) {
          ci_right_tail <- df(
            x[x <= up],
            ...
          )
          subtitle <- bquote(
            P(F <= .(up)) == .(1 - alpha)
          )
        } else {
          subtitle <- bquote(
            P(F >= .(up)) == .(alpha)
          )
        }
      }
      if (dist == "chisq") {
        up_density <- dchisq(
          x[x >= up],
          ...
        )
        if (ci) {
          ci_right_tail <- dchisq(
            x[x <= up],
            ...
          )
          subtitle <- bquote(
            P(chi^2 <= .(up)) == .(1 - alpha)
          )
        } else {
          subtitle <- bquote(
            P(chi^2 >= .(up)) == .(alpha)
          )
        }
      }
    } else { # left
      lo <- critical[1]
      x_axis <- c(
        lo,
        0
      )
      x <- seq(
        from = lo - 2,
        to = -1 * lo + 2,
        length = 1000
      )
      if (dist == "z") {
        lo_density <- dnorm(
          x[x <= lo],
          ...
        )
        if (ci) {
          ci_left_tail <- dnorm(
            x[x >= lo]
          )
          subtitle <- bquote(
            P(z >= .(lo)) == .(1 - alpha)
          )
        } else {
          subtitle <- bquote(
            P(z <= .(lo)) == .(alpha)
          )
        }
      }
      if (dist == "t") {
        lo_density <- dt(
          x[x <= lo],
          ...
        )
        if (ci) {
          ci_left_tail <- dt(
            x[x >= lo],
            ...
          )
          subtitle <- bquote(
            P(t >= .(lo)) == .(1 - alpha)
          )
        } else {
          subtitle <- bquote(
            P(t <= .(lo)) == .(alpha)
          )
        }
      }
    }
  }
  # For both one-tailed and two-tailed
  if (dist == "z") {
    y <- dnorm(x)
    mu <- 0
    sigma2 <- 1
    xlab <- "z"
    main <- bquote(
      "Normal Distribution" ~ (list(mu == .(mu), sigma^2 == .(sigma2)))
    )
  }
  if (dist == "t") {
    df <- list(...)[["df"]]
    y <- dt(
      x,
      ...
    )
    xlab <- "t"
    main <- bquote(
      "Student's t Distribution" ~ (nu == .(df))
    )
  }
  if (dist == "F") {
    df1 <- list(...)[["df1"]]
    df2 <- list(...)[["df2"]]
    y <- df(
      x,
      ...
    )
    xlab <- "F"
    main <- bquote(
      "Fisher-Snedecor Distribution" ~ (list(d[1] == .(df1), d[2] == .(df2)))
    )
  }
  if (dist == "chisq") {
    df <- list(...)[["df"]]
    y <- dchisq(
      x,
      ...
    )
    xlab <- expression(chi^2)
    main <- bquote(
      chi^2 ~ "Distribution" ~ (k == .(df))
    )
  }
  plot(
    x = x,
    y = y,
    type = "l",
    main = main,
    xlab = xlab,
    ylab = "Density",
    xaxt = "none"
  )
  axis(
    1,
    x_axis
  )
  mtext(
    text = subtitle
  )
  if (two.tailed) {
    if (ci) {
      polygon(
        c(
          lo,
          x[x >= lo & x <= up],
          up
        ),
        c(
          0,
          ci_center,
          0
        ),
        col = "blue",
        border = NA
      )
    } else {
      polygon(
        c(
          min(x),
          x[x <= lo],
          lo
        ),
        c(
          0,
          lo_density,
          0
        ),
        col = "red",
        border = NA
      )
      polygon(
        c(
          up,
          x[x >= up],
          max(x)
        ),
        c(
          0,
          up_density,
          0
        ),
        col = "red",
        border = NA
      )
    }
  } else {
    if (right.tail) {
      if (ci) {
        polygon(
          c(
            min(x),
            x[x <= up],
            up
          ),
          c(
            0,
            ci_right_tail,
            0
          ),
          col = "blue",
          border = NA
        )
      } else {
        polygon(
          c(
            up,
            x[x >= up],
            max(x)
          ),
          c(
            0,
            up_density,
            0
          ),
          col = "red",
          border = NA
        )
      }
    } else {
      if (ci) {
        polygon(
          c(
            lo,
            x[x >= lo],
            max(x)
          ),
          c(
            0,
            ci_left_tail,
            0
          ),
          col = "blue",
          border = NA
        )
      } else {
        polygon(
          c(
            min(x),
            x[x <= lo],
            lo
          ),
          c(
            0,
            lo_density,
            0
          ),
          col = "red",
          border = NA
        )
      }
    }
  }
  if (!is.null(statistic)) {
    abline(
      v = statistic,
      col = "red",
      lwd = 2,
      lty = 3
    )
  }
}
