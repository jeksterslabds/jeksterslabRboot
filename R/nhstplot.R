#' Null Hypothesis Significance Testing Plot
#'
#' Generates a plot to illustrate
#' null hypothesis significance testing.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family alpha functions
#' @keywords alpha
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
#' nhstplot(
#'   alpha = 0.05,
#'   dist = "z",
#'   statistic = 3
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   dist = "t",
#'   df = 5
#' )
#' # one-tailed
#' nhstplot(
#'   alpha = 0.05,
#'   two.tailed = FALSE,
#'   right.tail = TRUE
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = TRUE,
#'   df = 5
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   dist = "F",
#'   df1 = 10,
#'   df2 = 3
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   dist = "chisq",
#'   df = 3
#' )
#' # ci = TRUE
#' # two-tailed
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "z",
#'   statistic = 3
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "t",
#'   df = 5
#' )
#' # one-tailed
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   two.tailed = FALSE,
#'   right.tail = TRUE
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   two.tailed = FALSE,
#'   right.tail = FALSE
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = TRUE,
#'   df = 5
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "t",
#'   two.tailed = FALSE,
#'   right.tail = FALSE,
#'   df = 5
#' )
#' nhstplot(
#'   alpha = 0.05,
#'   ci = TRUE,
#'   dist = "F",
#'   df1 = 10,
#'   df2 = 3
#' )
#' nhstplot(
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
nhstplot <- function(alpha = 0.05,
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
            P(z < .(up)) == .(alpha)
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
            P(z < .(up)) == .(1 - alpha)
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
            P(F < .(up)) == .(1 - alpha)
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
            P(chi^2 < .(up)) == .(1 - alpha)
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
            P(z > .(lo)) == .(1 - alpha)
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
            P(t > .(lo)) == .(1 - alpha)
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
