#' Wald Test Statistic for a Single Parameter
#'
#' Calculates the Wald test statistic for a single parameter.
#'
#' The Wald test statistic for a single parameter is calculated
#' as follows:
#'   \deqn{
#'     W
#'     =
#'     \frac{
#'       \left(
#'         \hat{\theta} - \theta_0
#'       \right)^2
#'     }
#'     {
#'       \widehat{\mathrm{Var}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     }
#'   }
#' The associated `p`-value
#' from the `chi-square` distribution is calculated
#' using [`pchisq()`] with `df = 1`.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat Numeric.
#'   Parameter estimate
#'   \eqn{\left( \hat{\theta} \right)}.
#' @param varhat_thetahat Numeric.
#'   Estimated variance of `thetahat`
#'   \eqn{\left( \widehat{\mathrm{Var}} \left( \hat{\theta} \right) \right)}.
#' @param theta_null Numeric.
#'   Hypothesized value of `theta`
#'   \eqn{\left( \theta_{0} \right)}.
#'   Set to zero by default.
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Wald test statistic.}
#'     \item{p}{p-value.}
#'   }
#' @importFrom stats pchisq
#' @references
#' [Wikipedia: Wald test](https://en.wikipedia.org/wiki/Wald_test)
#' @family Wald test functions
#' @keywords Wald test
#' @examples
#' thetahat <- 0.860
#' sehat_thetahat <- 0.409
#' varhat_thetahat <- sehat_thetahat^2
#' wald_test(
#'   thetahat = thetahat,
#'   varhat_thetahat = varhat_thetahat
#' )
#' @export
wald_test <- function(thetahat,
                      varhat_thetahat,
                      theta_null = 0) {
  statistic <- (thetahat - theta_null)^2 / varhat_thetahat
  p <- 2 * pchisq(
    q = statistic,
    df = 1,
    lower.tail = FALSE
  )
  c(
    statistic = statistic,
    p = p
  )
}

#' Square Root of the Wald Test Statistic for a Single Parameter
#'
#' Calculates the square root of the Wald test statistic for a single parameter.
#'
#' The square root of the Wald test statistic for a single parameter is calculated
#' as follows:
#'   \deqn{
#'     \sqrt{W}
#'     =
#'     \frac{
#'       \hat{\theta} - \theta_0
#'     }
#'     {
#'       \widehat{\mathrm{se}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     } .
#'   }
#' The associated `p`-value
#' from the `z` or `t` distribution is calculated
#' using [`pnorm()`] or [`pt()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param sehat_thetahat Numeric.
#'   Estimated standard error of `thetahat`
#'   \eqn{\left( \widehat{\mathrm{se}} \left( \hat{\theta} \right) \right)}.
#' @param distribution Character string.
#'   `distribution = "z"` for the standard normal distribution.
#'   `distribution = "t"` for the t distribution.
#' @param df Numeric.
#'   Degrees of freedom (df) if `distribution = "t"`.
#'   Ignored if `distribution = "z"`.
#' @inheritParams wald_test
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic.}
#'     \item{p}{p-value.}
#'   }
#' @importFrom stats pnorm
#' @importFrom stats pt
#' @inherit wald_test references
#' @family Wald test functions
#' @keywords Wald test
#' @examples
#' thetahat <- 0.860
#' sehat_thetahat <- 0.409
#' sqrt_wald_test(
#'   thetahat = thetahat,
#'   sehat_thetahat = sehat_thetahat
#' )
#' sqrt_wald_test(
#'   thetahat = thetahat,
#'   sehat_thetahat = sehat_thetahat,
#'   distribution = "t",
#'   df = 1000
#' )
#' @export
sqrt_wald_test <- function(thetahat,
                           sehat_thetahat,
                           theta_null = 0,
                           distribution = "z",
                           df) {
  statistic <- (thetahat - theta_null) / sehat_thetahat
  if (distribution == "z") {
    p <- 2 * pnorm(
      q = statistic,
      lower.tail = FALSE
    )
  }
  if (distribution == "t") {
    p <- 2 * pt(
      q = statistic,
      df = df,
      lower.tail = FALSE
    )
  }
  c(
    statistic = statistic,
    p = p
  )
}

#' Confidence Interval - Wald
#'
#' Calculates symmetric Wald confidence intervals.
#'
#' As the sample size approaches infinity \eqn{\left( n \to \infty \right)},
#' the distribution of \eqn{\hat{\theta}}
#' approaches the normal distribution
#' with mean equal to \eqn{\theta}
#' and variance equal to \eqn{\widehat{\mathrm{Var}} \left( \hat{\theta} \right)}
#'   \deqn{
#'     \hat{\theta}
#'     \mathrel{\dot\sim}
#'     \mathcal{N}
#'     \left(
#'       \theta,
#'       \widehat{\mathrm{Var}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     \right) .
#'   }
#' As such,
#'   \eqn{
#'     \hat{\theta}
#'   }
#' can be expressed in terms of
#' \eqn{z}-scores
#' from a standard normal distribution
#'   \deqn{
#'     \frac{
#'       \hat{\theta} - \theta
#'     }
#'     {
#'       \widehat{\mathrm{se}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     }
#'     \mathrel{\dot\sim}
#'     \mathcal{N}
#'     \left(
#'       0,
#'       1
#'     \right)
#'   }
#' where
#'   \eqn{
#'     \widehat{\mathrm{se}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     =
#'     \sqrt{
#'       \widehat{\mathrm{Var}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     }
#'   }.
#'
#' To form a confidence interval around \eqn{\hat{\theta}},
#' the \eqn{z}-score associated with a particular
#' alpha level
#' can be plugged-in the equation below.
#'   \deqn{
#'     \hat{\theta}
#'     \pm
#'     z_{\frac{\alpha}{2}}
#'     \times
#'     \widehat{\mathrm{se}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'   }
#'
#' Note that this is valid only when \eqn{n \to \infty}.
#' In finite samples,
#' this is only an approximation.
#' Gosset derived a better approximation
#' in the context of \eqn{\hat{\theta} = \bar{x}}
#'   \deqn{
#'     \frac{
#'       \hat{\theta} - \theta
#'     }
#'     {
#'       \widehat{\mathrm{se}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     }
#'     \mathrel{\dot\sim}
#'     t \left( \nu \right)
#'   }
#' where
#' \eqn{t}
#' is the Student's
#' \eqn{t} distribution
#' and
#' \eqn{\nu} is the degrees of freedom
#' \eqn{n - 1}.
#' As such,
#' the symmetric Wald confidence interval is given by
#'   \deqn{
#'     \hat{\theta}
#'     \pm
#'     t_{\frac{\alpha}{2}}
#'     \times
#'     \widehat{\mathrm{se}}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'   }
#' from a \eqn{t} distribution with degrees of freedom \eqn{\nu = n - 1}.
#' Note that in large sample sizes,
#' \eqn{t} converges to \eqn{z}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param eval Logical.
#'   Evaluate confidence intervals using
#'   [`zero_hit()`],
#'   [`theta_hit()`],
#'   [`len()`],
#'   and
#'   [`shape()`].
#' @inheritParams sqrt_wald_test
#' @inheritParams alpha2prob
#' @inheritParams ci_eval
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic.}
#'     \item{p}{p-value.}
#'     \item{se}{Estimated d=standard error of thetahat \eqn{\left( \widehat{\mathrm{se}} \left( \hat{\theta} \right) \right)}.}
#'     \item{ci_}{Estimated confidence limits corresponding to alpha.}
#'   }
#' If `eval = TRUE`,
#' also returns
#'   \describe{
#'     \item{zero_hit_}{Logical. Tests if confidence interval contains zero.}
#'     \item{theta_hit_}{Logical. Tests if confidence interval contains theta.}
#'     \item{length_}{Length of confidence interval.}
#'     \item{shape_}{Shape of confidence interval.}
#'   }
#' @examples
#' #############################################################
#' # Generate sample data from a normal distribution
#' # with mu = 0 and sigma^2 = 1.
#' #############################################################
#' set.seed(42)
#' n <- 10000
#' x <- rnorm(n = n)
#'
#' #############################################################
#' # Estimate the population mean mu using the sample mean xbar.
#' #############################################################
#' # Parameter estimate xbar
#' thetahat <- mean(x)
#' thetahat
#' # Closed form solution for the standard error of the mean
#' sehat_thetahat <- sd(x) / sqrt(n)
#' sehat_thetahat
#'
#' #############################################################
#' # Generate Wald confidence intervals
#' # for alpha = c(0.001, 0.01, 0.05).
#' #############################################################
#' wald(
#'   thetahat = thetahat,
#'   sehat_thetahat = sehat_thetahat
#' )
#' @references
#' [Wikipedia: Confidence interval](https://en.wikipedia.org/wiki/Confidence_interval)
#' @family Wald confidence interval functions
#' @keywords confidence interval, Wald test
#' @export
wald <- function(thetahat,
                 sehat_thetahat,
                 theta_null = 0,
                 alpha = c(
                   0.001,
                   0.01,
                   0.05
                 ),
                 distribution = "z",
                 df,
                 eval = FALSE,
                 theta = 0) {
  sqrt_W <- sqrt_wald_test(
    thetahat = thetahat,
    sehat_thetahat = sehat_thetahat,
    theta_null = theta_null,
    distribution = distribution,
    df = df
  )
  prob <- alpha2prob(alpha = alpha)
  if (distribution == "z") {
    critical <- qnorm(
      p = prob
    )
  }
  if (distribution == "t") {
    critical <- qt(
      p = prob,
      df = df
    )
  }
  ci <- rep(x = NA, times = length(critical))
  for (i in seq_along(critical)) {
    ci[i] <- thetahat + (critical[i] * sehat_thetahat)
  }
  names(ci) <- paste0(
    "ci_",
    prob * 100
  )
  out <- c(
    sqrt_W,
    se = sehat_thetahat,
    ci
  )
  if (eval) {
    ci_eval <- ci_eval(
      ci = ci,
      thetahat = thetahat,
      theta = theta,
      label = alpha
    )
    out <- c(
      out,
      ci_eval
    )
  }
  out
}
