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
#' @export
zero_hit <- function(thetahat_lo,
                     thetahat_up) {
  theta_hit(
    thetahat_lo,
    theta = 0,
    thetahat_up
  )
}

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
#' @export
shape <- function(thetahat_lo,
                  thetahat,
                  thetahat_up) {
  (thetahat_up - thetahat) / (thetahat - thetahat_lo)
}

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
#'       \mathrm{Var}
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
#'   (\eqn{\hat{\theta}}).
#' @param var_thetahat Numeric.
#'   Variance of `thetahat`
#'   (\eqn{\mathrm{Var} \left( \hat{\theta} \right)}).
#' @param theta_null Numeric.
#'   Hypothesized value of `theta`
#'   (\eqn{\theta_{0}}).
#'   Set to zero by default.
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Wald test statistic.}
#'     \item{p}{p-value.}
#'   }
#' @importFrom stats pchisq
#' @references
#' [Wikipedia: Wald test](https://en.wikipedia.org/wiki/Wald_test)
#' @export
wald_test <- function(thetahat,
                      var_thetahat,
                      theta_null) {
  statistic <- (thetahat - theta_null)^2 / var_thetahat
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
#'       \mathrm{se}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     }
#'   }
#' The associated `p`-value
#' from the `z` or `t` distribution is calculated
#' using [`pnorm()`] or [`pt()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param se_thetahat Numeric.
#'   Standard error of `thetahat`
#'   (\eqn{\mathrm{se} \left( \hat{\theta} \right)}).
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
#' @export
sqrt_wald_test <- function(thetahat,
                           se_thetahat,
                           theta_null = 0,
                           distribution = "z",
                           df) {
  statistic <- (thetahat - theta_null) / se_thetahat
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
#' @export
ci2prob <- function(ci = c(
                      0.999,
                      0.99,
                      0.95
                    )) {
  alpha2prob(alpha = 1 - ci)
}

#' Confidence Interval - Wald
#'
#' Calculates symmetric Wald confidence intervals.
#'
#' As the sample size approaches infinity (\eqn{n \to \infty}),
#' the distribution of \eqn{\hat{\theta}}
#' approaches the normal distribution
#' with mean equal to \eqn{\theta}
#' and variance equal to \eqn{\mathrm{Var} \left( \hat{\theta} \right)}
#'   \deqn{
#'     \hat{\theta}
#'     \mathrel{\dot\sim}
#'     \mathcal{N}
#'     \left(
#'       \theta,
#'       \mathrm{Var}
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
#'       \mathrm{se}
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
#'     \mathrm{se}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     =
#'     \sqrt{
#'       \mathrm{Var}
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
#'     \mathrm{se}
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
#'       \mathrm{se}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'     }
#'     \mathrel{\dot\sim}
#'     t \left( \nu \right)
#'   }
#' where \eqn{\nu} is the degrees of freedom
#' \eqn{n - 1}.
#' As such,
#' the symmetric Wald confidence interval is given by
#'   \deqn{
#'     \hat{\theta}
#'     \pm
#'     t_{\frac{\alpha}{2}}
#'     \times
#'     \mathrm{se}
#'       \left(
#'         \hat{\theta}
#'       \right)
#'   }
#' from a \eqn{t} distribution with degrees of freedom \eqn{\nu = n - 1}.
#' Note that in large sample sizes,
#' \eqn{t} converges to \eqn{z}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams sqrt_wald_test
#' @inheritParams alpha2prob
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic.}
#'     \item{p}{p-value.}
#'     \item{se}{Standard error of thetahat (\eqn{\mathrm{se} \left( \hat{\theta} \right)}).}
#'     \item{ci_*}{Estimated confidence limits corresponding to alpha.}
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
#' se_thetahat <- sd(x) / sqrt(n)
#' se_thetahat
#'
#' #############################################################
#' # Generate Wald confidence intervals
#' # for alpha = c(0.001, 0.01, 0.05).
#' #############################################################
#' wald(
#'   thetahat = thetahat,
#'   se_thetahat = se_thetahat
#' )
#' @references
#' [Wikipedia: Confidence interval](https://en.wikipedia.org/wiki/Confidence_interval)
#' @export
wald <- function(thetahat,
                 se_thetahat,
                 theta_null = 0,
                 alpha = c(
                   0.001,
                   0.01,
                   0.05
                 ),
                 distribution = "z",
                 df) {
  sqrt_W <- sqrt_wald_test(
    thetahat = thetahat,
    se_thetahat = se_thetahat,
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
    ci[i] <- thetahat + (critical[i] * se_thetahat)
  }
  names(ci) <- paste0(
    "ci_",
    prob * 100
  )
  c(
    sqrt_W,
    se = se_thetahat,
    ci
  )
}

#' Confidence Interval - Percentile
#'
#' Calculates percentile confidence intervals.
#'
#' The standard error of the bootstrap sampling distribution
#' \eqn{\hat{\theta}^{*}}
#' is given by
#' \deqn{
#'   \mathrm{se}
#'   \left(
#'     \hat{\theta}^{*}
#'   \right)
#'   =
#'   \sqrt{
#'     \frac{1}{B - 1}
#'     \sum_{b = 1}^{B}
#'     \left[
#'       \hat{\theta}^{*} \left( b \right)
#'       -
#'       \hat{\theta}^{*} \left( \cdot \right)
#'     \right]^2
#'   }
#' }
#' where
#' the mean of the bootstrap sampling distribution
#' \eqn{\hat{\theta}^{*}}
#' is given by
#' \deqn{
#'   \hat{\theta}^{*}
#'   \left(
#'     \cdot
#'   \right)
#'   =
#'   \frac{1}{n}
#'   \sum_{b = 1}^{B}
#'   \hat{\theta}^{*} \left( b \right) .
#' }
#'
#' The percentile confidence interval is given by
#'   \deqn{
#'     \left[
#'       \hat{\theta}_{\mathrm{lo}},
#'       \hat{\theta}_{\mathrm{up}}
#'     \right]
#'     =
#'     \left[
#'       \hat{\theta}^{*}_{\frac{\alpha}{2}},
#'       \hat{\theta}^{*}_{1 - \frac{\alpha}{2}}
#'     \right] .
#'   }
#' The corresponding percentile value
#' from the bootstrap sampling distribution `thetahat_star`
#' (\eqn{\hat{\theta}^{*}}) is calculated using
#' [`quantile()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_star Numeric vector.
#'   The bootstrap sampling distribution (\eqn{\hat{\theta}^{*}}).
#' @param wald Logical.
#'   If `TRUE`,
#'   calculates the square root of the Wald test statistic and p-value.
#'   If `FALSE`,
#'   returns `statistic = NA`
#'   and
#'   `p = NA`.
#'   If `FALSE`,
#'   the arguments
#'   `thetahat`,
#'   `theta_null`,
#'   `distribution`,
#'   and
#'   `df`
#'   are ignored.
#' @inheritParams wald
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic. `NA` if `wald = FALSE`.}
#'     \item{p}{p-value. `NA` if `wald = FALSE`.}
#'     \item{se}{Standard error of thetahat_star (\eqn{\mathrm{se} \left( \hat{\theta}^{*} \right)}), i.e., standard deviation of thetahat_star (\eqn{\mathrm{sd} \left( \hat{\theta}^{*} \right)}).}
#'     \item{ci_*}{Estimated percentile confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star (\eqn{\hat{\theta}^{*}}).}
#'   }
#' @importFrom stats quantile
#' @importFrom stats sd
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
#' se_thetahat <- sd(x) / sqrt(n)
#' se_thetahat
#'
#' #############################################################
#' # Generate B = 2000 nonparametric bootstrap samples.
#' #############################################################
#' x_star <- nb(
#'   data = x,
#'   B = 2000L,
#'   par = FALSE,
#'   ncores = NULL
#' )
#'
#' #############################################################
#' # Estimate xbar for each of the nonparametric bootstrap samples.
#' # The B = 2000 xbars form
#' # the empirical sampling distribution thetahat_star.
#' #############################################################
#' thetahat_star <- lapply(
#'   X = x_star,
#'   FUN = mean
#' )
#' thetahat_star <- as.vector(
#'   do.call(
#'     what = "rbind",
#'     args = thetahat_star
#'   )
#' )
#' hist(thetahat_star)
#' summary(thetahat_star)
#'
#' #############################################################
#' # Generate percentile confidence intervals
#' # for alpha = c(0.001, 0.01, 0.05).
#' #############################################################
#' # With square root of Wald test statistic and p-value
#' # using sd(thetahat_star) as se_thetahat
#' pc(
#'   thetahat_star = thetahat_star,
#'   wald = TRUE,
#'   thetahat = thetahat
#' )
#'
#' # Without square root of Wald test statistic and p-value
#' pc(
#'   thetahat_star = thetahat_star
#' )
#' @references
#' Bradley Efron and Robert J. Tibshirani, 1993, "An Introduction to the Bootstrap"
#' @export
pc <- function(thetahat_star,
               alpha = c(
                 0.001,
                 0.01,
                 0.05
               ),
               wald = FALSE,
               thetahat,
               theta_null = 0,
               distribution = "z",
               df) {
  se_thetahat_star <- sd(thetahat_star)
  if (wald) {
    sqrt_W <- sqrt_wald_test(
      thetahat = thetahat,
      se_thetahat = se_thetahat_star,
      theta_null = theta_null,
      distribution = distribution,
      df = df
    )
  } else {
    sqrt_W <- c(
      statistic = NA,
      p = NA
    )
  }
  probs <- alpha2prob(alpha = alpha)
  ci <- quantile(
    x = thetahat_star,
    probs = probs,
    names = FALSE
  )
  names(ci) <- paste0(
    "ci_",
    probs * 100
  )
  c(
    sqrt_W,
    se = se_thetahat_star,
    ci
  )
}

#' Confidence Interval - Bias Corrected
#'
#' Calculates bias corrected confidence intervals.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param wald Logical.
#'   If `TRUE`,
#'   calculates the square root of the Wald test statistic and p-value.
#'   If `FALSE`,
#'   returns `statistic = NA`
#'   and
#'   `p = NA`.
#'   If `FALSE`,
#'   the arguments
#'   `theta_null`,
#'   `distribution`,
#'   and
#'   `df`
#'   are ignored.
#' @inheritParams pc
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic. `NA` if `wald = FALSE`.}
#'     \item{p}{p-value. `NA` if `wald = FALSE`.}
#'     \item{se}{Standard error of thetahat_star (\eqn{\mathrm{se} \left( \hat{\theta}^{*} \right)}), i.e., standard deviation of thetahat_star (\eqn{\mathrm{sd} \left( \hat{\theta}^{*} \right)}).}
#'     \item{ci_*}{Estimated bias corrected confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star (\eqn{\hat{\theta}^{*}}).}
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
#' se_thetahat <- sd(x) / sqrt(n)
#' se_thetahat
#'
#' #############################################################
#' # Generate B = 2000 nonparametric bootstrap samples.
#' #############################################################
#' x_star <- nb(
#'   data = x,
#'   B = 2000L,
#'   par = FALSE,
#'   ncores = NULL
#' )
#'
#' #############################################################
#' # Estimate xbar for each of the nonparametric bootstrap samples.
#' # The B = 2000 xbars form
#' # the empirical sampling distribution thetahat_star.
#' #############################################################
#' thetahat_star <- lapply(
#'   X = x_star,
#'   FUN = mean
#' )
#' thetahat_star <- as.vector(
#'   do.call(
#'     what = "rbind",
#'     args = thetahat_star
#'   )
#' )
#' hist(thetahat_star)
#' summary(thetahat_star)
#'
#' #############################################################
#' # Generate percentile confidence intervals
#' # for alpha = c(0.001, 0.01, 0.05).
#' #############################################################
#' # With square root of Wald test statistic and p-value
#' # using sd(thetahat_star) as se_thetahat
#' bc(
#'   thetahat_star = thetahat_star,
#'   wald = TRUE,
#'   thetahat = thetahat
#' )
#'
#' # Without square root of Wald test statistic and p-value
#' bc(
#'   thetahat_star = thetahat_star,
#'   thetahat = thetahat
#' )
#' @inherit pc references
#' @export
bc <- function(thetahat_star,
               thetahat,
               alpha = c(
                 0.001,
                 0.01,
                 0.05
               ),
               wald = FALSE,
               theta_null = 0,
               distribution = "z",
               df) {
  z0 <- qnorm(
    sum(thetahat_star < thetahat) / length(thetahat_star)
  )
  se_thetahat_star <- sd(thetahat_star)
  if (wald) {
    sqrt_W <- sqrt_wald_test(
      thetahat = thetahat,
      se_thetahat = se_thetahat_star,
      theta_null = theta_null,
      distribution = distribution,
      df = df
    )
  } else {
    sqrt_W <- c(
      statistic = NA,
      p = NA
    )
  }
  probs <- alpha2prob(alpha = alpha)
  bc_probs <- pnorm(
    q = 2 * z0 + qnorm(
      p = probs
    )
  )
  ci <- quantile(
    x = thetahat_star,
    probs = bc_probs,
    names = FALSE
  )
  names(ci) <- paste0(
    "ci_",
    probs * 100
  )
  c(
    sqrt_W,
    se = se_thetahat_star,
    ci
  )
}

#' Confidence Interval - Bias Corrected and Accelerated
#'
#' Calculates bias corrected and accelerated confidence intervals.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param data Vector, matrix, or data frame.
#'   Sample data.
#' @param fitFUN Function.
#'   Fit function to use on `data`.
#'   The first argument should correspond to `data`.
#'   Other arguments can be passed to `fitFUN`
#'   using `...`.
#'   `fitFUN` should return a single value.
#' @param ... Arguments to pass to `fitFUN`.
#' @inheritParams bc
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic. `NA` if `wald = FALSE`.}
#'     \item{p}{p-value. `NA` if `wald = FALSE`.}
#'     \item{se}{Standard error of thetahat_star (\eqn{\mathrm{se} \left( \hat{\theta}^{*} \right)}), i.e., standard deviation of thetahat_star (\eqn{\mathrm{sd} \left( \hat{\theta}^{*} \right)}).}
#'     \item{ci_*}{Estimated bias corrected and accelerated confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star (\eqn{\hat{\theta}^{*}}).}
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
#' se_thetahat <- sd(x) / sqrt(n)
#' se_thetahat
#'
#' #############################################################
#' # Generate B = 2000 nonparametric bootstrap samples.
#' #############################################################
#' x_star <- nb(
#'   data = x,
#'   B = 2000L,
#'   par = FALSE,
#'   ncores = NULL
#' )
#'
#' #############################################################
#' # Estimate xbar for each of the nonparametric bootstrap samples.
#' # The B = 2000 xbars form
#' # the empirical sampling distribution thetahat_star.
#' #############################################################
#' thetahat_star <- lapply(
#'   X = x_star,
#'   FUN = mean
#' )
#' thetahat_star <- as.vector(
#'   do.call(
#'     what = "rbind",
#'     args = thetahat_star
#'   )
#' )
#' hist(thetahat_star)
#' summary(thetahat_star)
#'
#' #############################################################
#' # Generate percentile confidence intervals
#' # for alpha = c(0.001, 0.01, 0.05).
#' #############################################################
#' # With square root of Wald test statistic and p-value
#' # using sd(thetahat_star) as se_thetahat
#' bca(
#'   thetahat_star = thetahat_star,
#'   wald = TRUE,
#'   thetahat = thetahat,
#'   data = x,
#'   fitFUN = mean
#' )
#'
#' # Without square root of Wald test statistic and p-value
#' bca(
#'   thetahat_star = thetahat_star,
#'   thetahat = thetahat,
#'   data = x,
#'   fitFUN = mean
#' )
#' @inherit pc references
#' @export
bca <- function(thetahat_star,
                thetahat,
                data,
                fitFUN,
                alpha = c(
                  0.001,
                  0.01,
                  0.05
                ),
                wald = FALSE,
                theta_null = 0,
                distribution = "z",
                df,
                ...) {
  n <- length(thetahat_star)
  z0 <- qnorm(
    sum(thetahat_star < thetahat) / length(thetahat_star)
  )
  probs <- alpha2prob(alpha = alpha)
  z1 <- qnorm(
    p = probs
  )
  se_thetahat_star <- sd(thetahat_star)
  if (wald) {
    sqrt_W <- sqrt_wald_test(
      thetahat = thetahat,
      se_thetahat = se_thetahat_star,
      theta_null = theta_null,
      distribution = distribution,
      df = df
    )
  } else {
    sqrt_W <- c(
      statistic = NA,
      p = NA
    )
  }
  u <- rep(x = 0, times = n)
  for (i in 1:n) {
    if (is.vector(data)) {
      data_minus_i <- data[-i]
    }
    if (is.matrix(data) | is.data.frame(data)) {
      data_minus_i <- data[-i, ]
    }
    u[i] <- fitFUN(
      data_minus_i,
      ...
    )
  }
  uu <- mean(u) - u
  acc <- sum(uu * uu * uu) / (6 * (sum(uu * uu))^1.5)
  bca_probs <- pnorm(
    z0 + (z0 + z1) / (1 - acc * (z0 + z1))
  )
  ci <- quantile(
    x = thetahat_star,
    probs = bca_probs,
    names = FALSE
  )
  names(ci) <- paste0(
    "ci_",
    probs * 100
  )
  c(
    sqrt_W,
    se = se_thetahat_star,
    ci
  )
}
