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
#'   \frac{1}{B}
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
#'     \item{ci_}{Estimated percentile confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star (\eqn{\hat{\theta}^{*}}).}
#'   }
#' If `eval = TRUE`,
#' also returns
#'   \describe{
#'     \item{zero_hit_}{Logical. Tests if confidence interval contains zero.}
#'     \item{theta_hit_}{Logical. Tests if confidence interval contains theta.}
#'     \item{length_}{Length of confidence interval.}
#'     \item{shape_}{Shape of confidence interval.}
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
#'   thetahat_star = thetahat_star,
#'   thetahat = thetahat
#' )
#' @references
#' Efron, B., & Tibshirani, R. J. (1993).
#' An introduction to the bootstrap. New York, N.Y: Chapman & Hall.
#' @family bootstrap confidence interval functions
#' @keywords confidence interval
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
               df,
               eval = FALSE,
               theta = 0) {
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
  out <- c(
    sqrt_W,
    se = se_thetahat_star,
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
#'     \item{ci_}{Estimated bias corrected confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star (\eqn{\hat{\theta}^{*}}).}
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
#' @family bootstrap confidence interval functions
#' @keywords confidence interval
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
               df,
               eval = FALSE,
               theta = 0) {
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
  out <- c(
    sqrt_W,
    se = se_thetahat_star,
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
#'     \item{ci_}{Estimated bias corrected and accelerated confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star (\eqn{\hat{\theta}^{*}}).}
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
#' @family bootstrap confidence interval functions
#' @keywords confidence interval
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
                eval = FALSE,
                theta = 0,
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
  out <- c(
    sqrt_W,
    se = se_thetahat_star,
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
