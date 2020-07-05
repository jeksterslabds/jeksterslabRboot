#' Confidence Interval - Percentile
#'
#' Calculates percentile confidence intervals.
#'
#' The estimated bootstrap standard error
#' is given by
#' \deqn{
#'   \widehat{\mathrm{se}}_{\mathrm{B}}
#'   \left(
#'     \hat{\theta}
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
#' Note that
#' \eqn{
#'   \widehat{\mathrm{se}}_{\mathrm{B}}
#'   \left(
#'     \hat{\theta}
#'   \right)
#' }
#' is the standard deviation of
#' \eqn{
#'   \boldsymbol{\hat{\theta}^{*}}
#' }
#' and
#' \eqn{
#'   \hat{\theta}^{*}
#'   \left(
#'     \cdot
#'   \right)
#' }
#' is the mean of
#' \eqn{
#'   \boldsymbol{\hat{\theta}^{*}}
#' } .
#'
#' The percentile confidence interval is given by
#'   \deqn{
#'     \left[
#'       \hat{\theta}_{\mathrm{lo}},
#'       \hat{\theta}_{\mathrm{up}}
#'     \right]
#'     =
#'     \left[
#'       \hat{\theta}^{*}_{z_{\left( \frac{\alpha}{2} \right)}},
#'       \hat{\theta}^{*}_{z_{\left( 1 - \frac{\alpha}{2} \right)}}
#'     \right] .
#'   }
#'
#' For more details and examples see the following vignettes:
#'
#' [Notes: Introduction to Nonparametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Introduction to Parametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_star Numeric vector.
#'   The bootstrap sampling distribution
#'   \eqn{\left( \boldsymbol{\hat{\theta}^{*}} \right)},
#'   that is,
#'   the sampling distribution of `thetahat`
#'   estimated for each `b` bootstrap sample.
#'   \eqn{
#'     \hat{\theta}_{\left( 1 \right)},
#'     \hat{\theta}_{\left( 2 \right)},
#'     \hat{\theta}_{\left( 3 \right)},
#'     \hat{\theta}_{\left( b \right)},
#'     \dots,
#'     \hat{\theta}_{\left( B \right)}
#'   } .
#' @param thetahat Numeric.
#'   Parameter estimate
#'   \eqn{\left( \hat{\theta} \right)}
#'   from the original sample data.
#' @param wald Logical.
#'   If `TRUE`,
#'   calculates the square root of the Wald test statistic and p-value.
#'   The estimated bootstrap standard error is used.
#'   The arguments
#'   `theta_null`,
#'   `distribution`,
#'   and
#'   `df`
#'   are used to calculate
#'   the square root of the Wald test statistic and p-value
#'   and are NOT used in constructing the bootstrap confidence interval.
#'   If `FALSE`,
#'   returns `statistic = NA`
#'   and
#'   `p = NA`
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
#'     \item{se}{Estimated bootstrap standard error \eqn{\left( \widehat{\mathrm{se}}_{\mathrm{B}} \left( \hat{\theta} \right) \right)}.}
#'     \item{ci_}{Estimated percentile confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star \eqn{\left( \boldsymbol{\hat{\theta}^{*}} \right)}.}
#'   }
#' If `eval = TRUE`,
#' appends the following to the results vector
#'   \describe{
#'     \item{zero_hit_}{Logical. Tests if confidence interval contains zero.}
#'     \item{theta_hit_}{Logical. Tests if confidence interval contains theta.}
#'     \item{length_}{Length of confidence interval.}
#'     \item{shape_}{Shape of confidence interval.}
#'   }
#' @importFrom stats quantile
#' @importFrom stats sd
#' @references
#' Efron, B., & Tibshirani, R. J. (1993).
#' *An introduction to the bootstrap*.
#' New York, N.Y: Chapman & Hall.
#' @family bootstrap confidence interval functions
#' @keywords confidence interval
#' @export
pc <- function(thetahat_star,
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
  # standard error
  sehat_B_thetahat <- sd(thetahat_star)
  # confidence interval
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
  # optional - wald
  if (wald) {
    sqrt_W <- sqrt_wald_test(
      thetahat = thetahat,
      sehat_thetahat = sehat_B_thetahat,
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
  # initial output
  out <- c(
    sqrt_W,
    se = sehat_B_thetahat,
    ci
  )
  # optional - ci evaluation
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

#' Confidence Interval - Bias-Corrected
#'
#' Calculates bias-corrected confidence intervals.
#'
#' The estimated bootstrap standard error
#' is given by
#' \deqn{
#'   \widehat{\mathrm{se}}_{\mathrm{B}}
#'   \left(
#'     \hat{\theta}
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
#' Bias-correction
#' \eqn{\hat{z}_{0}}
#' has to be computed first
#' and used to adjust the percentile ranks
#' used in generating confidence intervals.
#'
#' The bias-correction \eqn{\hat{z}_{0}}
#' is given by
#'
#' \deqn{
#'   \hat{z}_{0}
#'   =
#'   \Phi^{-1}
#'   \left(
#'     \frac{
#'       \#
#'       \left\{
#'         {\hat{\theta}}^{*}
#'         \left(
#'         b
#'       \right)
#'       <
#'      \hat{\theta}
#'      \right\}
#'     }
#'     {
#'       B
#'     }
#'   \right)
#' }
#'
#' where
#' the quantity inside the parenthesis
#' is the proportion of bootstrap replications
#' less than the original parameter estimate \eqn{\hat{\theta}}
#' and
#' \eqn{
#'   \Phi^{-1}
#'   \left(
#'     \cdot
#'   \right)
#' }
#' is the inverse function
#' of a standard normal cumulative distribution function
#' ([`qnorm()`]).
#'
#' \eqn{\hat{z}_{0}}
#' can then be used to obtain adjusted \eqn{z}-scores
#' for the lower limit and the upper limit of the confidence interval
#' as follows
#'
#' \deqn{
#'   z_{
#'     \mathrm{BC}_{
#'       \mathrm{lo}
#'     }
#'   }
#'   =
#'   2
#'   \hat{z_0}
#'   +
#'   z_{
#'     \left(
#'       \frac{\alpha}{2}
#'     \right)
#'   } ,
#' }
#'
#' \deqn{
#'   z_{
#'     \mathrm{BC}_{
#'       \mathrm{up}
#'     }
#'   }
#'   =
#'   2
#'   \hat{z_0}
#'   +
#'   z_{
#'     \left(
#'       1 - \frac{\alpha}{2}
#'     \right)
#'   } .
#' }
#'
#' The adjusted \eqn{z}-scores
#' are used to determine the adjusted
#' percentile ranks
#' to form the confidence interval.
#'
#' The bias-corrected confidence interval is given by
#' \deqn{
#'   \left[
#'     \hat{\theta}_{\mathrm{lo}},
#'     \hat{\theta}_{\mathrm{up}}
#'   \right]
#'   =
#'   \left[
#'     \hat{\theta}^{*}_{
#'       \left(
#'         z_{
#'           \mathrm{BC}_{
#'             \mathrm{lo}
#'           }
#'         }
#'       \right)
#'     },
#'     \hat{\theta}^{*}_{
#'       \left(
#'         z_{
#'           \mathrm{BC}_{
#'             \mathrm{up}
#'           }
#'         }
#'       \right)
#'     }
#'   \right] .
#' }
#'
#' For more details and examples see the following vignettes:
#'
#' [Notes: Introduction to Nonparametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Introduction to Parametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams pc
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic. `NA` if `wald = FALSE`.}
#'     \item{p}{p-value. `NA` if `wald = FALSE`.}
#'     \item{se}{Estimated bootstrap standard error \eqn{\left( \widehat{\mathrm{se}}_{\mathrm{B}} \left( \hat{\theta} \right) \right)}.}
#'     \item{ci_}{Estimated bias-corrected confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star \eqn{\left( \boldsymbol{\hat{\theta}^{*}} \right)}.}
#'   }
#' If `eval = TRUE`,
#' appends the following to the results vector
#'   \describe{
#'     \item{zero_hit_}{Logical. Tests if confidence interval contains zero.}
#'     \item{theta_hit_}{Logical. Tests if confidence interval contains theta.}
#'     \item{length_}{Length of confidence interval.}
#'     \item{shape_}{Shape of confidence interval.}
#'   }
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
  # standard error
  sehat_B_thetahat <- sd(thetahat_star)
  # bias-correction
  z0hat <- qnorm(
    sum(thetahat_star < thetahat) / length(thetahat_star)
  )
  # confidence interval
  probs <- alpha2prob(alpha = alpha)
  bc_probs <- pnorm(
    q = 2 * z0hat + qnorm(
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
  # optional - wald
  if (wald) {
    sqrt_W <- sqrt_wald_test(
      thetahat = thetahat,
      sehat_thetahat = sehat_B_thetahat,
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
  # initial output
  out <- c(
    sqrt_W,
    se = sehat_B_thetahat,
    ci
  )
  # optional - ci evaluation
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

#' Confidence Interval - Bias-Corrected and Accelerated
#' (from \eqn{\hat{\theta}^{*}} and \eqn{\hat{\theta}^{*}_{\mathrm{jack}}})
#'
#' Calculates bias-corrected and accelerated confidence intervals.
#'
#' The estimated bootstrap standard error
#' is given by
#' \deqn{
#'   \widehat{\mathrm{se}}_{\mathrm{B}}
#'   \left(
#'     \hat{\theta}
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
#' In addition to the bias-correction
#' \eqn{\hat{z}_{0}}
#' discussed in [`bc()`],
#' the acceleration
#' \eqn{\hat{a}},
#' which refers to the rate of change
#' of the standard error of
#' \eqn{\hat{\theta}}
#' with respect to the true parameter value
#' \eqn{\theta},
#' is factored in.
#'
#' The acceleration \eqn{\hat{a}}
#' is given by
#'
#' \deqn{
#'   \hat{a}
#'   =
#'   \frac{
#'     \sum_{i = 1}^{n}
#'     \left[
#'       \hat{\theta}_{
#'         \left(
#'           \cdot
#'         \right)
#'       }
#'       -
#'       \hat{\theta}_{
#'         \left(
#'           i
#'         \right)
#'       }
#'     \right]^{3}
#'   }
#'   {
#'     6
#'     \left\{
#'       \sum_{i = 1}^{n}
#'       \left[
#'         \hat{\theta}_{
#'           \left(
#'             \cdot
#'           \right)
#'         }
#'         -
#'         \hat{\theta}_{
#'           \left(
#'             i
#'           \right)}
#'       \right]^{2}
#'     \right\}^{3/2}
#'   }
#' }
#' where
#' \deqn{
#'   \hat{\theta}_{
#'     \left(
#'       i
#'     \right)
#'   }
#'   =
#'   \left\{
#'     \hat{\theta}_{\left( 1 \right)},
#'     \hat{\theta}_{\left( 2 \right)},
#'     \hat{\theta}_{\left( 3 \right)},
#'     \dots,
#'     \hat{\theta}_{\left( n \right)}
#'   \right\}
#' }
#' is the jackknife sampling distribution
#' and
#' \deqn{
#'   \hat{\theta}_{
#'     \left(
#'       \cdot
#'     \right)
#'   }
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \hat{\theta}_{
#'     \left(
#'       i
#'     \right)
#'   }
#' }
#' is the jackknife mean.
#' See
#' [`jack()`]
#' and
#' [`jack_hat()`] .
#'
#' Using
#' \eqn{\hat{z}_{0}}
#' and
#' \eqn{\hat{a}}
#' we can obtain the adjusted \eqn{z}-scores as follows
#'
#' \deqn{
#'   z_{
#'     \mathrm{BCa}_{\mathrm{lo}}
#'   }
#'   =
#'   \hat{z}_{0}
#'   +
#'   \frac{
#'     \hat{z}_{0}
#'     +
#'     z_{
#'       \left(
#'         \frac{
#'           \alpha
#'         }
#'         {
#'           2
#'         }
#'       \right)
#'     }
#'   }
#'   {
#'     1
#'     -
#'     \hat{a}
#'     \left[
#'       \hat{z}_{0}
#'       +
#'       z_{
#'         \left(
#'           \frac{
#'             \alpha
#'           }
#'           {
#'             2
#'           }
#'         \right)
#'       }
#'     \right]
#'   } ,
#' }
#'
#' \deqn{
#'   z_{
#'     \mathrm{BCa}_{\mathrm{up}}
#'   }
#'   =
#'   \hat{z}_{0}
#'   +
#'   \frac{
#'     \hat{z}_{0}
#'     +
#'     z_{
#'       \left[
#'         1
#'         -
#'         \frac{
#'           \alpha
#'         }
#'         {
#'           2
#'         }
#'       \right]
#'     }
#'   }
#'   {
#'     1
#'     -
#'     \hat{a}
#'     \left[
#'       \hat{z}_{0}
#'       +
#'       z_{
#'         \left(
#'           1
#'           -
#'           \frac{
#'             \alpha
#'           }
#'           {
#'             2
#'           }
#'         \right)
#'       }
#'     \right]
#'   } .
#' }
#'
#' The adjusted \eqn{z}-scores
#' are used to determine the adjusted
#' percentile ranks
#' to form the confidence interval.
#'
#' The bias-corrected and accelerated confidence interval is given by
#' \deqn{
#'   \left[
#'     \hat{\theta}_{\mathrm{lo}},
#'     \hat{\theta}_{\mathrm{up}}
#'   \right]
#'   =
#'   \left[
#'     \hat{\theta}^{*}_{
#'       \left(
#'         z_{
#'           \mathrm{BCa}_{
#'             \mathrm{lo}
#'           }
#'         }
#'       \right)
#'     },
#'     \hat{\theta}^{*}_{
#'       \left(
#'         z_{
#'           \mathrm{BCa}_{
#'             \mathrm{up}
#'           }
#'         }
#'       \right)
#'     }
#'   \right] .
#' }
#'
#' For more details and examples see the following vignettes:
#'
#' [Notes: Introduction to Nonparametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Introduction to Parametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat_star_jack Numeric vector.
#'   Jackknife sampling distribution,
#'   that is,
#'   the sampling distribution of `thetahat`
#'   estimated for each `i` jackknife sample.
#'   \eqn{
#'     \hat{\theta}_{\left( 1 \right)},
#'     \hat{\theta}_{\left( 2 \right)},
#'     \hat{\theta}_{\left( 3 \right)},
#'     \dots,
#'     \hat{\theta}_{\left( n \right)}
#'   } .
#'   If not provided,
#'   the jackknife sampling distribution is generated
#'   using `x_star_jack` and `fitFUN`
#'   if `x_star_jack` is provided.
#' @param x_star_jack Vector, matrix, or data frame.
#'   Jackknife sample.
#'   If not provided,
#'   the jackknife sample is generated
#'   using `data`, and the [`jack()`] function.
#'   Ignored if `thetahat_star_jack`
#'   is provided.
#' @param data Vector, matrix, or data frame.
#'   Sample data.
#'   Ignored if `thetahat_star_jack`
#'   is provided.
#' @param fitFUN Function.
#'   Fit function to use on `data`.
#'   The first argument should correspond to `data`.
#'   Other arguments can be passed to `fitFUN`
#'   using `...`.
#'   `fitFUN` should return a single value.
#'   Ignored if `thetahat_star_jack`
#'   is provided.
#' @param ... Arguments to pass to `fitFUN`.
#' @inheritParams bc
#' @inheritParams jack
#' @inheritParams jack_hat
#' @return Returns a vector with the following elements:
#'   \describe{
#'     \item{statistic}{Square root of Wald test statistic. `NA` if `wald = FALSE`.}
#'     \item{p}{p-value. `NA` if `wald = FALSE`.}
#'     \item{se}{Estimated bootstrap standard error \eqn{\left( \widehat{\mathrm{se}}_{\mathrm{B}} \left( \hat{\theta} \right) \right)}.}
#'     \item{ci_}{Estimated bias-corrected and accelerated confidence limits corresponding to alpha from the bootstrap sampling distribution thetahat_star \eqn{\left( \boldsymbol{\hat{\theta}^{*}} \right)}.}
#'   }
#' If `eval = TRUE`,
#' appends the following to the results vector
#'   \describe{
#'     \item{zero_hit_}{Logical. Tests if confidence interval contains zero.}
#'     \item{theta_hit_}{Logical. Tests if confidence interval contains theta.}
#'     \item{length_}{Length of confidence interval.}
#'     \item{shape_}{Shape of confidence interval.}
#'   }
#' @inherit pc references
#' @export
#' @family bootstrap confidence interval functions
#' @keywords confidence interval
.bca <- function(thetahat_star,
                 thetahat_star_jack = NULL,
                 x_star_jack = NULL,
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
                 par = FALSE,
                 ncores = NULL,
                 ...) {
  # standard error
  sehat_B_thetahat <- sd(thetahat_star)
  # bias-correction
  z0hat <- qnorm(
    sum(thetahat_star < thetahat) / length(thetahat_star)
  )
  # acceleration
  ## generate thetahat_star_jack if not provided
  if (is.null(thetahat_star_jack)) {
    # generate x_star_jack if not provided
    if (is.null(x_star_jack)) {
      x_star_jack <- jack(
        data = data,
        par = par,
        ncores = ncores
      )
    }
    ## generate thetahat_star_jack from generated x_star_jack or argument
    thetahat_star_jack <- util_lapply(
      FUN = fitFUN,
      args = list(
        data = x_star_jack
      ),
      par = par,
      ncores = ncores
    )
    thetahat_star_jack <- as.vector(
      do.call(
        what = "rbind",
        args = thetahat_star_jack
      )
    )
  }
  parenthesis <- mean(thetahat_star_jack) - thetahat_star_jack
  numerator <- sum(parenthesis^3)
  denominator <- 6 * ((sum(parenthesis^2))^(3 / 2))
  ahat <- numerator / denominator
  # confidence interval
  probs <- alpha2prob(alpha = alpha)
  z1 <- qnorm(
    p = probs
  )
  bca_probs <- pnorm(
    z0hat + (z0hat + z1) / (1 - ahat * (z0hat + z1))
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
  # optional - wald
  if (wald) {
    sqrt_W <- sqrt_wald_test(
      thetahat = thetahat,
      sehat_thetahat = sehat_B_thetahat,
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
  # initial output
  out <- c(
    sqrt_W,
    se = sehat_B_thetahat,
    ci
  )
  # optional - ci evaluation
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

#' Confidence Interval - Bias-Corrected and Accelerated
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams .bca
#' @inherit .bca description details return references
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
                par = FALSE,
                ncores = NULL,
                ...) {
  .bca(
    thetahat_star = thetahat_star,
    thetahat_star_jack = NULL,
    x_star_jack = NULL,
    thetahat = thetahat,
    data = data,
    fitFUN = fitFUN,
    alpha = alpha,
    wald = wald,
    theta_null = theta_null,
    distribution = distribution,
    df = df,
    eval = eval,
    theta = theta,
    par = par,
    ncores = ncores,
    ...
  )
}
