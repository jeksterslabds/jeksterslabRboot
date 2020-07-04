#' Monte Carlo
#'
#' Generates `R` number of parameter estimates
#' from a vector of `thetahat`
#' \eqn{\left( \boldsymbol{\hat{\theta} }\right)}
#' and
#' sampling variances and covariances of `thetahat`
#' \eqn{\hat{\mathrm{Cov}} \left( \boldsymbol{\hat{\theta}} \right)}.
#' If length of `thetahat` is greater than 1,
#' a multivariate normal distribution is assumed.
#' If length of `thetahat` is equal 1,
#' a univariate normal distribution is assumed.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat Numeric vector.
#'   Parameter estimates
#'   \eqn{\boldsymbol{\hat{\theta}}}.
#' @param covhat_thetahat Symmetric matrix.
#'   Estimated sampling variances and covariances of
#'   parameter estimates
#'   \eqn{\hat{\mathrm{Cov}} \left( \boldsymbol{\hat{\theta}} \right)}.
#' @param R Integer.
#'   Number of Monte Carlo replications.
#' @examples
#' R <- 20000L
#' # length 1
#' thetahat <- 100
#' covhat_thetahat <- 1.5
#' mc_star_length_1 <- mc(
#'   thetahat = thetahat,
#'   covhat_thetahat = covhat_thetahat,
#'   R = R
#' )
#' head(mc_star_length_1)
#' hist(
#'   mc_star_length_1,
#'   main = expression(
#'     paste(
#'       "Histogram of ",
#'       hat(theta),
#'       "*"
#'     )
#'   ),
#'   xlab = expression(
#'     paste(
#'       hat(theta),
#'       "*"
#'     )
#'   )
#' )
#' qqnorm(mc_star_length_1)
#' qqline(mc_star_length_1)
#' # length greater than 1
#' alphahat <- 0.3386
#' betahat <- 0.4510
#' alphahat_betahat <- alphahat * betahat
#' varhat_alphahat <- 0.1224^2
#' varhat_betahat <- 0.1460^2
#' thetahat <- c(
#'   alphahat,
#'   betahat
#' )
#' covhat_thetahat <- matrix(
#'   data = c(
#'     varhat_alphahat,
#'     0.00,
#'     0.00,
#'     varhat_betahat
#'   ),
#'   ncol = 2
#' )
#' mc_star_length_2 <- mc(
#'   thetahat = thetahat,
#'   covhat_thetahat = covhat_thetahat,
#'   R = R
#' )
#' head(mc_star_length_2)
#' alphahat_betahat_star <- mc_star_length_2[, 1] * mc_star_length_2[, 2]
#' hist(
#'   alphahat_betahat_star,
#'   main = expression(
#'     paste(
#'       "Histogram of ",
#'       hat(alpha),
#'       hat(beta),
#'       "*"
#'     )
#'   ),
#'   xlab = expression(
#'     paste(
#'       hat(alpha),
#'       hat(beta),
#'       "*"
#'     )
#'   )
#' )
#' qqnorm(alphahat_betahat_star)
#' qqline(alphahat_betahat_star)
#' wald(
#'   thetahat = alphahat_betahat,
#'   sehat_thetahat = sd(alphahat_betahat_star)
#' )
#' @export
mc <- function(thetahat,
               covhat_thetahat,
               R = 20000L) {
  if (length(thetahat) == 1) {
    # if (length(covhat_thetahat) != 1) {
    #  stop(
    #    "Dimensions of thetahat and covhat_thetahat are incompatible."
    #  )
    # }
    return(
      rnorm(
        n = R,
        mean = thetahat,
        sd = sqrt(covhat_thetahat)
      )
    )
  } else {
    # dimensions <- c(
    #  length(thetahat),
    #  nrow(covhat_thetahat),
    #  ncol(covhat_thetahat)
    # )
    # if (mean(dimensions) != length(thetahat)) {
    #  stop(
    #    "Dimensions of thetahat and covhat_thetahat are incompatible."
    #  )
    # }
    # symmetric <- is.symmetric(
    #  X = covhat_thetahat,
    #  stop = FALSE
    # )
    # if (!symmetric) {
    #  stop(
    #    "covhat_thetahat is not a symmetric matrix."
    #  )
    # }
    return(
      mvrnorm(
        n = R,
        Sigma = covhat_thetahat,
        mu = thetahat
      )
    )
  }
}
