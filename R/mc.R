#' Monte Carlo
#'
#' Generates `R` number of parameter estimates
#' from a vector of `thetahat`
#' \eqn{\left( \boldsymbol{\hat{\theta} }\right)}
#' and
#' sampling variances and covariances of `thetahat`
#' \eqn{\widehat{\mathrm{Cov}} \left( \boldsymbol{\hat{\theta}} \right)}.
#' If length of `thetahat` is greater than 1,
#' a multivariate normal distribution is assumed.
#' If length of `thetahat` is equal 1,
#' a univariate normal distribution is assumed.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param thetahat Numeric vector.
#'   Parameter estimates
#'   \eqn{\boldsymbol{\hat{\theta}}}.
#' @param vcovhat Symmetric matrix.
#'   Estimated sampling variances and covariances of
#'   parameter estimates
#'   \eqn{\widehat{\mathrm{Cov}} \left( \boldsymbol{\hat{\theta}} \right)}.
#' @param R Integer.
#'   Number of Monte Carlo replications.
#' @examples
#' R <- 20000L
#' # length 1
#' thetahat <- 100
#' vcovhat <- 1.5
#' mc_star_length_1 <- mc(
#'   thetahat = thetahat,
#'   vcovhat = vcovhat,
#'   R = R
#' )
#' str(mc_star_length_1)
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
#' vcovhat <- matrix(
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
#'   vcovhat = vcovhat,
#'   R = R
#' )
#' str(mc_star_length_2)
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
#'   sehat = sd(alphahat_betahat_star)
#' )
#' @export
mc <- function(thetahat,
               vcovhat,
               R = 20000L) {
  if (length(thetahat) == 1) {
    # if (length(vcovhat) != 1) {
    #  stop(
    #    "Dimensions of thetahat and vcovhat are incompatible."
    #  )
    # }
    return(
      rnorm(
        n = R,
        mean = thetahat,
        sd = sqrt(vcovhat)
      )
    )
  } else {
    # dimensions <- c(
    #  length(thetahat),
    #  nrow(vcovhat),
    #  ncol(vcovhat)
    # )
    # if (mean(dimensions) != length(thetahat)) {
    #  stop(
    #    "Dimensions of thetahat and vcovhat are incompatible."
    #  )
    # }
    # symmetric <- is.symmetric(
    #  X = vcovhat,
    #  stop = FALSE
    # )
    # if (!symmetric) {
    #  stop(
    #    "vcovhat is not a symmetric matrix."
    #  )
    # }
    return(
      mvrnorm(
        n = R,
        Sigma = vcovhat,
        mu = thetahat
      )
    )
  }
}
