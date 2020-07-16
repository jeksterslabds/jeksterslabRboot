#' Jackknife
#'
#' Generates jackknife samples.
#'
#' \deqn{
#'   \mathbf{
#'     x
#'   }_i
#'   =
#'   \left\{
#'     x_1,
#'     x_2,
#'     x_3,
#'     \dots,
#'     x_n
#'   \right\}
#' }
#' for
#' \eqn{
#'   i
#'   =
#'   \left\{
#'     1,
#'     2,
#'     3,
#'     \dots,
#'     n
#'   \right\}
#' }
#' the \eqn{i}th jackknife sample
#' consists of the original sample data
#' with the \eqn{i}th observation removed.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family jackknife functions
#' @keywords jackknife
#' @inheritParams nb
#' @param data Vector, matrix, or data frame.
#'   Sample data.
#' @return Returns a list of jackknife samples
#'   of length `n`,
#'   where `n` is the sample size of `data`.
#'   Each item in the list,
#'   will have a sample size of
#'   `n - 1`.
#' @examples
#' n <- 5
#' #################################
#' # vector
#' #################################
#' x <- rnorm(n = n)
#' xstar <- jack(
#'   data = x
#' )
#' str(xstar)
#' #################################
#' # matrix
#' #################################
#' x1 <- rnorm(n = n)
#' x2 <- rnorm(n = n)
#' x3 <- rnorm(n = n)
#' X <- cbind(x1, x2, x3)
#' Xstar <- jack(
#'   data = X
#' )
#' str(Xstar)
#' #################################
#' # data frame
#' #################################
#' X <- as.data.frame(X)
#' Xstar <- jack(
#'   data = X
#' )
#' str(Xstar)
#' @references
#' [Wikipedia: Jackknife resampling](https://en.wikipedia.org/wiki/Jackknife_resampling)
#' @export
jack <- function(data,
                 par = FALSE,
                 ncores = NULL,
                 mc = TRUE,
                 lb = FALSE,
                 cl_eval = FALSE,
                 cl_export = FALSE,
                 cl_expr,
                 cl_vars) {
  if (is.vector(data)) {
    n <- length(data)
  }
  if (is.matrix(data) | is.data.frame(data)) {
    n <- nrow(data)
  }
  foo <- function(i, data) {
    if (is.data.frame(data) | is.matrix(data)) {
      return(data[-i, ])
    }
    if (is.vector(data)) {
      return(data[-i])
    }
  }
  par_lapply(
    X = 1:n,
    FUN = foo,
    data = data,
    par = par,
    ncores = ncores,
    mc = mc,
    lb = lb,
    cl_eval = cl_eval,
    cl_export = cl_export,
    cl_expr = cl_expr,
    cl_vars = cl_vars,
    rbind = NULL
  )
}

#' Jackknife Estimates
#'
#' Calculates jackknife estimates.
#'
#' The jackknife estimate of bias is given by
#' \deqn{
#'   \widehat{
#'     \mathrm{
#'       bias
#'     }
#'   }_{
#'     \mathrm{
#'       jack
#'     }
#'   }
#'   \left(
#'     \theta
#'   \right)
#'   =
#'   \left(
#'     n
#'     -
#'     1
#'   \right)
#'   \left(
#'     \hat{
#'       \theta
#'     }_{
#'       \left(
#'         \cdot
#'       \right)
#'     }
#'     -
#'     \hat{
#'       \theta
#'     }
#'   \right)
#' }
#'
#' where
#'
#' \deqn{
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       \cdot
#'     \right)
#'   }
#'   =
#'   \frac{
#'     1
#'   }
#'   {
#'     n
#'   }
#'   \sum_{
#'     i
#'     =
#'     1
#'   }^{
#'     n
#'   }
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       i
#'     \right)
#'   } .
#' }
#'
#' The jackknife estimate of standard error
#' is given by
#'
#' \deqn{
#'   \widehat{
#'     \mathrm{
#'       se
#'     }
#'   }_{
#'     \mathrm{
#'       jack
#'     }
#'   }
#'   \left(
#'     \hat{
#'       \theta
#'     }
#'   \right)
#'   =
#'   \sqrt{
#'     \frac{
#'       n
#'       -
#'       1
#'     }
#'     {
#'       n
#'     }
#'     \sum_{
#'       i
#'       =
#'       1
#'     }^{
#'       n
#'     }
#'     \left(
#'       \hat{
#'         \theta
#'       }_{
#'         \left(
#'           i
#'         \right)
#'       }
#'       -
#'       \hat{
#'         \theta
#'       }_{
#'         \left(
#'           \cdot
#'         \right)
#'       }
#'     \right)^2
#'   } .
#' }
#'
#' The bias-corrected jackknife estimate
#' is given by
#'
#' \deqn{
#'   \hat{
#'     \theta
#'   }_{
#'     \mathrm{
#'       jack
#'     }
#'   }
#'   =
#'   \hat{
#'     \theta
#'   }
#'   -
#'   \hat{
#'     \mathrm{
#'       bias
#'     }
#'   }_{
#'     \mathrm{
#'       jack
#'     }
#'   }
#'   \left(
#'     \theta
#'   \right)
#'   =
#'   n
#'   \hat{
#'     \theta
#'   }
#'   -
#'   \left(
#'     n
#'     -
#'     1
#'   \right)
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       \cdot
#'     \right)
#'   } .
#' }
#'
#' Pseudo-values can be computed using
#'
#' \deqn{
#'   \tilde{
#'     \theta
#'   }_{
#'     i
#'   }
#'   =
#'   n
#'   \hat{
#'     \theta
#'   }
#'   -
#'   \left(
#'     n
#'     -
#'     1
#'   \right)
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       i
#'     \right)
#'   } .
#' }
#'
#' The standard error can be estimated using the pseudo-values
#'
#' \deqn{
#'   \widehat{
#'     \mathrm{
#'       se
#'     }
#'   }_{
#'     \mathrm{
#'       jack
#'     }
#'   }
#'   \left(
#'     \tilde{
#'       \theta
#'     }
#'   \right)
#'   =
#'   \sqrt{
#'     \sum_{
#'       i
#'       =
#'       1
#'     }^{
#'       n
#'     }
#'     \frac{
#'       \left(
#'         \tilde{
#'           \theta
#'         }_{
#'           i
#'         }
#'         -
#'         \tilde{
#'           \theta
#'         }
#'       \right)^2
#'     }
#'     {
#'       \left(
#'         n
#'         -
#'         1
#'       \right)
#'       n
#'     }
#'   }
#' }
#'
#' where
#'
#' \deqn{
#'   \tilde{
#'     \theta
#'   }
#'   =
#'   \frac{
#'     1
#'   }
#'   {
#'     n
#'   }
#'   \sum_{
#'     i
#'     =
#'     1
#'   }^{
#'     n
#'   }
#'   \tilde{
#'     \theta
#'   }_{
#'     i
#'   } .
#' }
#'
#' An interval can be generated using
#'
#' \deqn{
#'   \tilde{
#'     \theta
#'   }
#'   \pm
#'   t_{
#'     \frac{
#'       \alpha
#'     }
#'     {
#'       2
#'     }
#'   }
#'   \times
#'   \widehat{
#'     \mathrm{
#'       se
#'     }
#'   }_{
#'     \mathrm{
#'       jack
#'     }
#'   }
#'   \left(
#'     \tilde{
#'       \theta
#'     }
#'   \right)
#' }
#' with degrees of freedom
#' \eqn{
#'   \nu
#'   =
#'   n
#'   -
#'   1
#' } .
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family jackknife functions
#' @keywords jackknife
#' @inheritParams wald
#' @inherit jack references
#' @param thetahatstarjack Numeric vector.
#' Jackknife sampling distribution,
#' that is,
#' the sampling distribution of `thetahat`
#' estimated for each `i` jackknife sample.
#' \eqn{
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       1
#'     \right)
#'   },
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       2
#'     \right)
#'   },
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       3
#'     \right)
#'   },
#'   \dots,
#'   \hat{
#'     \theta
#'   }_{
#'     \left(
#'       n
#'     \right)
#'   }
#' } .
#' @param thetahat Numeric.
#' Parameter estimate
#' \eqn{
#'   \left(
#'     \hat{
#'       \theta
#'     }
#'   \right)
#' }
#' from the original sample data.
#' @return Returns a list with the following elements:
#'   \describe{
#'     \item{hat}{Jackknife estimates.}
#'     \item{ps}{Pseudo-values.}
#'     \item{ci}{Confidence intervals using pseudo-values.}
#'   }
#' The first list element `hat` contains the following:
#'   \describe{
#'     \item{mean}{Mean of `thetahatstarjack` \eqn{\left( \hat{\theta}_{\left( \cdot \right) } \right)}.}
#'     \item{bias}{Jackknife estimate of bias \eqn{\left( \widehat{\mathrm{bias}}_{\mathrm{jack}} \left( \theta \right) \right)}.}
#'     \item{se}{Jackknife estimate of standard error \eqn{\left( \widehat{\mathrm{se}}_{\mathrm{jack}} \left( \hat{\theta} \right) \right)}.}
#'     \item{thetahatjack}{Bias-corrected jackknife estimate \eqn{\left( \hat{\theta}_{\mathrm{jack}} \right)}.}
#'   }
#' @examples
#' n <- 100
#' x <- rnorm(n = n)
#' thetahat <- mean(x)
#' xstar <- jack(
#'   data = x
#' )
#' thetahatstarjack <- sapply(
#'   X = xstar,
#'   FUN = mean
#' )
#' str(xstar, list.len = 6)
#' hist(
#'   thetahatstarjack,
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
#' jack_hat(
#'   thetahatstarjack = thetahatstarjack,
#'   thetahat = thetahat
#' )
#' @export
jack_hat <- function(thetahatstarjack,
                     thetahat,
                     alpha = c(
                       0.001,
                       0.01,
                       0.05
                     ),
                     eval = FALSE,
                     theta = 0) {
  n <- length(thetahatstarjack)
  mean_thetahat <- mean(thetahatstarjack)
  bias <- (n - 1) * (mean_thetahat - thetahat)
  se <- sqrt(((n - 1) / n) * sum(thetahatstarjack - mean_thetahat)^2)
  thetahatjack <- mean_thetahat - bias
  # pseudo-values
  pseudo_values <- n * thetahat - (n - 1) * thetahatstarjack
  mean_pseudo_values <- mean(pseudo_values)
  # research more on formula for se_pseudo_values
  se_pseudo_values <- sqrt(sum(((pseudo_values - mean_pseudo_values)^2) / ((n - 1) * n)))
  ci <- wald(
    thetahat = mean_pseudo_values,
    sehat = se_pseudo_values,
    null = 0,
    alpha = alpha,
    dist = "t",
    df = n - 1,
    eval = eval,
    theta = theta
  )
  ci[1] <- NA
  ci[2] <- NA
  hat <- c(
    mean = mean_thetahat,
    bias = bias,
    se = se,
    thetahatjack = thetahatjack,
    mean_ps = mean_pseudo_values,
    se_ps = se_pseudo_values
  )
  list(
    hat = hat,
    ps = pseudo_values,
    ci = ci
  )
}
