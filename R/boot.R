#' Nonparametric Bootstrap
#'
#' @description Generates `B` number of nonparametric bootstrap
#' samples from the original sample `data`.
#' `data` is referred to as the empirical distribution
#' \eqn{
#'   \hat{
#'     F
#'   }
#'   %(\#eq:boot-ecdf)
#' }.
#'
#' @details For more details and examples see the following vignettes:
#'
#' [Notes: Intro to NB](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Intro to PB](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams jeksterslabRpar::par_lapply
#' @param data Vector, matrix or data frame.
#'   Sample data to bootstrap.
#'   The empirical distribution
#'   \eqn{
#'     \hat{
#'       F
#'     }
#'     %(\#eq:boot-ecdf)
#'   }.
#' @param B Integer.
#'   Number of bootstrap samples.
#' @return Returns a list of length `B` of nonparametric bootstrap samples.
#' @examples
#' B <- 5L
#' n <- 5
#' # vector----------------------------------------------------------------------
#' x <- rnorm(n = n)
#' xstar <- nb(
#'   data = x,
#'   B = B
#' )
#' str(xstar)
#' # matrix----------------------------------------------------------------------
#' x1 <- rnorm(n = n)
#' x2 <- rnorm(n = n)
#' x3 <- rnorm(n = n)
#' X <- cbind(x1, x2, x3)
#' Xstar <- nb(
#'   data = X,
#'   B = B
#' )
#' str(Xstar)
#' # data frame------------------------------------------------------------------
#' X <- as.data.frame(X)
#' Xstar <- nb(
#'   data = X,
#'   B = B
#' )
#' str(Xstar)
#' @references
#'   Efron, B., & Tibshirani, R. J. (1993).
#'   *An introduction to the bootstrap*.
#'   New York, N.Y: Chapman & Hall.
#'
#'   [Wikipedia: Bootstrapping (statistics)](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
#' @importFrom jeksterslabRpar par_lapply
#' @export
nb <- function(data,
               B = 2000L,
               # par_lapply
               par = FALSE,
               ncores = NULL,
               mc = TRUE,
               lb = FALSE,
               cl_eval = FALSE,
               cl_export = FALSE,
               cl_expr,
               cl_vars) {
  if (is.vector(data)) {
    len <- 1:length(data)
  }
  if (is.data.frame(data) | is.matrix(data)) {
    len <- 1:nrow(data)
  }
  foo <- function(iter,
                  data) {
    i <- sample(
      x = len,
      replace = TRUE
    )
    if (is.data.frame(data) | is.matrix(data)) {
      data <- data[i, ]
      rownames(data) <- len
      return(data)
    }
    if (is.vector(data)) {
      return(data[i])
    }
  }
  par_lapply(
    X = 1:B,
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

#' Parametric Bootstrap
#' (Multivariate Normal)
#' from
#' \eqn{
#'   \boldsymbol{
#'     \hat{
#'       \mu
#'     }
#'   }
#'   \left(
#'     \boldsymbol{
#'       \hat{
#'         \theta
#'       }
#'     }
#'   \right)
#'   %(\#eq:boot-pb-mvn-mu)
#' }
#' and
#' \eqn{
#'   \boldsymbol{
#'     \hat{
#'       \Sigma
#'     }
#'   }
#'   \left(
#'     \boldsymbol{
#'       \hat{
#'         \theta
#'       }
#'     }
#'   \right)
#'   %(\#eq:boot-pb-mvn-Sigma)
#' } .
#'
#' @description Generates `B` number of parametric bootstrap
#' samples using estimated parameters
#' from the original sample `data`.
#' `data` is referred to as the empirical distribution
#' with the following distributional assumption
#' \deqn{
#'   \hat{
#'     F
#'   }_{
#'     \mathcal{
#'       N
#'     }_{k}
#'     \left(
#'       \boldsymbol{
#'         \hat{
#'           \mu
#'         }
#'       }
#'       \left(
#'         \boldsymbol{
#'           \hat{
#'             \theta
#'           }
#'         }
#'       \right) ,
#'       \boldsymbol{
#'         \hat{
#'           \Sigma
#'         }
#'       }
#'       \left(
#'         \boldsymbol{
#'           \hat{
#'             \theta
#'           }
#'         }
#'       \right)
#'     \right)
#'   } .
#'   %(\#eq:boot-pb-mvn)
#' }
#' Bootstrap samples are generated from a multivariate normal distribution
#' using the fitted model-implied mean vector
#' and variance-covariance matrix.
#'
#' @details For more details and examples see the following vignettes:
#'
#' [Notes: Intro to NB](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Intro to PB](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams nb
#' @inheritParams jeksterslabRdata::mvn
#' @inherit nb references
#' @param n Integer.
#'   Sample size.
#' @param muhatthetahat Vector.
#'   Mean vector as a function of estimated parameters
#'   or the fitted model-implied mean vector
#'   \eqn{
#'     \boldsymbol{
#'       \hat{
#'         \mu
#'       }
#'     }
#'     \left(
#'       \boldsymbol{
#'         \hat{
#'           \theta
#'         }
#'       }
#'     \right)
#'     %(\#eq:boot-pb-mvn-mu)
#'   } .
#' @param Sigmahatthetahat Matrix.
#'   Variance-covariance matrix as a function of estimated parameters
#'   or the fitted model-implied variance-covariance matrix
#'   \eqn{
#'     \boldsymbol{
#'       \hat{
#'         \Sigma
#'       }
#'     }
#'     \left(
#'       \boldsymbol{
#'         \hat{
#'           \theta
#'         }
#'       }
#'     \right)
#'     %(\#eq:boot-pb-mvn-Sigma)
#'   } .
#' @return Returns a list of length `B` of parametric bootstrap samples.
#' @examples
#' B <- 5L
#' Sigmahatthetahat <- matrix(
#'   data = c(
#'     82.37344,
#'     70.55922,
#'     17.83930,
#'     70.55922,
#'     112.57145,
#'     -75.98558,
#'     17.83930,
#'     -75.98558,
#'     338.46263
#'   ),
#'   nrow = 3
#' )
#' muhatthetahat <- c(
#'   108.3060,
#'   105.3324,
#'   103.4009
#' )
#' Xstar <- pbmvn(
#'   n = 5,
#'   Sigmahatthetahat = Sigmahatthetahat,
#'   muhatthetahat = muhatthetahat,
#'   B = B
#' )
#' str(Xstar)
#' @importFrom jeksterslabRdata mvn
#' @export
pbmvn <- function(n, # mvrnorm
                  muhatthetahat,
                  Sigmahatthetahat,
                  tol = 1e-6,
                  empirical = FALSE,
                  # iter
                  B = 2000L,
                  # par_lapply
                  par = FALSE,
                  ncores = NULL,
                  mc = TRUE,
                  lb = FALSE,
                  cl_eval = FALSE,
                  cl_export = FALSE,
                  cl_expr,
                  cl_vars) {
  mvn(
    # mvrnorm
    n = n,
    mu = muhatthetahat,
    Sigma = Sigmahatthetahat,
    tol = tol,
    empirical = empirical,
    # iter
    R = B,
    # par_lapply
    par = par,
    ncores = ncores,
    mc = mc,
    lb = lb,
    cl_eval = cl_eval,
    cl_export = cl_export,
    cl_expr = cl_expr,
    cl_vars = cl_vars
  )
}

#' Parametric Bootstrap (Univariate)
#'
#' @description Generates `B` number of parametric bootstrap
#' samples from estimated parameters of original univariate sample `data`.
#' `data` is referred to as the empirical distribution
#' \deqn{
#'   \hat{
#'     F
#'   }_{
#'     \mathcal{
#'       Distribution
#'     }
#'   }
#' }.
#' The default distribution is
#' \deqn{
#'   \mathcal{
#'     N
#'   }
#'   \left(
#'     \hat{
#'       \mu
#'     },
#'     \hat{
#'       \sigma
#'     }^2
#'   \right)
#'   %(\#eq:boot-pb-norm)
#' }.
#' The univariate distribution and parameters used in the
#' data generating process can be specified using `rFUN` and `...`.
#'
#' @details For more details and examples see the following vignettes:
#'
#' [Notes: Intro to NB](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Intro to PB](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams nb
#' @inheritParams jeksterslabRpar::par_lapply
#' @inherit pbmvn references return
#' @param n Integer.
#'   Sample size.
#' @param rFUN Function.
#'   Data generating function to generate univariate data.
#' @param ... Arguments to pass to `rFUN`.
#' @examples
#' n <- 5
#' B <- 5
#' # normal distribution---------------------------------------------------------
#' mu <- 100
#' sigma2 <- 225
#' sigma <- sqrt(sigma2)
#' x <- rnorm(
#'   n = n,
#'   mean = mu,
#'   sd = sigma
#' )
#' muhatthetahat <- mean(x)
#' Sigmahatthetahat <- sd(x)
#' xstar <- pbuniv(
#'   n = n,
#'   rFUN = rnorm,
#'   B = B,
#'   mean = muhatthetahat,
#'   sd = Sigmahatthetahat
#' )
#' str(xstar)
#' # binomial distribution-------------------------------------------------------
#' n_trials <- 1
#' p <- 0.50
#' x <- rbinom(
#'   n = n,
#'   size = n_trials,
#'   prob = p
#' )
#' phat <- mean(x) / n_trials
#' xstar <- pbuniv(
#'   n = n,
#'   rFUN = rbinom,
#'   B = B,
#'   size = n_trials,
#'   prob = phat
#' )
#' str(xstar)
#' @importFrom jeksterslabRdata univ
#' @importFrom stats rnorm
#' @export
pbuniv <- function(n, # rFUN
                   rFUN = rnorm,
                   ...,
                   # iter
                   B = 2000L,
                   # par_lapply
                   par = FALSE,
                   ncores = NULL,
                   mc = TRUE,
                   lb = FALSE,
                   cl_eval = FALSE,
                   cl_export = FALSE,
                   cl_expr,
                   cl_vars,
                   rbind = NULL) {
  univ(
    # rFUN
    n = n,
    rFUN = rFUN,
    ...,
    # iter
    R = B,
    # par_lapply
    par = par,
    ncores = ncores,
    mc = mc,
    lb = lb,
    cl_eval = cl_eval,
    cl_export = cl_export,
    cl_expr = cl_expr,
    cl_vars = cl_vars,
    rbind = rbind
  )
}
