#' Nonparametric Bootstrap
#'
#' Generates `B` number of nonparametric bootstrap
#'   samples from the original sample data.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param data Vector, matrix or data frame.
#'   Sample data to bootstrap.
#' @param B Integer. Number of bootstrap samples.
#' @inheritParams jeksterslabRutils::util_lapply
#' @return Returns a list of nonparametric bootstrap samples.
#' @family bootstrap functions
#' @keywords bootstraping
#' @examples
#' B <- 5L
#' # vector
#' data_vector <- rnorm(n = 5)
#' nb(
#'   data = data_vector,
#'   B = B
#' )
#' # matrix
#' X <- rnorm(n = 5)
#' M <- rnorm(n = 5)
#' Y <- rnorm(n = 5)
#' data_matrix <- cbind(X, M, Y)
#' nb(
#'   data = data_matrix,
#'   B = B
#' )
#' # data frame
#' data_dataframe <- data.frame(X, M, Y)
#' nb(
#'   data = data_dataframe,
#'   B = B
#' )
#' @references
#' [Wikipedia: Bootstrapping (statistics)](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
#' @importFrom jeksterslabRutils util_lapply
#' @export
nb <- function(data,
               B = 2000L,
               par = FALSE,
               ncores = NULL) {
  if (is.vector(data)) {
    len <- 1:length(data)
  }
  if (is.data.frame(data) | is.matrix(data)) {
    len <- 1:nrow(data)
  }
  foo <- function(x, data) {
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
  util_lapply(
    FUN = foo,
    args = list(
      x = 1:B,
      data = data
    ),
    par = par,
    ncores = ncores
  )
}

#' Parametric Bootstrap
#' (Multivariate Normal)
#' from
#' \eqn{\boldsymbol{\hat{\Sigma}}}
#' and
#' \eqn{\boldsymbol{\hat{\mu}}}
#'
#' Generates `B` number of parametric bootstrap
#'   samples from the original sample data.
#'   Data is generated from a multivariate normal distribution
#'   using the estimated variance-covariance matrix and mean vector.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param n Integer.
#'   Sample size.
#' @param Sigmahat Matrix.
#'   Estimated variance-covariance matrix
#'   from the original sample data.
#' @param muhat Vector.
#'   Estimated mean vector from the original sample data.
#' @inheritParams nb
#' @return Returns a list of parametric bootstrap samples.
#' @importFrom MASS mvrnorm
#' @family bootstrap functions
#' @keywords bootstraping
#' @examples
#' B <- 5L
#' Sigmahat <- matrix(
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
#' muhat <- c(
#'   108.3060,
#'   105.3324,
#'   103.4009
#' )
#' .pb_mvn(
#'   n = 5,
#'   Sigmahat = Sigmahat,
#'   muhat = muhat,
#'   B = B
#' )
#' @inherit nb references
#' @export
.pb_mvn <- function(n,
                    Sigmahat,
                    muhat,
                    B = 2000L,
                    par = FALSE,
                    ncores = NULL) {
  foo <- function(x,
                  n,
                  Sigmahat,
                  muhat) {
    mvrnorm(
      n = n,
      Sigma = Sigmahat,
      mu = muhat
    )
  }
  util_lapply(
    FUN = foo,
    args = list(
      x = 1:B,
      n = n,
      Sigmahat = Sigmahat,
      muhat = muhat
    ),
    par = par,
    ncores = ncores
  )
}

#' Parametric Bootstrap
#' (Multivariate Normal)
#' from Data
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams nb
#' @examples
#' data <- matrix(
#'   data = c(
#'     121.03482,
#'     99.74874,
#'     101.44951,
#'     114.31910,
#'     104.97803,
#'     118.31483,
#'     101.39711,
#'     90.20035,
#'     105.31186,
#'     111.43784,
#'     95.68517,
#'     109.52585,
#'     105.11345,
#'     128.37871,
#'     78.30152
#'   ),
#'   nrow = 5
#' )
#' pb_mvn(
#'   data = data,
#'   B = 5L
#' )
#' @inherit .pb_mvn references description return
#' @importFrom stats cov
#' @export
pb_mvn <- function(data,
                   B = 2000L,
                   par = FALSE,
                   ncores = NULL) {
  n <- nrow(data)
  Sigmahat <- cov(data)
  muhat <- colMeans(data)
  .pb_mvn(
    n = n,
    Sigmahat = Sigmahat,
    muhat = muhat,
    B = B,
    par = par,
    ncores = ncores
  )
}

#' Parametric Bootstrap (Univariate)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param rFUN Function.
#'   Data generating function to generate univariate data.
#' @param n Integer.
#'   Sample size.
#' @param ... Arguments to pass to rFUN.
#' @inheritParams nb
#' @importFrom stats rnorm
#' @examples
#' pb_univ(
#'   rFUN = rexp,
#'   n = 5,
#'   B = 5,
#'   rate = 1
#' )
#' @export
pb_univ <- function(rFUN = rnorm,
                    n,
                    B = 2000L,
                    par = FALSE,
                    ncores = NULL,
                    ...) {
  args <- list(
    x = 1:B,
    rFUN = rFUN,
    n = n,
    ...
  )
  foo <- function(x,
                  rFUN,
                  n,
                  ...) {
    rFUN(
      n = n,
      ...
    )
  }
  util_lapply(
    FUN = foo,
    args = args,
    par = par,
    ncores = ncores
  )
}
