#' Nonparametric Bootstrap
#'
#' Generates `B` number of nonparametric bootstrap
#' samples from the original sample `data`.
#' `data` is referred to as the empirical distribution \eqn{\hat{F}}.
#'
#' For more details and examples see the following vignettes:
#'
#' [Notes: Introduction to Nonparametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Introduction to Parametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param data Vector, matrix or data frame.
#'   Sample data to bootstrap.
#'   The empirical distribution \eqn{\hat{F}}.
#' @param B Integer. Number of bootstrap samples.
#' @inheritParams jeksterslabRutils::util_lapply
#' @return Returns a list of length `B` of nonparametric bootstrap samples.
#' @examples
#' B <- 5L
#' n <- 5
#' # vector
#' data_vector <- rnorm(n = n)
#' nb(
#'   data = data_vector,
#'   B = B
#' )
#' # matrix
#' X <- rnorm(n = n)
#' Y <- rnorm(n = n)
#' Z <- rnorm(n = n)
#' data_matrix <- cbind(X, Y, Z)
#' nb(
#'   data = data_matrix,
#'   B = B
#' )
#' # data frame
#' data_dataframe <- data.frame(X, Y, Z)
#' nb(
#'   data = data_dataframe,
#'   B = B
#' )
#' @references
#' Efron, B., & Tibshirani, R. J. (1993).
#' An introduction to the bootstrap. New York, N.Y: Chapman & Hall.
#'
#' [Wikipedia: Bootstrapping (statistics)](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
#' @family bootstrap functions
#' @keywords bootstrap
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
#' samples from the original sample `data`.
#' `data` is referred to as the empirical distribution
#' \eqn{
#'   \hat{F}_{
#'     \mathrm{MVN}
#'     \left(
#'       \boldsymbol{\hat{\Sigma}},
#'       \boldsymbol{\hat{\mu}}
#'     \right)
#'   }
#' } .
#' Data is generated from a multivariate normal distribution
#' using the estimated variance-covariance matrix
#' \eqn{\boldsymbol{\hat{\Sigma}}}
#' and
#' mean vector
#' \eqn{\boldsymbol{\hat{\mu}}}
#' .
#'
#' For more details and examples see the following vignettes:
#'
#' [Notes: Introduction to Nonparametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Introduction to Parametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
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
#' @return Returns a list of length `B` of parametric bootstrap samples.
#' @importFrom MASS mvrnorm
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
#' @family bootstrap functions
#' @keywords bootstrap
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
#' @inherit .pb_mvn references description details return
#' @importFrom stats cov
#' @family bootstrap functions
#' @keywords bootstrap
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
#' Generates `B` number of parametric bootstrap
#' samples from estimated parameters of original univariate sample `data`.
#' `data` is referred to as the empirical distribution
#' \eqn{\hat{F}_{\mathcal{Distribution}}}.
#' The default distribution is
#' \eqn{
#'   \mathcal{N}
#'   \left(
#'     \hat{\mu},
#'     \hat{\sigma}^2
#'   \right)
#' }.
#' The univariate distribution and parameters used in the
#' data generating process can be specified using `rFUN` and `...`.
#'
#' For more details and examples see the following vignettes:
#'
#' [Notes: Introduction to Nonparametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_nb.html)
#'
#' [Notes: Introduction to Parametric Bootstrapping](https://jeksterslabds.github.io/jeksterslabRboot/articles/notes/notes_intro_pb.html)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param rFUN Function.
#'   Data generating function to generate univariate data.
#' @param n Integer.
#'   Sample size.
#' @param ... Arguments to pass to rFUN.
#' @inheritParams nb
#' @examples
#' # Normal distribution
#' n <- 5
#' B <- 5
#' mu <- 100
#' sigma2 <- 225
#' sigma <- sqrt(sigma2)
#' data <- rnorm(
#'   n = n,
#'   mean = mu,
#'   sd = sigma
#' )
#' muhat <- mean(data)
#' sigmahat <- sd(data)
#' pb_univ(
#'   rFUN = rnorm,
#'   n = n,
#'   B = B,
#'   mean = muhat,
#'   sd = sigmahat
#' )
#' # Binomial distribution
#' n_trials <- 1
#' p <- 0.50
#' data <- rbinom(
#'   n = n,
#'   size = n_trials,
#'   prob = p
#' )
#' phat <- mean(data) / n_trials
#' pb_univ(
#'   rFUN = rbinom,
#'   n = n,
#'   B = B,
#'   size = n_trials,
#'   prob = phat
#' )
#' @inherit .pb_mvn references return
#' @family bootstrap functions
#' @keywords bootstrap
#' @importFrom stats rnorm
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
