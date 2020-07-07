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
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams jeksterslabRutils::util_lapply
#' @param data Vector, matrix or data frame.
#'   Sample data to bootstrap.
#'   The empirical distribution \eqn{\hat{F}}.
#' @param B Integer. Number of bootstrap samples.
#' @return Returns a list of length `B` of nonparametric bootstrap samples.
#' @examples
#' B <- 5L
#' n <- 5
#' #################################
#' # vector
#' #################################
#' x <- rnorm(n = n)
#' xstar <- nb(
#'   data = x,
#'   B = B
#' )
#' str(xstar)
#' #################################
#' # matrix
#' #################################
#' x1 <- rnorm(n = n)
#' x2 <- rnorm(n = n)
#' x3 <- rnorm(n = n)
#' X <- cbind(x1, x2, x3)
#' Xstar <- nb(
#'   data = X,
#'   B = B
#' )
#' str(Xstar)
#' #################################
#' # data frame
#' #################################
#' X <- as.data.frame(X)
#' Xstar <- nb(
#'   data = X,
#'   B = B
#' )
#' str(Xstar)
#' @references
#' Efron, B., & Tibshirani, R. J. (1993).
#' *An introduction to the bootstrap*.
#' New York, N.Y: Chapman & Hall.
#'
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
  foo <- function(iter,
                  data) {
    i <- sample(
      iter = len,
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
      iter = 1:B,
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
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams nb
#' @inherit nb references
#' @param n Integer.
#'   Sample size.
#' @param Sigmahat Matrix.
#'   Estimated variance-covariance matrix
#'   from the original sample data.
#' @param muhat Vector.
#'   Estimated mean vector from the original sample data.
#' @return Returns a list of length `B` of parametric bootstrap samples.
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
#' Xstar <- .pb_mvn(
#'   n = 5,
#'   Sigmahat = Sigmahat,
#'   muhat = muhat,
#'   B = B
#' )
#' str(Xstar)
#' @importFrom MASS mvrnorm
#' @export
.pb_mvn <- function(n,
                    Sigmahat,
                    muhat,
                    B = 2000L,
                    par = FALSE,
                    ncores = NULL) {
  foo <- function(iter,
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
      iter = 1:B,
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
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams nb
#' @inherit .pb_mvn references description details return
#' @examples
#' X <- matrix(
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
#' Xstar <- pb_mvn(
#'   data = X,
#'   B = 5L
#' )
#' str(Xstar)
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
#' @family bootstrap functions
#' @keywords bootstrap
#' @inheritParams nb
#' @inherit .pb_mvn references return
#' @param rFUN Function.
#'   Data generating function to generate univariate data.
#' @param n Integer.
#'   Sample size.
#' @param ... Arguments to pass to rFUN.
#' @examples
#' n <- 5
#' B <- 5
#' #################################
#' # normal distribution
#' #################################
#' mu <- 100
#' sigma2 <- 225
#' sigma <- sqrt(sigma2)
#' x <- rnorm(
#'   n = n,
#'   mean = mu,
#'   sd = sigma
#' )
#' muhat <- mean(x)
#' sigmahat <- sd(x)
#' xstar <- pb_univ(
#'   rFUN = rnorm,
#'   n = n,
#'   B = B,
#'   mean = muhat,
#'   sd = sigmahat
#' )
#' str(xstar)
#' #################################
#' # binomial distribution
#' #################################
#' n_trials <- 1
#' p <- 0.50
#' x <- rbinom(
#'   n = n,
#'   size = n_trials,
#'   prob = p
#' )
#' phat <- mean(x) / n_trials
#' xstar <- pb_univ(
#'   rFUN = rbinom,
#'   n = n,
#'   B = B,
#'   size = n_trials,
#'   prob = phat
#' )
#' str(xstar)
#' @importFrom stats rnorm
#' @export
pb_univ <- function(rFUN = rnorm,
                    n,
                    B = 2000L,
                    par = FALSE,
                    ncores = NULL,
                    ...) {
  args <- list(
    iter = 1:B,
    rFUN = rFUN,
    n = n,
    ...
  )
  foo <- function(iter,
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
