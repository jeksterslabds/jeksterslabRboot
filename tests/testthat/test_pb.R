#' ---
#' title: "Test: Parametric Bootstrap"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Parametric Bootstrap}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ knitr_options, include=FALSE, cache=FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
#'
#+ setup
library(testthat)
library(MASS)
library(jeksterslabRboot)
context("Test Parametric Bootstrap.")
#'
#' ## Parameters
#'
#+ parameters
n <- 1000
theta <- runif(
  n = 1,
  min = -1,
  max = 1
)
Sigma <- matrix(
  data = c(
    1,
    theta,
    theta,
    1
  ),
  nrow = 2
)
mu <- c(0, 0)
var_thetahat <- ((1 - theta^2)^2) / n
Variable <- c(
  "`n`",
  "`theta`",
  "`var_thetahat`"
)
Description <- c(
  "Sample size ($n$).",
  "Population correlation ($\\theta$).",
  "Variance of the sampling distribution of sample correlation ($\\mathrm{Var} \\left( \\hat{\\theta} \\right) = \\frac{\\left( 1 - \\theta^2 \\right)^2}{n}$)."
)
Value <- c(
  n,
  theta,
  var_thetahat
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Value
  ),
  row.names = FALSE,
  caption = "Population Parameters"
)
#'
#' ## Generate Data
#'
#+ generate_data
X <- mvrnorm(
  n = n,
  Sigma = Sigma,
  mu = mu
)
#'
#' ## Estimate Correlation
#'
#+ estimate
thetahat <- cor(X)[2, 1]
varhat_thetahat <- ((1 - thetahat^2)^2) / (n - 2)
ci <- cor.test(X[, 1], X[, 2])
ci_ll <- ci$conf.int[1]
ci_ul <- ci$conf.int[2]
Variable <- c(
  "`n`",
  "`thetahat`",
  "`varhat_thetahat`"
)
Description <- c(
  "Sample size ($n$).",
  "Sample correlation ($\\hat{\\theta}$).",
  "Estimate of the variance of the sampling distribution of sample correlation ($\\hat{\\mathrm{Var}} \\left( \\hat{\\theta} \\right) = \\frac{\\left( 1 - \\hat{\\theta}^2 \\right)^2}{n}$)."
)
Value <- c(
  n,
  thetahat,
  varhat_thetahat
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Value
  ),
  row.names = FALSE,
  caption = "Sample Statistics"
)
#'
#' ## Bootstrap
#'
#+ bootstrap
B <- 10000
X_star <- pb_mvn(
  data = X,
  B = B
)
thetahat_star <- lapply(
  X = X_star,
  FUN = function(X) cor(X)[2, 1]
)
thetahat_star <- as.vector(
  do.call(
    what = "rbind",
    args = thetahat_star
  )
)
mean_thetahat_star <- mean(thetahat_star)
var_thetahat_star <- var(thetahat_star)
Variable <- c(
  "`B`",
  "`mean_thetahat_star`",
  "`var_thetahat_star`"
)
Description <- c(
  "Bootstrap samples($B$).",
  "Mean of $B$ sample correlations.",
  "Variance of $B$ sample correlations."
)
Value <- c(
  B,
  mean_thetahat_star,
  var_thetahat_star
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Value
  ),
  row.names = FALSE,
  caption = "Parametric Bootstrapping Results"
)
#'
#' ### Bootstrap Sampling Distribution
#'
#+ plot
hist(
  thetahat_star,
  breaks = 100
)
#'
#' ### Confidence Intervals
#'
#+ wald
wald <- wald(
  thetahat = thetahat,
  sehat_thetahat = sqrt(varhat_thetahat),
  theta_null = 0,
  alpha = c(
    0.001,
    0.01,
    0.05
  ),
  distribution = "z"
)
#'
#+ pc
pc <- pc(
  thetahat_star = thetahat_star,
  alpha = c(
    0.001,
    0.01,
    0.05
  ),
  wald = TRUE,
  thetahat = thetahat,
  theta_null = 0,
  distribution = "z"
)
#'
#+ bc
bc <- bc(
  thetahat_star = thetahat_star,
  thetahat = thetahat,
  alpha = c(
    0.001,
    0.01,
    0.05
  ),
  wald = TRUE,
  theta_null = 0,
  distribution = "z"
)
#'
#+ bca
bca <- bca(
  thetahat_star = thetahat_star,
  thetahat = thetahat,
  data = X,
  fitFUN = function(X) cor(X)[2, 1],
  alpha = c(
    0.001,
    0.01,
    0.05
  ),
  wald = TRUE,
  theta_null = 0,
  distribution = "z"
)
#'
#+ ci
knitr::kable(
  x = as.data.frame(
    rbind(
      wald,
      pc,
      bc,
      bca
    )
  ),
  caption = "Confidence Intervals"
)
#'
#' ## testthat
#'
#+ testthat_01, echo = TRUE
test_that("pb mean_thetahat_star is equal to thetahat", {
  expect_equivalent(
    round(
      x = thetahat,
      digits = 1
    ),
    round(
      x = mean_thetahat_star,
      digits = 1
    )
  )
})
#'
#+ testthat_02, echo = TRUE
test_that("pb var_thetahat_star is equal to varhat_thetahat", {
  expect_equivalent(
    round(
      x = varhat_thetahat,
      digits = 1
    ),
    round(
      x = var_thetahat_star,
      digits = 1
    )
  )
})
#'
#+ testthat_03, echo = TRUE
test_that("wald alpha 0.05 is equal to cor.test", {
  expect_equivalent(
    round(
      x = ci_ll,
      digits = 1
    ),
    round(
      x = wald["ci_2.5"],
      digits = 1
    )
  )
  expect_equivalent(
    round(
      x = ci_ul,
      digits = 1
    ),
    round(
      x = wald["ci_97.5"],
      digits = 1
    )
  )
})

#+ univ
# Normal distribution
n <- 1000
B <- 5000
mu <- 100
sigma2 <- 225
sigma <- sqrt(sigma2)
data <- rnorm(
  n = n,
  mean = mu,
  sd = sigma
)
muhat <- mean(data)
sigmahat <- sd(data)
x_star_normal <- pb_univ(
  rFUN = rnorm,
  n = n,
  B = B,
  mean = muhat,
  sd = sigmahat
)
thetahat_star_normal <- sapply(
  X = x_star_normal,
  FUN = mean
)
# Binomial distribution
n_trials <- 1
p <- 0.50
data <- rbinom(
  n = n,
  size = n_trials,
  prob = p
)
phat <- mean(data) / n_trials
x_star_normal_binomial <- pb_univ(
  rFUN = rbinom,
  n = n,
  B = B,
  size = n_trials,
  prob = phat
)
thetahat_star_binomial <- sapply(
  X = x_star_normal_binomial,
  FUN = mean
)
# test_that("normal", {
#  expect_equivalent(
#    round(
#      x = mu,
#      digits = 0
#    ),
#    round(
#      x = mean(thetahat_star_normal),
#      digits = 0
#    )
#  )
# })
# test_that("binomial", {
#  expect_equivalent(
#    round(
#      x = p,
#      digits = 1
#    ),
#    round(
#      x = mean(thetahat_star_binomial) / n_trials,
#      digits = 1
#    )
#  )
# })

#+ for_coverage
 R <- 10L
 # length 1
 thetahat <- 100
 covhat_thetahat <- 1.5
 mc_star_length_1 <- mc(
   thetahat = thetahat,
   covhat_thetahat = covhat_thetahat,
   R = R
 )
 head(mc_star_length_1)
 # length greater than 1
 alphahat <- 0.3386
 betahat <- 0.4510
 alphahat_betahat <- alphahat * betahat
 varhat_alphahat <- 0.1224^2
 varhat_betahat <- 0.1460^2
 thetahat <- c(
   alphahat,
   betahat
 )
 covhat_thetahat <- matrix(
   data = c(
     varhat_alphahat,
     0.00,
     0.00,
     varhat_betahat
   ),
   ncol = 2
 )
 mc_star_length_2 <- mc(
   thetahat = thetahat,
   covhat_thetahat = covhat_thetahat,
   R = R
 )
