#' ---
#' title: "Test: Nonparametric Bootstrap - Bivariate"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Nonparametric Bootstrap - Bivariate}
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
context("Test Nonparametric Bootstrap - Bivariate.")
#'
#' ## Parameters
#'
#+ parameters
n <- 1000
B <- 1000
rho <- runif(
  n = 1,
  min = .5,
  max = .8
)
theta <- rho
Sigma <- matrix(
  data = c(
    1,
    rho,
    rho,
    1
  ),
  nrow = 2
)
mu <- c(0, 0)
var_thetahat <- ((1 - rho^2)^2) / n
se_thetahat <- sqrt(var_thetahat)
# Parameters
Variable <- c(
  "`rho`"
)
Description <- c(
  "Population correlation."
)
Notation <- c(
  "$\\rho$"
)
Value <- c(
  theta
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  row.names = FALSE,
  caption = "Population Parameter"
)
# Sampling distribution of theta
Variable <- c(
  "`n`",
  "`theta`",
  "`var_thetahat`",
  "`se_thetahat`"
)
Description <- c(
  "Sample size.",
  "Population correlation.",
  "Variance of the sampling distribution of the sample correlation.",
  "Standard error of the sampling distribution of the sample correlation."
)
Notation <- c(
  "$n$",
  "$\\theta = \\rho$",
  "$\\mathrm{Var} \\left( \\hat{\\theta} \\right) = \\frac{ \\left( 1 - \\rho^2 \\right)^2}{n}$",
  "$\\mathrm{se} \\left( \\hat{\\theta} \\right) = \\frac{ 1 - \\rho^2 }{\\sqrt{n}}$"
)
Value <- c(
  n,
  theta,
  var_thetahat,
  se_thetahat
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  row.names = FALSE,
  caption = "Sampling Distribution of $\\hat{\\theta}$ with Known Parameter"
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
str(X)
hist(X[, 1])
qqnorm(X[, 1])
qqline(X[, 1])
hist(X[, 2])
qqnorm(X[, 2])
qqline(X[, 2])
plot(X[, 1], X[, 2])
#'
#' ## Estimate Correlation
#'
#+ estimate
rhohat <- cor(X)[2, 1]
thetahat <- rhohat
varhat_thetahat <- ((1 - rhohat^2)^2) / (n - 2)
sehat <- sqrt(varhat_thetahat)
ci <- cor.test(X[, 1], X[, 2])
ci_ll <- ci$conf.int[1]
ci_ul <- ci$conf.int[2]
Variable <- c(
  "`n`",
  "`thetahat`",
  "`varhat_thetahat`",
  "`sehat`"
)
Description <- c(
  "Sample size.",
  "Sample correlation.",
  "Estimate of the variance of the sampling distribution of the sample correlation.",
  "Estimate of the standard error of the sampling distribution of the sample correlation."
)
Notation <- c(
  "$n$",
  "$\\hat{\\theta} = \\hat{\\rho}_{x, y} = \\frac{\\sum_{i = 1}^{n} \\left(x_i - \\hat{\\mu}_{x} \\right) \\left(y_i - \\hat{\\mu}_{y} \\right)}{\\sqrt{\\sum_{i = 1}^{n} \\left( x_i -  \\hat{\\mu}_{x} \\right)^2} \\sqrt{\\sum_{i = 1}^{n} \\left( y_i -  \\hat{\\mu}_{y} \\right)^2}}$",
  "$\\widehat{\\mathrm{Var}} \\left( \\hat{\\theta} \\right) = \\frac{ \\left( 1 - \\hat{\\rho}^2 \\right)^2}{n - 2}$",
  "$\\widehat{\\mathrm{se}} \\left( \\hat{\\theta} \\right) = \\frac{ 1 - \\hat{\\rho}^2 }{\\sqrt{n - 2}}$"
)
Value <- c(
  n,
  thetahat,
  varhat_thetahat,
  sehat
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  row.names = FALSE,
  caption = "Sample Statistics (Parameter Estimates)"
)
#'
#' ## Bootstrap
#'
#+ bootstrap
X_star <- nb(
  data = X,
  B = B
)
thetahatstar <- sapply(
  X = X_star,
  FUN = function(X) cor(X)[2, 1]
)
mean_thetahatstar <- mean(thetahatstar)
var_thetahatstar <- var(thetahatstar)
sd_thetahatstar <- sqrt(var_thetahatstar)
Variable <- c(
  "`B`",
  "`mean_thetahatstar`",
  "`var_thetahatstar`",
  "`sd_thetahatstar`"
)
Description <- c(
  "Bootstrap samples.",
  "Mean of $B$ sample correlations.",
  "Variance of $B$ sample correlations.",
  "Standard deviation of $B$ sample correlations."
)
Notation <- c(
  "$B$",
  "$\\hat{\\theta}^{*} \\left( \\cdot \\right) = \\frac{1}{B} \\sum_{b = 1}^{B} \\hat{\\theta}^{*} \\left( b \\right)$",
  "$\\widehat{\\mathrm{Var}}_{\\mathrm{B}} \\left( \\hat{\\theta} \\right) = \\frac{1}{B - 1} \\sum_{b = 1}^{B} \\left[ \\hat{\\theta}^{*} \\left( b \\right) - \\hat{\\theta}^{*} \\left( \\cdot \\right) \\right]^2$",
  "$\\widehat{\\mathrm{se}}_{\\mathrm{B}} \\left( \\hat{\\theta} \\right) = \\sqrt{ \\frac{1}{B - 1} \\sum_{b = 1}^{B} \\left[ \\hat{\\theta}^{*} \\left( b \\right) - \\hat{\\theta}^{*} \\left( \\cdot \\right) \\right]^2 }$"
)
Value <- c(
  B,
  mean_thetahatstar,
  var_thetahatstar
)
Value <- c(
  B,
  mean_thetahatstar,
  var_thetahatstar,
  sd_thetahatstar
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  row.names = FALSE,
  caption = "Nonparametric Bootstrapping Results"
)
#'
#' ### Bootstrap Sampling Distribution
#'
#+ plot
hist(
  thetahatstar,
  breaks = 100,
  main = expression(
    paste(
      "Histogram of ",
      hat(theta),
      "*"
    )
  ),
  xlab = expression(
    paste(
      hat(theta),
      "*"
    )
  )
)
qqnorm(thetahatstar)
qqline(thetahatstar)
#'
#'
#' ### Confidence Intervals
#'
#+ wald
wald_out <- wald(
  thetahat = thetahat,
  sehat = sehat,
  eval = TRUE
)
#'
#+ pc
pc_out <- pc(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  wald = TRUE,
  eval = TRUE
)
#'
#+ bc
bc_out <- bc(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  wald = TRUE,
  eval = TRUE
)
#'
#+ bca
bca_out <- bca(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  data = X,
  fitFUN = function(X) cor(X)[2, 1],
  wald = TRUE,
  eval = TRUE
)
#'
#+ ci
knitr::kable(
  x = as.data.frame(
    rbind(
      wald = wald_out,
      pc = pc_out,
      bc = bc_out,
      bca = bca_out
    )
  ),
  caption = "Confidence Intervals"
)
#'
#' ## testthat
#'
#+ testthat_01
test_that("mean", {
  expect_equivalent(
    thetahat,
    mean_thetahatstar,
    tolerance = 0.01
  )
})
#'
#+ testthat_02
test_that("se", {
  expect_equivalent(
    se_thetahat,
    sd_thetahatstar,
    tolerance = 0.01
  )
})
#'
#+ testthat_03
test_that("2.5 is equal to cor.test", {
  expect_equivalent(
    round(
      x = ci_ll,
      digits = 2
    ),
    round(
      x = wald_out["ci_2.5"],
      digits = 2
    ),
    tolerance = 0.02
  )
})
#'
#+ testthat_04
test_that("97.5 is equal to cor.test", {
  expect_equivalent(
    round(
      x = ci_ul,
      digits = 2
    ),
    round(
      x = wald_out["ci_97.5"],
      digits = 2
    ),
    tolerance = 0.02
  )
})
#'
#+ testthat_05
test_that("2.5", {
  expect_equivalent(
    round(
      x = wald_out["ci_2.5"],
      digits = 2
    ),
    round(
      x = pc_out["ci_2.5"],
      digits = 2
    ),
    round(
      x = bc_out["ci_2.5"],
      digits = 2
    ),
    round(
      x = bca_out["ci_2.5"],
      digits = 2
    ),
    tolerance = 0.02
  )
})
#'
#+ testthat_06
test_that("97.5", {
  expect_equivalent(
    round(
      x = wald_out["ci_97.5"],
      digits = 2
    ),
    round(
      x = pc_out["ci_97.5"],
      digits = 2
    ),
    round(
      x = bc_out["ci_97.5"],
      digits = 2
    ),
    round(
      x = bca_out["ci_97.5"],
      digits = 2
    ),
    tolerance = 0.02
  )
})
