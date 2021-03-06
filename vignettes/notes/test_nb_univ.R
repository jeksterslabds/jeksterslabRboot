#' ---
#' title: "Test: Nonparametric Bootstrap - Univariate"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Nonparametric Bootstrap - Univariate}
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
library(jeksterslabRboot)
context("Test Nonparametric Bootstrap - Univariate.")
#'
#' ## Parameters
#'
#+ parameters
n <- 1000
df <- n - 1
B <- 1000
mu <- runif(
  n = 1,
  min = .5,
  max = 1
)
theta <- mu
sigma2 <- runif(
  n = 1,
  min = .5,
  max = 1
)
sigma <- sqrt(sigma2)
var_thetahat <- sigma2 / n
se_thetahat <- sqrt(var_thetahat)
# Parameters
Variable <- c(
  "`mu`",
  "`sigma2`",
  "`sigma2`"
)
Description <- c(
  "Population mean.",
  "Population variance.",
  "Population standard deviation."
)
Notation <- c(
  "$\\mu$",
  "$\\sigma^2$",
  "$\\sigma$"
)
Value <- c(
  mu,
  sigma2,
  sigma
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  row.names = FALSE,
  caption = "Population Parameters"
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
  "Population mean.",
  "Variance of the sampling distribution of the mean.",
  "Standard error of the mean."
)
Notation <- c(
  "$n$",
  "$\\theta = \\mu$",
  "$\\mathrm{Var} \\left( \\hat{\\theta} \\right) = \\frac{ \\sigma^2 }{n}$",
  "$\\mathrm{se} \\left( \\hat{\\theta} \\right) = \\frac{ \\sigma }{\\sqrt{n}}$"
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
  caption = "Sampling Distribution of $\\hat{\\theta}$ with Known Parameters"
)
#'
#' ## Generate Data
#'
#+ generate_data
x <- rnorm(
  n = n,
  mean = mu,
  sd = sigma
)
str(x)
hist(x)
qqnorm(x)
qqline(x)
#'
#' ## Estimate Correlation
#'
#+ estimate
muhat <- mean(x)
thetahat <- muhat
sigma2hat <- var(x)
sigmahat <- sqrt(sigma2hat)
varhat_thetahat <- sigma2hat / n
sehat <- sqrt(varhat_thetahat)
Variable <- c(
  "`n`",
  "`thetahat`",
  "`sigma2hat`",
  "`sigmahat`",
  "`varhat_thetahat`",
  "`sehat`"
)
Description <- c(
  "Sample size.",
  "Sample mean.",
  "Sample variance.",
  "Sample standard deviation.",
  "Estimate of the variance of the sampling distribution of the mean.",
  "Estimate of the standard error of the mean."
)
Notation <- c(
  "$n$",
  "$\\hat{\\theta} = \\hat{\\mu} = \\frac{1}{n} \\sum_{i = 1}^{n} x_i$",
  "$\\hat{\\sigma}^2 = \\frac{1}{n - 1} \\sum_{i = 1}^{n} \\left( x_i - \\hat{\\mu} \\right)^2$",
  "$\\hat{\\sigma} = \\sqrt{\\frac{1}{n - 1} \\sum_{i = 1}^{n} \\left( x_i - \\hat{\\mu} \\right)^2}$",
  "$\\widehat{\\mathrm{Var}} \\left( \\hat{\\theta} \\right) = \\frac{ \\hat{\\sigma}^2 }{n}$",
  "$\\widehat{\\mathrm{se}} \\left( \\hat{\\theta} \\right) = \\frac{ \\hat{\\sigma} }{\\sqrt{n}}$"
)
Value <- c(
  n,
  thetahat,
  sigma2hat,
  sigmahat,
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
x_star <- nb(
  data = x,
  B = B
)
thetahatstar <- sapply(
  X = x_star,
  FUN = mean
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
  "Number of bootstrap samples.",
  "Mean of $B$ sample means.",
  "Variance of $B$ sample means.",
  "Standard deviation of $B$ sample means."
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
wald_out_t <- wald(
  thetahat = thetahat,
  sehat = sehat,
  dist = "t",
  df = df,
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
# for coverage
pc_out_t <- pc(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  dist = "t",
  df = df
)
#'
#+ bc
bc_out <- bc(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  wald = TRUE,
  eval = TRUE
)
# for coverage
bc_out_t <- bc(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  dist = "t",
  df = df
)
#'
#+ bca
bca_out <- bca(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  data = x,
  fitFUN = mean,
  wald = TRUE,
  eval = TRUE
)
# for coverage
bca_out_t <- bca(
  thetahatstar = thetahatstar,
  thetahat = thetahat,
  data = x,
  fitFUN = mean,
  dist = "t",
  df = df
)
#'
#+ ci
knitr::kable(
  x = as.data.frame(
    rbind(
      wald = wald_out,
      wald_t = wald_out_t,
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
