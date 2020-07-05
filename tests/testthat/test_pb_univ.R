#' ---
#' title: "Test: Parametric Bootstrap - Univariate"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Parametric Bootstrap - Univariate}
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
context("Test Parametric Bootstrap - Univariate.")
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
Variable <- c(
  "`theta`",
  "`sigma2`",
  "`sigma2`"
)
Description <- c(
  "Population mean.",
  "Population variance.",
  "Population standard deviation."
)
Notation <- c(
  "$\\theta = \\mu$",
  "$\\sigma^2$",
  "$\\sigma$"
)
Value <- c(
  theta,
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
#'
#' ## Estimate Correlation
#'
#+ estimate
muhat <- mean(x)
thetahat <- muhat
sigma2hat <- var(x)
sigmahat <- sqrt(sigma2hat)
varhat_thetahat <- sigma2hat / n
sehat_thetahat <- sqrt(varhat_thetahat)
Variable <- c(
  "`n`",
  "`thetahat`",
  "`sigma2hat`",
  "`sigmahat`",
  "`varhat_thetahat`",
  "`sehat_thetahat`"
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
  "$\\hat{\\mathrm{Var}} \\left( \\hat{\\theta} \\right) = \\frac{ \\hat{\\sigma}^2 }{n}$",
  "$\\hat{\\mathrm{se}} \\left( \\hat{\\theta} \\right) = \\frac{ \\hat{\\sigma} }{\\sqrt{n}}$"
)
Value <- c(
  n,
  thetahat,
  sigma2hat,
  sigmahat,
  varhat_thetahat,
  sehat_thetahat
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
x_star <- pb_univ(
  rFUN = rnorm,
  n = n,
  B = B,
  mean = thetahat,
  sd = sigmahat
)
thetahat_star <- sapply(
  X = x_star,
  FUN = mean
)
mean_thetahat_star <- mean(thetahat_star)
var_thetahat_star <- var(thetahat_star)
sd_thetahat_star <- sqrt(var_thetahat_star)
Variable <- c(
  "`B`",
  "`mean_thetahat_star`",
  "`var_thetahat_star`",
  "`sd_thetahat_star`"
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
  "$\\hat{\\mathrm{Var}}_{\\mathrm{B}} \\left( \\hat{\\theta} \\right) = \\frac{1}{B - 1} \\sum_{b = 1}^{B} \\left[ \\hat{\\theta}^{*} \\left( b \\right) - \\hat{\\theta}^{*} \\left( \\cdot \\right) \\right]^2$",
  "$\\hat{\\mathrm{se}}_{\\mathrm{B}} \\left( \\hat{\\theta} \\right) = \\sqrt{ \\frac{1}{B - 1} \\sum_{b = 1}^{B} \\left[ \\hat{\\theta}^{*} \\left( b \\right) - \\hat{\\theta}^{*} \\left( \\cdot \\right) \\right]^2 }$"
)
Value <- c(
  B,
  mean_thetahat_star,
  var_thetahat_star,
  sd_thetahat_star
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
  thetahat_star,
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
#'
#'
#' ### Confidence Intervals
#'
#+ wald
wald_out <- wald(
  thetahat = thetahat,
  sehat_thetahat = sehat_thetahat,
  eval = TRUE
)
wald_out_t <- wald(
  thetahat = thetahat,
  sehat_thetahat = sehat_thetahat,
  distribution = "t",
  df = df,
  eval = TRUE
)
#'
#+ pc
pc_out <- pc(
  thetahat_star = thetahat_star,
  thetahat = thetahat,
  wald = TRUE,
  eval = TRUE
)
#'
#+ bc
bc_out <- bc(
  thetahat_star = thetahat_star,
  thetahat = thetahat,
  wald = TRUE,
  eval = TRUE
)
#'
#+ bca
bca_out <- bca(
  thetahat_star = thetahat_star,
  thetahat = thetahat,
  data = x,
  fitFUN = mean,
  wald = TRUE,
  eval = TRUE
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
    mean_thetahat_star,
    tolerance = 0.01
  )
})
#'
#+ testthat_02
test_that("se", {
  expect_equivalent(
    se_thetahat,
    sd_thetahat_star,
    tolerance = 0.01
  )
})
