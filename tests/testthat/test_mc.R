#' ---
#' title: "Test: mc"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: mc}
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
context("Test mc.")
#'
#' ## Length 1
#'
#' ### Parameters
#'
#+ parameters_01
R <- 20000L
thetahat <- 100
vcovhat <- 1.5
#'
#' ## Results
#'
#+ length_01
mc_star_length_1 <- mc(
  thetahat = thetahat,
  vcovhat = vcovhat,
  R = R
)
str(mc_star_length_1)
hist(
  mc_star_length_1,
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
qqnorm(mc_star_length_1)
qqline(mc_star_length_1)
#'
#' ## Length Greater than 1
#'
#' ### Parameters
#'
#+ parameters_02
alphahat <- 0.3386
betahat <- 0.4510
alphahat_betahat <- alphahat * betahat
alphahat_betahat_ci_2.5 <- 0.0033
alphahat_betahat_ci_97.5 <- 0.2979
varhat_alphahat <- 0.1224^2
varhat_betahat <- 0.1460^2
thetahat <- c(
  alphahat,
  betahat
)
vcovhat <- matrix(
  data = c(
    varhat_alphahat,
    0.00,
    0.00,
    varhat_betahat
  ),
  ncol = 2
)
#'
#' ## Results
#'
#+ mediation
mc_star_length_2 <- mc(
  thetahat = thetahat,
  vcovhat = vcovhat,
  R = R
)
str(mc_star_length_2)
alphahat_betahat_star <- mc_star_length_2[, 1] * mc_star_length_2[, 2]
hist(
  alphahat_betahat_star,
  main = expression(
    paste(
      "Histogram of ",
      hat(alpha),
      hat(beta),
      "*"
    )
  ),
  xlab = expression(
    paste(
      hat(alpha),
      hat(beta),
      "*"
    )
  )
)
qqnorm(alphahat_betahat_star)
qqline(alphahat_betahat_star)
wald_out <- wald(
  thetahat = alphahat_betahat,
  sehat = sd(alphahat_betahat_star)
)
wald_out
#'
#' ## testthat
#'
#+ testthat_01
test_that("ci_2.5", {
  expect_equivalent(
    alphahat_betahat_ci_2.5,
    wald_out["ci_2.5"],
    tolerance = 0.05
  )
})
#'
#+ testthat_02
test_that("ci_97.5", {
  expect_equivalent(
    alphahat_betahat_ci_97.5,
    wald_out["ci_97.5"],
    tolerance = 0.05
  )
})
