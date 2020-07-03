#' ---
#' title: "Test: wald"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: wald}
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
context("Test wald.")

#+ parameters
thetahat <- 0.860
sehat_thetahat <- 0.409
sqrt_wald <- 2.10
p <- 0.036

#+ results
result_01 <- sqrt_wald_test(
  thetahat = thetahat,
  sehat_thetahat = sehat_thetahat,
  theta_null = 0,
  distribution = "z"
)
result_02 <- sqrt_wald_test(
  thetahat = thetahat,
  sehat_thetahat = sehat_thetahat,
  theta_null = 0,
  distribution = "t",
  df = 1000
)
result_03 <- wald_test(
  thetahat = thetahat,
  varhat_thetahat = sehat_thetahat^2
)
wald(
  thetahat = thetahat,
  sehat_thetahat = sehat_thetahat,
  distribution = "z",
  eval = TRUE,
  theta = 0
)
wald(
  thetahat = thetahat,
  sehat_thetahat = sehat_thetahat,
  distribution = "t",
  df = 1000,
  eval = TRUE,
  theta = 0
)

#+ testthat_01
test_that("statistic", {
  expect_equivalent(
    result_01[1],
    result_02[1],
    sqrt(result_03[1]),
    sqrt_wald
  )
})

#+ testthat_02
test_that("p", {
  expect_equivalent(
    round(
      result_01[2],
      digits = 3
    ),
    round(
      result_02[2],
      digits = 3
    ),
    p
  )
})