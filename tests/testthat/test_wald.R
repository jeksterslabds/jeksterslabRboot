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
#'
#' ## Parameters
#'
#+ parameters
thetahat <- 0.860
sehat <- 0.409
sqrt_wald <- 2.10
p <- 0.036
#'
#' ## Results
#'
#+ results
result_01 <- sqrt_wald_test(
  thetahat = thetahat,
  sehat = sehat,
  null = 0,
  dist = "z"
)
result_01
result_02 <- sqrt_wald_test(
  thetahat = thetahat,
  sehat = sehat,
  null = 0,
  dist = "t",
  df = 1000
)
result_02
result_03 <- wald_test(
  thetahat = thetahat,
  varhat = sehat^2
)
result_03
#'
#' ## testthat
#'
#+ testthat_01
test_that("statistic", {
  expect_equivalent(
    result_01[1],
    result_02[1],
    sqrt(result_03[1]),
    sqrt_wald
  )
})
#'
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
