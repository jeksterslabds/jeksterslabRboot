#' ---
#' title: "Test: eval"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: eval}
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
context("Test eval.")

#+ results
result <- ci_eval(
  ci = c(1, 3),
  thetahat = 2,
  theta = 2
)

#+ testthat_01
test_that("TRUE", {
  expect_true(
    zero_hit(
      thetahat_lo = -1,
      thetahat_up = 1
    ),
    theta_hit(
      thetahat_lo = -1,
      theta = 0,
      thetahat_up = 1
    )
  )
})

#+ testthat_02
test_that("FALSE", {
  expect_false(
    zero_hit(
      thetahat_lo = 1,
      thetahat_up = 2
    ),
    theta_hit(
      thetahat_lo = 1,
      theta = 0,
      thetahat_up = 2
    )
  )
})

#+ testthat_03
test_that("1", {
  expect_equal(
    len(
      thetahat_lo = 1,
      thetahat_up = 2
    ),
    shape(
      thetahat_lo = -1,
      thetahat = 0,
      thetahat_up = 1
    ),
    1
  )
})

#+ testthat_04
test_that("zero_hit", {
  expect_false(
    as.logical(result[1])
  )
})

#+ testthat_05
test_that("theta_hit", {
  expect_true(
    as.logical(result[2])
  )
})

#+ testthat_06
test_that("len", {
  expect_equivalent(
    result[3],
    2
  )
})

#+ testthat_07
test_that("shape", {
  expect_equivalent(
    result[4],
    1
  )
})
