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
library(MASS)
library(jeksterslabRboot)
context("Test eval.")

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
