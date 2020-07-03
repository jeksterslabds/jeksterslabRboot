#' ---
#' title: "Test: alpha2prob"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: alpha2prob}
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
context("Test alpha2prob.")

#+ testthat_01
test_that("vector", {
  expect_equivalent(
    alpha2prob(
      alpha = c(
        0.001,
        0.01,
        0.05
      )
    ),
    alpha2prob(
      alpha = c(
        0.05,
        0.01,
        0.001
      )
    ),
    ci2prob(
      ci = c(
        0.999,
        0.99,
        0.95
      )
    ),
    ci2prob(
      ci = c(
        0.95,
        0.99,
        0.999
      )
    ),
    c(
      0.0005,
      0.0050,
      0.0250,
      0.9750,
      0.9950,
      0.9995
    )
  )
})

#+ testthat_02
test_that("vector", {
  expect_equivalent(
    alpha2prob(
      alpha = 0.05
    ),
    ci2prob(
      ci = 0.95
    ),
    c(
      0.025,
      0.975
    )
  )
})
