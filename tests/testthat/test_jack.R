#' ---
#' title: "Test: jack"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: jack}
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
context("Test jack.")
#'
#' ## Parameters
#'
#+ parameters
data <- c(
  11,
  22,
  33
)
thetahat <- mean(data)
jack_samples <- jack(data = data)
thetahatstarjack <- sapply(
  X = jack_samples,
  FUN = mean
)
#'
#' ## Results
#'
#+ results
results <- jack_hat(
  thetahatstarjack = thetahatstarjack,
  thetahat = thetahat
)
#'
#' ## testthat
#'
#+ testthat_01
test_that("mean", {
  expect_equivalent(
    results[["hat"]][["mean"]],
    22
  )
})
#'
#+ testthat_02
test_that("bias", {
  expect_equivalent(
    results[["hat"]][["bias"]],
    0
  )
})
#'
#+ testthat_03
test_that("se", {
  expect_equivalent(
    results[["hat"]][["se"]],
    0
  )
})
#'
#+ testthat_04
test_that("thetahatjack", {
  expect_equivalent(
    results[["hat"]][["thetahatjack"]],
    22
  )
})

#' ## Example
#+ example
data <- c(
  17.23,
  13.93,
  15.78,
  14.91,
  18.21,
  14.28,
  18.83,
  13.45,
  18.71,
  18.81,
  11.29,
  13.39,
  11.57,
  10.94,
  15.52,
  15.25
)
ps <- c(
  1.605,
  1.151,
  0.998,
  0.942,
  2.416,
  1.043,
  3.122,
  1.362,
  2.972,
  3.097,
  3.308,
  1.393,
  2.951,
  3.806,
  0.958,
  0.937
)
mean_ps <- 2.00389
ci_ll <- 1.45
ci_ul <- 2.56
example <- jack_hat(
  thetahatstarjack = sapply(X = jack(data), FUN = function(x) log(var(x))),
  thetahat = log(var(data))
)
#'
#+ testthat_05
test_that("mean_ps", {
  expect_equivalent(
    round(
      x = mean_ps,
      digits = 3
    ),
    round(
      x = example[["hat"]][["mean_ps"]],
      digits = 3
    )
  )
})
#'
#+ testthat_06
test_that("ps", {
  expect_equivalent(
    round(
      x = ps,
      digits = 3
    ),
    round(
      x = example[["ps"]],
      digits = 3
    )
  )
})
#'
#+ testthat_07
test_that("ll", {
  expect_equivalent(
    ci_ll,
    round(
      x = example[["ci"]][["ci_2.5"]],
      digits = 2
    )
  )
})
#'
#+ testthat_08
test_that("ul", {
  expect_equivalent(
    ci_ul,
    round(
      x = example[["ci"]][["ci_97.5"]],
      digits = 2
    )
  )
})
