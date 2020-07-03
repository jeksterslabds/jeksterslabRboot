#' ---
#' title: "Test: jackknife"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: jackknife}
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
context("Test jackknife.")

data <- c(
  11,
  22,
  33
)
thetahat <- mean(data)
jack_samples <- jack(data = data)
thetahat_star <- sapply(
  X = jack_samples,
  FUN = mean
)
results <- jack_hat(
  thetahat_star = thetahat_star,
  thetahat = thetahat
)

#+ testthat_01
test_that("mean", {
  expect_equivalent(
    results["mean"],
    22
  )
})

#+ testthat_02
test_that("bias", {
  expect_equivalent(
    results["bias"],
    0
  )
})

#+ testthat_03
test_that("se", {
  expect_equivalent(
    results["se"],
    0
  )
})

#+ testthat_04
test_that("thetahat_jack", {
  expect_equivalent(
    results["thetahat_jack"],
    22
  )
})
