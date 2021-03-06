#' ---
#' title: "Test: nhstplot"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: nhstplot}
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
context("Test nhstplot.")
#'
#' ## Alpha Level
#'
#' ### Two-Tailed
#'
nhstplot(
  alpha = 0.05,
  dist = "z",
  statistic = 3
)
nhstplot(
  alpha = 0.05,
  dist = "t",
  df = 5
)
#'
#' ### One-Tailed
#'
nhstplot(
  alpha = 0.05,
  two.tailed = FALSE,
  right.tail = TRUE
)
nhstplot(
  alpha = 0.05,
  two.tailed = FALSE,
  right.tail = FALSE
)
nhstplot(
  alpha = 0.05,
  dist = "t",
  two.tailed = FALSE,
  right.tail = TRUE,
  df = 5
)
nhstplot(
  alpha = 0.05,
  dist = "t",
  two.tailed = FALSE,
  right.tail = FALSE,
  df = 5
)
nhstplot(
  alpha = 0.05,
  dist = "F",
  df1 = 10,
  df2 = 3
)
nhstplot(
  alpha = 0.05,
  dist = "chisq",
  df = 3
)
#'
#' ## Confidence Interval
#'
#' ### Two-Tailed
#'
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  dist = "z",
  statistic = 3
)
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  dist = "t",
  df = 5
)
#'
#' ### One-Tailed
#'
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  two.tailed = FALSE,
  right.tail = TRUE
)
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  two.tailed = FALSE,
  right.tail = FALSE
)
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  dist = "t",
  two.tailed = FALSE,
  right.tail = TRUE,
  df = 5
)
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  dist = "t",
  two.tailed = FALSE,
  right.tail = FALSE,
  df = 5
)
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  dist = "F",
  df1 = 10,
  df2 = 3
)
nhstplot(
  alpha = 0.05,
  ci = TRUE,
  dist = "chisq",
  df = 3
)

#+ testthat_01
test_that("error", {
  expect_error(
    nhstplot(
      alpha = c(
        0.001,
        0.01,
        0.05
      )
    )
  )
})
