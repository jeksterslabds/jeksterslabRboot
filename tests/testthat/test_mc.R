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
context("Test mc.")

#+ coverage
R <- 20000L
# length 1
thetahat <- 100
covhat_thetahat <- 1.5
mc_star_length_1 <- mc(
  thetahat = thetahat,
  covhat_thetahat = covhat_thetahat,
  R = R
)
head(mc_star_length_1)
# length greater than 1
alphahat <- 0.3386
betahat <- 0.4510
alphahat_betahat <- alphahat * betahat
varhat_alphahat <- 0.1224^2
varhat_betahat <- 0.1460^2
thetahat <- c(
  alphahat,
  betahat
)
covhat_thetahat <- matrix(
  data = c(
    varhat_alphahat,
    0.00,
    0.00,
    varhat_betahat
  ),
  ncol = 2
)
mc_star_length_2 <- mc(
  thetahat = thetahat,
  covhat_thetahat = covhat_thetahat,
  R = R
)
head(mc_star_length_2)
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
wald(
  thetahat = alphahat_betahat,
  sehat_thetahat = sd(alphahat_betahat_star)
)