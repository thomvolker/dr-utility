
###############################################################################
## functions.R containing code to run simulation 1
###############################################################################


# R-script accompanying the paper "A density ratio framework for evaluating the
# utility of synthetic data"
# 24-04-2024
# Thom Benjamin Volker

###############################################################################
## Functions for generating and evaluating synthetic data
###############################################################################


## Laplace distribution (density function)
dlaplace <- function(x, mu, sigma) {
  exp(-abs(x-mu)/(sigma / sqrt(2))) / (2 * (sigma / sqrt(2)))
}

## Laplace distribution (random number generator)
rlaplace <- function(n, mu, sigma) {
  p <- runif(n, min = -0.5, max = 0.5)
  mu - sigma / sqrt(2) * sign(p) * log(1 - 2 * abs(p))
}

## Ratio between laplace and normal distribution with same mean and variance
dratio_lap_norm <- function(x, mu, sigma) {
  dlaplace(x, mu, sigma) / dnorm(x, mu, sigma)
}

## Ratio between lognormal and normal distribution with same mean and variance
dratio_lnorm_norm <- function(x, mu, sigma) {
  mu_rlnorm    <- log(mu^2 / sqrt(mu^2 + sigma^2))
  sigma_rlnorm <- sqrt(log(1 + sigma^2 / mu^2))
  dlnorm(x, mu_rlnorm, sigma_rlnorm) / dnorm(x, mu, sigma)
}

## Ratio between t-distribution and normal distribution with same mean and variance
dratio_t_norm <- function(x, mu, sigma) {
  df <- 2 / (1 - 1/sigma^2)
  dt(x - mu, df) / dnorm(x - mu, 0, sigma)
}

## Ratio between normal and normal distribution with same mean and variance
## (hence, the ratio equals 1 of course, but I'll use the function for
## consistency nevertheless).
dratio_norm_norm <- function(x, mu, sigma) {
  dnorm(x, mu, sigma) / dnorm(x, mu, sigma)
}

## k-nearest neighbor density ratio estimation

knn_ratio <- function(nu, de, k) {
  nu <- as.matrix(nu)
  de <- as.matrix(de)
  nnu <- nrow(nu)
  nde <- nrow(de)
  p <- ncol(nu)
  Dnu <- densityratio::distance(nu, nu) |> sqrt()
  Dde <- densityratio::distance(nu, de) |> sqrt()
  nn_nu <- apply(Dnu, 1, \(x) order(x)[k+1])
  nn_de <- apply(Dde, 1, \(x) order(x)[k])
  Dnn_nu <- t(Dnu)[0:(nnu-1)*ncol(Dnu) + nn_nu]
  Dnn_de <- t(Dde)[0:(nnu-1)*ncol(Dde) + nn_de]

  nde/nnu * Dnn_de^p / Dnn_nu^p
}


