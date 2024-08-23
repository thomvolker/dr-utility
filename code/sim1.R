
###############################################################################
## sim1.R containing code to run simulation 1
###############################################################################

# R-script accompanying the paper "A density ratio framework for evaluating the
# utility of synthetic data"
# 24-04-2024
# Thom Benjamin Volker

###############################################################################
## Preliminaries: packages, seed and parallel processing
###############################################################################

# load required packages
library(purrr)
library(furrr)
library(densityratio)
library(dplyr)
library(synthpop)
library(tidyr)

source("code/functions.R")

# set seed for reproducibility
set.seed(123)

# set up parallel processing
plan(multisession, workers = 16)

###############################################################################
## Simulations
###############################################################################

# parameters
nsim  <- 1000
nobs  <- nsyn <- 250
mu    <- 1
sigma <- sqrt(2)

# Simulate data
sims <- expand_grid(model = c("laplace", "lognormal", "t", "normal"),
                    sim = 1:nsim)

sim_dat <- function(model, n, mu, sigma) {
  switch(model,
         laplace = rlaplace(n, mu, sigma),
         lognormal = rlnorm(n,
                            log(mu^2 / sqrt(mu^2 + sigma^2)),
                            sqrt(log(1 + sigma^2 / mu^2))),
         t = rt(n, 2 / (1 - 1/sigma^2)) + mu,
         normal = rnorm(n, mu, sigma))
}

# Create evaluation data for plotting
x_eval <- seq(-3, 5, length.out = 1600)

# Run simulations
sims <- sims |>
  mutate(obs = map(model, ~sim_dat(.x, nobs, mu, sigma)),
         syn = map(obs, ~rnorm(nsyn, mean = mean(.x), sd = sd(.x))),
         r = future_map2(obs, syn, ~densityratio::ulsif(.x, .y, progressbar = FALSE),
                         .options = furrr_options(seed = TRUE),
                         .progress = TRUE),
         summary = future_map(r, ~summary(.x, test = FALSE, parallel = FALSE),
                              .options = furrr_options(seed = TRUE),
                              .progress = TRUE),
         model = factor(model,
                        levels = c("laplace", "lognormal", "t", "normal"),
                        labels = c("Laplace", "Log-normal", "italic(t)", "Normal"),
                        ordered = TRUE),
         xpreds = list(x_eval),
         ypreds = map2(r, xpreds, ~predict(.x, .y)))

plan(sequential)

# True density ratios
true_r <- function(model, x, mu, sigma) {
  switch(model,
         laplace = dratio_lap_norm(x, mu, sigma),
         lognormal = dratio_lnorm_norm(x, mu, sigma),
         t = dratio_t_norm(x, mu, sigma),
         normal = dratio_norm_norm(x, mu, sigma)
         )
}

# True density ratio data for plotting
true_ratio <- data.frame(model = c("laplace", "lognormal", "t", "normal")) |>
  mutate(xpreds = list(x_eval),
         ypreds = map2(model, xpreds, ~true_r(.x, .y, mu, sigma)),
         model  = factor(model,
                         levels = c("laplace", "lognormal", "t", "normal"),
                         labels = c("Laplace", "Log-normal", "italic(t)", "Normal"),
                         ordered = TRUE))

sims <- sims |>
  mutate(rknn = map2(obs, syn, ~knn_ratio(.x, .y, floor(sqrt(length(.x))))))

save(list = c("sims", "true_ratio"),
        file = "results/sim1.RData")
