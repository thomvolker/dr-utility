
## Packages
library(ggplot2)
library(densityratio)
library(synthpop)
library(purrr)
library(tidyr)
library(dplyr)
library(pbapply)


################################################################################
## Data generation functions
################################################################################


gen_real <- function(n, P, rho = 0.5) {

  if (P < 3) stop("P must be at least 4")
  Pnonlin <- ifelse(P > 5, 5, 1)
  Plin    <- P - Pnonlin

  Xlin  <- matrix(rnorm(n*Plin), n) %*% chol(rho + (1-rho) * diag(Plin))
  Xpoly <- sapply(1:Pnonlin, function(i) {
    Xi <- Xlin[,i]^{i+1}
    Xi + rnorm(n, 0, sd(Xi))
  })
  cbind(Xlin, Xpoly) |> unname()
}

gen_syn1 <- function(nsyn, obs) {

  P <- ncol(obs)
  means <- colMeans(obs)
  vars  <- diag(diag(var(obs)))

  matrix(rnorm(nsyn*P), nsyn) %*% chol(vars) + matrix(means, nsyn, P, byrow = TRUE)
}

gen_syn2 <- function(nsyn, obs) {

  P <- ncol(obs)
  means <- colMeans(obs)
  vars <- var(obs)

  matrix(rnorm(nsyn*P), nsyn) %*% chol(vars) + matrix(means, nsyn, P, byrow = TRUE)
}

gen_syn3 <- function(nsyn, obs) {

  P         <- ncol(obs)
  Pnonlin   <- ifelse(P > 5, 5, 1)
  Plin      <- P - Pnonlin
  means_lin <- colMeans(obs[,1:Plin])
  vars_lin  <- var(obs[,1:Plin])

  Xlin <- matrix(rnorm(nsyn*Plin), nsyn) %*% chol(vars_lin) + matrix(means_lin, nsyn, Plin, byrow = TRUE)
  Xpoly <- sapply(1:Pnonlin, function(i) {
    fit <- lm(obs[,Plin + i] ~ I(obs[,i]^{i+1}))
    cbind(1, Xlin[,i]^{i+1}) %*% coef(fit) + rnorm(nsyn, 0, sd(fit$residuals))
  })

  cbind(Xlin, Xpoly) |> unname()
}

################################################################################
## Evaluation functions
################################################################################

pearson <- function(real, syn1, syn2, syn3) {
  std_means <- colMeans(real)
  std_sds  <- apply(real, 2, sd)
  real_std <- scale(real, center = std_means, scale = std_sds)

  R1 <- ulsif(real_std,
              scale(syn1, std_means, std_sds),
              scale = NULL,
              progressbar = FALSE)

  R2 <- ulsif(real_std,
              scale(syn2, std_means, std_sds),
              scale = NULL,
              progressbar = FALSE)

  R3 <- ulsif(real_std,
              scale(syn3, std_means, std_sds),
              scale = NULL,
              progressbar = FALSE)

  sapply(list(R1, R2, R3), \(x) summary(x)$PE) |>
    t() |>
    set_names(c("syn1", "syn2", "syn3")) |>
    as.data.frame()
}

spmse_cart <- function(real, syn1, syn2, syn3) {

  sapply(list(syn1, syn2, syn3), \(syn) {
    synthpop::utility.gen(
      as.data.frame(syn),
      as.data.frame(real),
      print.flag = FALSE
    )$S_pMSE
  }) |>
    t() |>
    set_names(c("syn1", "syn2", "syn3")) |>
    as.data.frame()
}

spmse_logit <- function(real, syn1, syn2, syn3, maxorder = 0) {

  sapply(list(syn1, syn2, syn3), \(syn) {
    synthpop::utility.gen(
      as.data.frame(syn),
      as.data.frame(real),
      print.flag = FALSE,
      method = "logit",
      maxorder = maxorder
    )$S_pMSE
  }) |>
    t() |>
    set_names(c("syn1", "syn2", "syn3")) |>
    as.data.frame()
}

kl <- function(real, syn1, syn2, syn3, k) {

  realmu <- colMeans(real)
  realsd <- apply(real, 2, sd)
  real_std <- scale(real, center = realmu, scale = realsd)

  sapply(list(syn1, syn2, syn3), \(syn) {
    kldest::kld_est_nn(real_std, scale(syn, realmu, realsd), k = k)
  }) |>
    t() |>
    set_names(c("syn1", "syn2", "syn3")) |>
    as.data.frame()
}

################################################################################
## Synthetic data evaluation
################################################################################

set.seed(123)

nsim <- 1000
N <- c(100, 1000, 2000)
P <- c(5, 25, 50)

cl <- parallel::makeCluster(18)
parallel::clusterCall(cl, \(x) {
  library(densityratio)
  library(synthpop)
  library(kldest)
  library(tidyverse)
})
parallel::clusterExport(cl, c(
  "pearson", "spmse_cart", "spmse_logit", "kl"
))

out <- expand_grid(
  sim = 1:nsim,
  N = N,
  P = P
)

out <- out |>
  mutate(
    Xreal = map2(N, P, ~gen_real(.x, .y)),
    Xsyn1 = map2(Xreal, N, ~gen_syn1(.y, .x)),
    Xsyn2 = map2(Xreal, N, ~gen_syn2(.y, .x)),
    Xsyn3 = map2(Xreal, N, ~gen_syn3(.y, .x))
  )

out$out <- pbapply(out, 1, \(row) {
  list(
    pe = pearson(row$Xreal, row$Xsyn1, row$Xsyn2, row$Xsyn3),
    spmse_cart = spmse_cart(row$Xreal, row$Xsyn1, row$Xsyn2, row$Xsyn3),
    spmse_logit = spmse_logit(row$Xreal, row$Xsyn1, row$Xsyn2, row$Xsyn3, maxorder = 0),
    kl1 = kl(row$Xreal, row$Xsyn1, row$Xsyn2, row$Xsyn3, k = 1),
    kln = kl(row$Xreal, row$Xsyn1, row$Xsyn2, row$Xsyn3, k = floor(sqrt(row$N)))
  )
}, cl = cl, simplify = FALSE)

parallel::stopCluster(cl)

saveRDS(out, "results/sim2.rds")

