
################################################################################
## Load data and packages
################################################################################

load("data/cps5000.RData")

library(synthpop)
library(densityratio)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)

plan(multisession, workers = 16)

################################################################################
## Data wrangling and creating synthetic data
################################################################################

df <- cps |>
  select(-csp) |>
  mutate(educ = as.numeric(as.character(educ)),
         educ = case_when(educ < 39 ~ 1,
                          educ < 40 ~ 2,
                          educ < 44 ~ 3,
                          educ < 47 ~ 4) |>
           factor(labels = c("NoHS", "HS", "AoBD", "MoH")),
         race = factor(race, labels = c("White", "Non-white", "Non-white", "Non-white")),
         marital = factor(marital, labels = c("Married", "Married", "Separated",
                                              "Separated", "Widowed", "Single",
                                              "WidowedORDivorced")),
         sex = factor(sex, labels = c("Male", "Female")))

method.ini <- c("norm", "norm", "norm", "polyreg", "polyreg", "logreg", "logreg", "norm")
visit <- c(7, 1, 2, 3, 4, 5, 6, 8)
m <- 5
synlist <- list(
  trans = syn(
    df |> mutate(across(c(income, tax, age, ss), ~ .x^{1/3})),
    m = m,
    method = method.ini,
    visit.sequence = visit,
    seed = 1234,
    print.flag = FALSE
  ),
  semi = syn(
    df |> mutate(across(c(income, tax, age, ss), ~ .x^{1/3})),
    m = m,
    method = method.ini,
    visit.sequence = visit,
    semicont = list(tax = "0", ss = "0"),
    seed = 1234,
    print.flag = FALSE
  )
)

synlist$trans$syn <- synlist$trans$syn |>
  map(~ .x |> mutate(across(c(income, tax, age, ss), ~.x^3)))
synlist$semi$syn <- synlist$semi$syn |>
  map(~ .x |> mutate(across(c(income, tax, age, ss), ~.x^3)))


################################################################################
## Calculate density ratio estimates (Pearson divergence after ulsif)
################################################################################

vars <- c("age", "income", "ss", "tax")

## Calculate density ratio estimates for each variable

PE <- future_map(
  synlist, function(syndat) {
    future_map(vars, function(var) {
      future_map_dbl(syndat$syn, function(syndati) {
        ulsif(
          df[[var]],
          syndati[[var]],
          progressbar = FALSE
        ) |>
          summary(test = FALSE) |>
          (\(x) x$PE)()
      }, .options = furrr_options(seed = TRUE), .progress = TRUE)
    }, .options = furrr_options(seed = TRUE), .progress = TRUE
    )
  }, .options = furrr_options(seed = TRUE), .progress = TRUE
)

# Calculate density ratio estimates over entire data set (all variables
# simultaneously)

PE_allvars <- synlist |>
  future_map(function(x) {
    future_map(x$syn, ~ulsif(df, .x, progressbar = FALSE),
               .options = furrr_options(seed = TRUE),
               .progress = TRUE)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

## Check values
PE |> map(\(syn) map_dbl(syn, mean)) |> data.frame()
PE_allvars |> map(~map_dbl(.x, ~summary(.x) |> {\(x) x$PE}()) |> mean())

################################################################################
## Reweighting using density ratio estimates
################################################################################

true_fit <- lm(
  log(income) ~ . - tax - ss + log(tax+1) + log(ss+1), df
)
summary(true_fit)


true_coefs <- coef(true_fit)
true_ses <- vcov(true_fit) |> diag() |> sqrt()

syn_coefs_semi <- lapply(
  synlist$semi$syn,
  \(syn) lm(
    log(income) ~ . - tax - ss + log(tax+1) + log(ss+1), syn) |> coef()
)

syn_coefs_adj <- map2(PE_allvars$semi, synlist$semi$syn, \(w, syn) {
  lm(log(income) ~ . - tax - ss + log(tax+1) + log(ss+1),
     syn,
     weights = pmax(0, predict(w, syn))
     ) |> coef()
})

syn_coefs_cart <- lapply(
  synlist$smoothed$syn,
  \(syn) lm(log(income) ~ . - tax - ss + log(tax+1) + log(ss+1), syn) |> coef()
)

coefs <- data.frame(
  true = true_coefs,
  syn_semi = do.call(cbind, syn_coefs_semi) |> rowMeans(),
  syn_adj = do.call(cbind, syn_coefs_adj) |> rowMeans()
) |>
  round(5)


bias <- data.frame(
  syn_semi = (do.call(cbind, syn_coefs_semi) |> rowMeans() - true_coefs) / true_ses,
  syn_adj = (do.call(cbind, syn_coefs_adj) |> rowMeans() - true_coefs) / true_ses
)


save(list = c("df", "synlist", "PE", "PE_allvars",
               "coefs", "bias"),
     file = "results/application.RData")

