library(simsurv)
stopifnot(require(testthat))
stopifnot(require(survival))
stopifnot(require(eha))

source(file.path("helpers", "sim_run.R"))
source(file.path("helpers", "check_bias.R"))
source(file.path("helpers", "recover_params.R"))
source(file.path("helpers", "recover_lambdas.R"))
source(file.path("helpers", "recover_gammas.R"))
source(file.path("helpers", "recover_betas.R"))
source(file.path("helpers", "is_exponential.R"))
source(file.path("helpers", "is_weibull.R"))
source(file.path("helpers", "is_gompertz.R"))

run_sims <- FALSE
nsims <- 100
npat <- 1000
tol <- 0.05
set.seed(9898)

#------- Standard parametric models

test_that("exponential model returns unbiased estimates", {
  testthat::skip_if_not(run_sims, "Not running simulation tests.")
  true <- list(lambdas = 0.1, X1 = -0.5, X2 = 0.2)
  sims <- sapply(seq(nsims), sim_run, npat = npat, true = true,
                 dist = "exponential", interval = c(1E-8, 10000))
  check_bias(sims = sims, true = true, tol = tol)
})

test_that("weibull model returns unbiased estimates", {
  testthat::skip_if_not(run_sims, "Not running simulation tests.")
  true <- list(lambdas = 0.1, gammas = 1.5, X1 = -0.5, X2 = 0.2)
  sims <- sapply(seq(nsims), sim_run, npat = npat, true = true,
                 dist = "weibull")
  check_bias(sims = sims, true = true, tol = tol)
})

test_that("gompertz model returns unbiased estimates", {
  testthat::skip_if_not(run_sims, "Not running simulation tests.")
  true <- list(lambdas = 0.1, gammas = .7, X1 = -0.5, X2 = 0.2)
  sims <- sapply(seq(nsims), sim_run, npat = npat, true = true,
                 dist = "gompertz")
  check_bias(sims = sims, true = true, tol = tol)
})

#------- Standard parametric models with tde (interaction with time)

tdefunction <- NULL

test_that("tde (NULL) exponential model returns unbiased estimates", {
  testthat::skip_if_not(run_sims, "Not running simulation tests.")
  true <- list(lambdas = 0.1, X1 = -0.5, X2 = 0.2,
               X1tt = 0.1, X2tt = -0.1)
  sims <- sapply(seq(nsims), sim_run, npat = npat, true = true,
                 dist = "exponential", tdefunction = tdefunction)
  check_bias(sims = sims, true = true, tol = tol, type = "bias")
})

test_that("tde (NULL) weibull model returns unbiased estimates", {
  testthat::skip_if_not(run_sims, "Not running simulation tests.")
  true <- list(lambdas = 0.1, gammas = 1.5, X1 = -0.5, X2 = 0.2,
               X1tt = 0.1, X2tt = -0.1)
  sims <- sapply(seq(nsims), sim_run, npat = npat, true = true,
                 dist = "weibull", tdefunction = tdefunction)
  check_bias(sims = sims, true = true, tol = tol)
})

test_that("tde (NULL) gompertz model returns unbiased estimates", {
  testthat::skip_if_not(run_sims, "Not running simulation tests.")
  true <- list(lambdas = 0.1, gammas = .7, X1 = -0.5, X2 = 0.2,
               X1tt = 0.1, X2tt = -0.1)
  sims <- sapply(seq(nsims), sim_run, npat = npat, true = true,
                 dist = "gompertz", tdefunction = tdefunction)
  check_bias(sims = sims, true = true, tol = tol)
})

