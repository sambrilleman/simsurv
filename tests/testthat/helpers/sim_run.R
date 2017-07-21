sim_run <- function(i, npat, true, dist, tdefunction = NULL) {
  covs <- data.frame(id = 1:npat,
                     X1 = stats::rbinom(npat, 1L, 0.5),
                     X2 = stats::rnorm(npat, 2, 5))
  betas <- c(X1 = true$X1, X2 = true$X2)
  tde <- c(X1 = true$X1tt, X2 = true$X2tt)
  if (is.null(tde)) { # standard parametric models
    s1 <- simsurv(dist = dist,
                  lambdas = true$lambdas,
                  gammas = true$gammas,
                  betas = betas, x = covs, maxt = 10)
    phreg_data <- merge(covs, s1, by = "id")
    phreg_dist <- if (dist == "exponential") "weibull" else dist
    phreg_shape <- if (dist == "exponential") 1 else 0
    mod <- phreg(Surv(eventtime, status) ~ X1 + X2, data = phreg_data,
                 dist = phreg_dist, shape = phreg_shape, param = "rate")
  } else { # cox model with tde (to check tde log hazard ratios)
    s1 <- simsurv(dist = dist,
                  lambdas = true$lambdas,
                  gammas = true$gammas,
                  betas = betas, x = covs, maxt = 10,
                  tde = tde, tdefunction = tdefunction)
    if (is.null(tdefunction)) {
      ttfunction <- function(x, t, ...) x * t
    } else {
      ttfunction <- function(x, t, ...) x * log(t)
    }
    cox_data <- merge(covs, s1, by = "id")
    mod <- coxph(Surv(eventtime, status) ~ X1 + X2 + tt(X1) + tt(X2),
                 data = cox_data, tt = ttfunction)
  }
  recover_params(mod)
}
