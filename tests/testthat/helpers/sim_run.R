sim_run <- function(i, npat, true, dist) {
  covs <- data.frame(id = 1:npat,
                     X1 = stats::rbinom(npat, 1L, 0.5),
                     X2 = stats::rnorm(npat, 2, 5))
  betas <- c(X1 = true$X1, X2 = true$X2)
  s1 <- simsurv(dist = dist,
                lambdas = true$lambdas,
                gammas = true$gammas,
                betas = betas, x = covs, maxt = 10)
  phreg_data <- merge(covs, s1, by = "id")
  phreg_dist <- if (dist == "exponential") "weibull" else dist
  phreg_shape <- if (dist == "exponential") 1 else 0
  mod <- phreg(Surv(eventtime, status) ~ X1 + X2, data = phreg_data,
               dist = phreg_dist, shape = phreg_shape, param = "rate")
  recover_params(mod)
}
