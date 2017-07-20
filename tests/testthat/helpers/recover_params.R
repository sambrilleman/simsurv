recover_params <- function(mod, true) {
  if (class(mod) == "phreg") {
    res <- list(lambdas = recover_lambdas(mod),
                gammas = recover_gammas(mod),
                X1 = recover_betas(mod)$X1,
                X2 = recover_betas(mod)$X2)
    res <- Filter(function(x) !is.null(x), res)
  }
  return(res)
}
