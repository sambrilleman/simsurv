recover_params <- function(mod, true) {
  if (class(mod) == "phreg") {
    res <- list(lambdas = recover_lambdas(mod),
                gammas = recover_gammas(mod),
                X1 = recover_betas(mod)$X1,
                X2 = recover_betas(mod)$X2)
    res <- Filter(function(x) !is.null(x), res)
  } else if (class(mod) == "coxph") {
    res <- list(X1 = recover_betas(mod)$X1,
                X2 = recover_betas(mod)$X2,
                X1tt = recover_betas(mod)$`tt(X1)`,
                X2tt = recover_betas(mod)$`tt(X2)`)
    res <- Filter(function(x) !is.null(x), res)
  }
  return(res)
}
