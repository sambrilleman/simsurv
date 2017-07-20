recover_gammas <- function(mod) {
  if (is_exponential(mod)) {
    return(NULL)
  } else if (is_weibull(mod)) {
    return(exp(mod$coefficients[["log(shape)"]]))
  } else if (is_gompertz(mod)) {
    return(mod$coefficients[["rate"]])
  } else return(NULL)
}
