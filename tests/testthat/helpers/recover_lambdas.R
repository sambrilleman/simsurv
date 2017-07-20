recover_lambdas <- function(mod) {
  if (is_exponential(mod)) {
    scale <- exp(mod$coefficients[["log(scale)"]])
    return(scale ^ (-1)) # transform scale parameter for exponential
  } else if (is_weibull(mod)) {
    scale <- exp(mod$coefficients[["log(scale)"]])
    shape <- exp(mod$coefficients[["log(shape)"]])
    return(scale ^ (-shape)) # transform scale parameter for weibull
  } else if (is_gompertz(mod)) {
    level <- exp(mod$coefficients[["log(level)"]])
    return(level)
  } else return(NULL)
}
