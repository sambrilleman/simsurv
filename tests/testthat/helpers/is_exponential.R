is_exponential <- function(mod) {
  (mod$dist == "weibull" && mod$shape == 1)
}
