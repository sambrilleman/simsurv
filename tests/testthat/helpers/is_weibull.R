is_weibull <- function(mod) {
  (mod$dist == "weibull" && mod$shape == 0)
}
