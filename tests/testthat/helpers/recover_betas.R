recover_betas <- function(mod) {
  nms <- names(mod$coefficients)
  nms_drop <- c("log(scale)", "log(shape)")
  sel <- setdiff(nms, nms_drop)
  if (length(sel)) as.list(mod$coefficients[sel]) else NULL
}
