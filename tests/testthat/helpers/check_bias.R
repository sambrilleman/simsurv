check_bias <- function(sims, true, type = "relbias", tol = 0.05) {
  for (i in rownames(sims)) {
    sims_i <- unlist(sims[i,])
    true_i <- true[[i]]
    if (type == "bias") {
      bias <- mean(sims_i - true_i)
    } else if (type == "relbias") {
      bias <- mean((sims_i - true_i) / true_i)
    }
    expect_equal(bias, 0.0, tol = tol)
  }
}
