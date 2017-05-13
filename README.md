
<!-- README.md is generated from README.Rmd. Please edit that file -->
simsurv
=======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/simsurv)](http://www.r-pkg.org/pkg/simsurv) [![License](https://img.shields.io/badge/License-GPL%20%28%3E=%203%29-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html)

**simsurv** is an R package that allows the user to simulate survival times from any user-specified hazard function. The hazard function is allowed to be time-dependent, and the resulting cumulative hazard function does not need to have a closed-form solution. The cumulative hazard is evaluated using Gauss-Kronrod quadrature and survival times are generated using a combination of the method in Bender et al. (2005) and Brent's (1973) univariate root finder. Not requiring a closed form solution to the cumulative hazard function has the benefit that survival times can be generated for complex models such as joint longitudinal and survival models; the package documentation provides an example of this.

**Note:** Please note that the version available on GitHub is the most up-to-date *development* version of the package. A stable version of the package will be available from CRAN once it is released.

Getting Started
---------------

### Installation

You can install **simsurv** directly from GitHub using the **devtools** package. To do this you should first check you have devtools installed by executing the following commands from within your R session:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}
```

Then execute the following commands to install **simsurv**:

``` r
library(devtools)
install_github("sambrilleman/simsurv")
```

Example
-------

In this section we present an example showing how the **simsurv** package's main function can be used to simulate survival times based on a joint longitudinal and survival model. (In most instances, simulating survival times from a joint longitudinal and survival model can pose difficulties since time-dependency in the hazard function leads to there not being a general closed form solution to the integral when evaluating the cumulative hazard).

First we define the hazard function to pass to `simsurv`. We will simulate data based on a Weibull proportional hazards regression submodel from a joint longitudinal and survival model with a so-called "current value" association structure, where the hazard at time *t* depends on the expected value of the longitudinal outcome also evaluated at time *t*.

``` r
weibull_ph_hazfn <- function(t, x, pars) {
  pars[["shape"]] * (t ^ (pars[["shape"]] - 1)) * exp(
    pars[["betaEvent_intercept"]] +
    pars[["betaEvent_binary"]] * x[["Z1"]] +
    pars[["betaEvent_continuous"]] * x[["Z2"]] +
    pars[["betaEvent_assoc"]] * (
      pars[["betaLong_intercept"]] +
      pars[["betaLong_slope"]] * t +
      pars[["betaLong_binary"]] * x[["Z1"]] +
      pars[["betaLong_continuous"]] * x[["Z2"]]
    )
  )
}
```

Then we construct data frames with the true parameter values and the covariate data for each individual:

``` r
set.seed(5454) # set seed before simulating data
N <- 20        # number of individuals

# Population (fixed effect) parameters
betas <- data.frame(
  shape                = rep(2,    N),
  betaEvent_intercept  = rep(-11.9,N),
  betaEvent_binary     = rep(0.6,  N),
  betaEvent_continuous = rep(0.08, N),
  betaEvent_assoc      = rep(0.03, N),
  betaLong_binary      = rep(-1.5, N),
  betaLong_continuous  = rep(1,    N),
  betaLong_intercept   = rep(90,   N),
  betaLong_slope       = rep(2.5,  N)
)

# Individual-specific (random effect) parameters
b_corrmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
b_sds     <- c(20, 3)
b_means   <- rep(0, 2)
b_z       <- MASS::mvrnorm(n = N, mu = b_means, Sigma = b_corrmat)
b         <- sapply(1:length(b_sds), function(x) b_sds[x] * b_z[,x])
betas$betaLong_intercept <- betas$betaLong_intercept + b[,1]
betas$betaLong_slope     <- betas$betaLong_slope     + b[,2]

# Covariate data
covs <- data.frame(
  Z1 = rbinom(N, 1, 0.45), # a binary covariate
  Z2 = rnorm(N, 44, 8.5)   # a continuous covariate
)
```

Then we simulate the survival times based on the hazard function (passed to the `hazfn` argument), the covariate values for each individual (passed to the `x` argument), and true parameter values for each individual (passed to the `pars` argument). We will also censor individuals with survival times greater than 10 years (or whatever time units we have assumed to be using) by specifying the `maxt` argument:

``` r
s1 <- simsurv(hazfn = weibull_ph_hazfn, x = covs, pars = betas, maxt = 10)
head(s1)
#>   eventtime status
#> 1  5.504233      1
#> 2  5.946418      1
#> 3  6.618139      1
#> 4  2.753683      1
#> 5  3.884597      1
#> 6  9.295474      1
```

Bug Reports
-----------

If you find any bugs, please report them via email to [Sam Brilleman](mailto:sam.brilleman@monash.edu).

References
----------

1.  Bender R, Augustin T, and Blettner M. Generating survival times to simulate Cox proportional hazards models. *Statistics in Medicine* 2005; **24(11)**, 1713-1723.

2.  Brent R. (1973) *Algorithms for Minimization without Derivatives*. Englewood Cliffs, NJ: Prentice-Hall.
