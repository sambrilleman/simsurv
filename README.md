
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simsurv

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/simsurv)](http://www.r-pkg.org/pkg/simsurv)
[![License](https://img.shields.io/badge/License-GPL%20%28%3E=%203%29-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html)
[![](https://cranlogs.r-pkg.org/badges/simsurv)](https://cran.r-project.org/package=simsurv)
[![](https://cranlogs.r-pkg.org/badges/grand-total/simsurv)](https://CRAN.R-project.org/package=simsurv)

**simsurv** is an R package that allows the user to simulate survival
times from standard parametric survival distributions (exponential,
Weibull, and Gompertz), 2-component mixture distributions, or a
user-defined hazard function, log hazard function, cumulative hazard
function, or log cumulative hazard function. Baseline covariates can be
included under a proportional hazards assumption. Time dependent effects
(i.e. non-proportional hazards) can be included by interacting
covariates with linear time or some transformation of time.

Under the 2-component mixture distributions the baseline survival at
time *t* is taken to be *S(t) = p x S\_1(t) + (1 - p) x S\_2(t)* where
*S\_1(t)* and *S\_2(t)* are the baseline survival under each component
of the mixture distribution and *p* is the mixing parameter. Each
component of the mixture distribution is assumed to be either
exponential, Weibull or Gompertz. The 2-component mixture distributions
can allow for a variety of flexible baseline hazard functions (see
Crowther and Lambert (2013) for some examples).

If the user wishes to provide a user-defined \[log\] \[cumulative\]
hazard function (instead of using one of the standard parametric
survival distributions) then this is also possible. If a user-defined
hazard or log hazard function is specified, then this is allowed to be
time-dependent, and the resulting cumulative hazard function does not
need to have a closed-form solution. The survival times are generated
using the approach described in Crowther and Lambert (2013), whereby the
cumulative hazard is evaluated using numerical quadrature and survival
times are generated using an iterative algorithm which nests the
quadrature-based evaluation of the cumulative hazard inside Brent’s
(1973) univariate root finder. Not requiring a closed form solution to
the cumulative hazard function has the benefit that survival times can
be generated for complex models such as flexible (spline-based) baseline
hazards or under joint longitudinal and survival models; the package
documentation and vignettes provide examples of this.

For further details on the underlying methodology and examples of usage,
please see the package vignettes (available
[here](https://cran.r-project.org/web/packages/simsurv/index.html)).

Note that this package is modelled on the user-written **survsim**
package available in the Stata software (see Crowther and Lambert
(2012)).

**Note:** Please note that the version available on GitHub is the most
up-to-date *development* version of the package. A stable version of the
package is available from the Comprehensive R Archive Network (CRAN);
type `install.packages("simsurv")` into your R session to download it.

## Getting Started

### Installation: stable release

You can install the **simsurv** package directly from CRAN. Just type
the following into your R session console:

``` r
install.packages("simsurv")
```

### Installation: development version

If, for some reason, you wish to download the development version from
GitHub, then this can be done easily using the **devtools** package. To
do this you should first check you have devtools installed by executing
the following commands from within your R session:

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

## Examples

Please also see the package
[vignettes](https://cran.r-project.org/web/packages/simsurv/index.html)
for more detailed examples and a description of the methodology
underpinning the **simsurv** package.

### Simpler examples

Generate times from a Weibull model including a binary treatment
variable, with log hazard ratio of -0.5, and censoring after 5 years:

``` r
set.seed(13579)
covs <- data.frame(id = 1:1000, trt = stats::rbinom(1000, 1L, 0.5))
s1 <- simsurv(lambdas = 0.1, gammas = 1.5,
              x = covs, betas = c(trt = -0.5), maxt = 5)
head(s1)
#>   id eventtime status
#> 1  1  5.000000      0
#> 2  2  5.000000      0
#> 3  3  2.805715      1
#> 4  4  5.000000      0
#> 5  5  5.000000      0
#> 6  6  2.276330      1
```

Generate times from a Gompertz model:

``` r
s2 <- simsurv(dist = "gompertz", lambdas = 0.1, gammas = 0.05, x = covs)
```

Generate times from a 2-component mixture Weibull model:

``` r
s3 <- simsurv(lambdas = c(0.1, 0.05), gammas = c(1, 1.5),
              mixture = TRUE, pmix = 0.5, x = covs, maxt = 5)
```

Generate times from user-defined log hazard function:

``` r
fn <- function(t, x, betas, ...)
  (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3)
s4 <- simsurv(loghazard = fn, x = covs, maxt = 1.5)
```

### Complex examples

In this section we present an example showing how the **simsurv**
package’s main function can be used to simulate survival times based on
a joint longitudinal and survival model. (In most instances, simulating
survival times from a joint longitudinal and survival model can pose
difficulties since time-dependency in the hazard function leads to there
not being a general closed form solution to the integral when evaluating
the cumulative hazard).

First we define the hazard function to pass to `simsurv`. This function
must be defined by the user and must use three named arguments: `t`,
which is the time variable; `x`, which is a data frame with each row
containing the covariate values for one individual; `betas`, which is a
data frame with each row containing the “true” parameter values for one
individual. The function must return the hazard evaluated at time `t`.
We will define the hazard function for a Weibull proportional hazards
regression submodel from a joint longitudinal and survival model with a
so-called “current value” association structure, where the hazard at
time *t* depends on the expected value of the longitudinal outcome also
evaluated at time *t*.

``` r
haz <- function(t, x, betas, ...) {
  betas[["shape"]] * (t ^ (betas[["shape"]] - 1)) * exp(
    betas[["betaEvent_intercept"]] +
    betas[["betaEvent_binary"]] * x[["Z1"]] +
    betas[["betaEvent_continuous"]] * x[["Z2"]] +
    betas[["betaEvent_assoc"]] * (
      betas[["betaLong_intercept"]] +
      betas[["betaLong_slope"]] * t +
      betas[["betaLong_binary"]] * x[["Z1"]] +
      betas[["betaLong_continuous"]] * x[["Z2"]]
    )
  )
}
```

Then we construct data frames with the true parameter values and the
covariate data for each individual:

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

Then we simulate the survival times based on the hazard function (passed
to the `hazfn` argument), the covariate values for each individual
(passed to the `x` argument), and true parameter values for each
individual (passed to the `betas` argument). We will also censor
individuals with survival times greater than 10 years (or whatever time
units we have assumed to be using) by specifying the `maxt` argument:

``` r
s1 <- simsurv(hazard = haz, x = covs, betas = betas, maxt = 10)
head(s1)
#>   id eventtime status
#> 1  1  4.328229      1
#> 2  2  8.630500      1
#> 3  3  4.608717      1
#> 4  4  1.638580      1
#> 5  5  3.183835      1
#> 6  6  7.715658      1
```

## Bug Reports

If you find any bugs, please [file an issue on
GitHub](https://github.com/sambrilleman/simsurv/issues).

## References

1.  Brilleman SL, Wolfe R, Moreno-Betancur M, and Crowther MJ.
    Simulating survival data using the simsurv R package. *Journal of
    Statistical Software* 2021; **96**(9), 1–27.
    <doi:10.18637/jss.v097.i03>.
    <https://www.jstatsoft.org/article/view/v097i03>

2.  Crowther MJ, and Lambert PC. Simulating biologically plausible
    complex survival data. *Statistics in Medicine* 2013; **32**,
    4118–4134. <doi:10.1002/sim.5823>

3.  Bender R, Augustin T, and Blettner M. Generating survival times to
    simulate Cox proportional hazards models. *Statistics in Medicine*
    2005; **24(11)**, 1713–1723.

4.  Brent R. (1973) *Algorithms for Minimization without Derivatives*.
    Englewood Cliffs, NJ: Prentice-Hall.

5.  Crowther MJ, and Lambert PC. Simulating complex survival data. *The
    Stata Journal* 2012; **12**(4), 674–687.
    <https://www.stata-journal.com/sjpdf.html?articlenum=st0275>
