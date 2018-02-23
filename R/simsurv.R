#' Simulate survival data
#'
#' Simulate survival times from standard parametric survival distributions,
#' 2-component mixture distributions, or a user-defined hazard or log hazard
#' function.
#'
#' @export
#' @importFrom methods is
#'
#' @param dist Character string specifying the parametric survival distribution.
#'   Can be \code{"weibull"} (the default), \code{"exponential"}, or
#'   \code{"gompertz"}. Note that this is ignored if a user-defined baseline
#'   hazard function is provided via the \code{hazard} or \code{loghazard}
#'   argument.
#' @param lambdas A numeric vector corresponding to the scale parameters for
#'   the exponential, Weibull or Gompertz distributions. This vector should be
#'   length one for any of the standard parametric distributions, or length two
#'   for a 2-component mixture distribution when \code{mixture = TRUE}.
#' @param gammas A numeric vector corresponding to the shape parameters for
#'   the Weibull or Gompertz distributions. This vector should be length
#'   one for the standard Weibull or Gompertz distributions, or length two
#'   for a 2-component mixture distribution when \code{mixture = TRUE}.
#' @param x A data frame containing the covariate values for each individual.
#'   Each row of the data frame should supply covariate data for one individual.
#'   The column names should correspond to the named elements in \code{betas}.
#'   If no covariates are being used to simulate the event times, then this data
#'   frame should just include an ID value for each individual; this is necessary
#'   since the number of individuals to simulate event times for is taken to be the
#'   number of rows in \code{x} (unless \code{idvar} and \code{ids} are
#'   specified, in which case the number of individuals to simulate event times
#'   for is taken to be \code{length(ids)}).
#' @param betas A named vector, a data frame, or a list of data frames,
#'   containing the "true" parameters (e.g. log hazard ratios).
#'   If a standard exponential, Weibull, or Gompertz distribution, or a 2-component
#'   mixture distribution is being used then \code{betas} should provide the log hazard
#'   ratios for each of the baseline covariates that is to be included in the linear
#'   predictor of the proportional hazards model. This is most easily specified
#'   as a named vector (see the \strong{Examples}). Alternatively, if a
#'   user-defined baseline hazard (or log hazard) function is provided via
#'   the \code{hazard} (or \code{loghazard}) argument, then \code{betas} can be
#'   specified as a vector, a data frame, or a list of data frames, and the
#'   of user-defined hazard (or log hazard) function should extract named elements
#'   \code{betas} however necessary in order to calculate the hazard (or log
#'   hazard) for each individual. See the \strong{Examples}.
#' @param tde A named vector, containing the "true" parameters that will be used
#'   to create time dependent effects (i.e. non-proportional hazards). The
#'   values specified in \code{tde} are used as coefficients (in the linear
#'   predictor of the proportional hazards model) on an interaction term between
#'   the correponding covariate and time (or some function of time, for example
#'   log time, if \code{tdefunction} is not \code{NULL}).
#' @param tdefunction An optional function of time to which covariates specified
#'   in \code{tde} will be interacted, in order to generate time dependent
#'   effects (i.e. non-proportional hazards). If \code{NULL} then the covariates
#'   specified in \code{tde} will be interacted with linear time. This can be the
#'   character string corresponding to a standard function (such as "log") or it
#'   can be a user-defined function with a single argument (for example
#'   \code{function(x) x ^ 2}).
#' @param mixture Logical specifying whether to use a 2-component mixture
#'   model for the survival distribution. If \code{TRUE}, then the distribution
#'   of the mixture components is determined by the \code{dist} argument.
#' @param pmix Scalar between 0 and 1 defining the mixing parameter when
#'   \code{mixture = TRUE}. The baseline survival at time t is taken to be
#'   \eqn{S(t) = p * S_1(t) + (1 - p) * S_2(t)}
#'   where \eqn{S_1(t)} and \eqn{S_2(t)} are the baseline survival under each
#'   component of the mixture distribution.
#' @param hazard Optionally, a user-defined hazard function, with arguments
#'   \code{t}, \code{x}, and \code{betas}. This function should return the
#'   hazard at time \code{t} for an individual with covariates supplied via \code{x}
#'   and parameters supplied via \code{betas}. See the \strong{Examples}.
#' @param loghazard Optionally, a user-defined log hazard function, with arguments
#'   \code{t}, \code{x}, and \code{betas}. This function should return the
#'   log hazard at time \code{t} for an individual with covariates supplied via
#'   \code{x} and parameters supplied via \code{betas}. See the \strong{Examples}.
#' @param cumhazard Optionally, a user-defined cumulative hazard function, with
#'   arguments \code{t}, \code{x}, and \code{betas}. This function should return the
#'   cumulative hazard at time \code{t} for an individual with covariates supplied
#'   via \code{x} and parameters supplied via \code{betas}. See the \strong{Examples}.
#' @param logcumhazard Optionally, a user-defined log cumulative hazard function, with
#'   arguments \code{t}, \code{x}, and \code{betas}. This function should return the
#'   log cumulative hazard at time \code{t} for an individual with covariates supplied
#'   via \code{x} and parameters supplied via \code{betas}. See the \strong{Examples}.
#' @param idvar The name of the ID variable identifying individuals. This is
#'   only required when \code{x} (and \code{betas} if it is supplied as a
#'   data frame) contains multiple rows per individual. Otherwise, if
#'   \code{idvar = NULL} then each row of \code{x} (and \code{betas} if it is
#'   supplied as a data frame) is assumed to correspond to a different individual.
#' @param ids A vector containing the unique values of \code{idvar} (i.e. the
#'   unique individual IDs). This is only required when \code{x} (and \code{betas}
#'   if it is supplied as a data frame) contain multiple rows per individual. Otherwise,
#'   if \code{idvar = NULL} then each row of \code{x} (and \code{betas} if it is
#'   supplied as a data frame) is assumed to correspond to a different individual.
#' @param maxt The maximum event time. For simulated event times greater than
#'   \code{maxt}, the event time (\code{"eventtime"}) returned in the data frame
#'   will be truncated at \code{maxt} and the event indicator (\code{"status"})
#'   will be set to zero indicating that the individual was right-censored.
#' @param nodes Integer specifying the number of quadrature nodes to use for
#'   the Gauss-Kronrod quadrature. Can be 7, 11, or 15.
#' @param interval The interval over which to search for the
#'   \code{\link{uniroot}} corresponding to each simulated event time.
#' @param seed The \code{\link[=set.seed]{seed}} to use.
#' @param ... Other arguments passed to \code{hazard} or \code{loghazard}.
#'
#' @details The \code{simsurv} function simulates survival times from
#' standard parametric survival distributions (exponential, Weibull, Gompertz),
#' 2-component mixture distributions, or a user-defined hazard or log hazard function.
#' Baseline covariates can be included under a proportional hazards assumption.
#' Time dependent effects (i.e. non-proportional hazards) can be included by
#' interacting covariates with time (by specifying them in the \code{tde}
#' argument); the default behaviour is to interact the covariates with linear
#' time, however, they can be interacted with some other function of time simply
#' by using the \code{tdefunction} argument.
#'
#' Under the 2-component mixture distributions (obtained by setting
#' \code{mixture = TRUE}) the baseline survival at time \eqn{t} is taken to be
#' \eqn{S(t) = p * S_1(t) + (1 - p) * S_2(t)}
#' where \eqn{S_1(t)} and \eqn{S_2(t)} are the baseline survival under each
#' component of the mixture distribution and \eqn{p} is the mixing parameter
#' specified via the argument \code{pmix}. Each component of the mixture
#' distribution is assumed to be either exponential, Weibull or Gompertz.
#' The 2-component mixture distributions can allow for a variety of flexible
#' baseline hazard functions (see Crowther and Lambert (2013) for some examples).
#'
#' If the user wishes to provide a user-defined hazard or log hazard function
#' (instead of using one of the standard parametric survival distributions) then
#' this is also possible via the \code{hazard} or \code{loghazard} argument.
#' If a user-defined hazard or log hazard function is specified, then this is
#' allowed to be time-dependent, and the resulting cumulative hazard function
#' does not need to have a closed-form solution. The survival times are
#' generated using the approach described in Crowther and Lambert (2013),
#' whereby the cumulative hazard is evaluated using numerical quadrature and
#' survival times are generated using an iterative algorithm which nests the
#' quadrature-based evaluation of the cumulative hazard inside Brent's (1973)
#' univariate root finder (for the latter the \code{\link{uniroot}}
#' function is used). Not requiring a closed form solution to the cumulative
#' hazard function has the benefit that survival times can be generated for
#' complex models such as joint longitudinal and survival models; the
#' \strong{Examples} section provides an example of this.
#'
#' \subsection{Parameterisations for the Weibull distribution}{
#' For the exponential and Weibull distributions, with scale parameter
#' \eqn{lambda} and shape parameter \eqn{gamma} (with \eqn{gamma} fixed equal
#' to 1 for the exponential distribution) the baseline hazard and survival
#' functions used by \code{simsurv} are:
#' \eqn{h(t) = gamma * lambda * t ^ {gamma - 1}} and
#' \eqn{S(t) = exp(-lambda * t ^ {gamma})}.
#'
#' Note that this parameterisation differs from the one used by
#' \code{\link{dweibull}} or the \code{\link[eha]{phreg}} modelling
#' function in the \pkg{eha} package. The parameterisation used in those
#' functions can be achieved by transforming the scale parameter via the
#' relationship \eqn{b = lambda ^ {-1 / gamma}}, or equivalently
#' \eqn{lambda = b ^ {-gamma}} where \eqn{b} is the scale parameter under
#' the parameterisation of the Weibull distribution used by
#' \code{\link{dweibull}} or \code{\link[eha]{phreg}}.
#' }
#'
#' @note This package is modelled on the user-written \code{survsim} package
#'   available in the Stata software (see Crowther and Lambert (2012)).
#'
#' @return A data frame with a row for each individual, and the following three
#'   columns:
#'   \itemize{
#'     \item \code{id}  The individual identifier
#'     \item \code{eventtime} The simulated event (or censoring) time
#'     \item \code{status} The event indicator, 1 for failure, 0 for censored
#'   }
#'
#' @author Sam Brilleman (\email{sam.brilleman@@monash.edu})
#'
#' @references
#'   Crowther MJ, and Lambert PC. (2013) Simulating biologically plausible
#'   complex survival data. \emph{Statistics in Medicine} \strong{32},
#'   4118--4134. \doi{10.1002/sim.5823}
#'
#'   Bender R, Augustin T, and Blettner M. (2005) Generating survival times to
#'   simulate Cox proportional hazards models. \emph{Statistics in Medicine}
#'   \strong{24}(11), 1713--1723.
#'
#'   Brent R. (1973) \emph{Algorithms for Minimization without Derivatives}.
#'   Englewood Cliffs, NJ: Prentice-Hall.
#'
#'   Crowther MJ, and Lambert PC. (2012) Simulating complex survival data.
#'   \emph{The Stata Journal} \strong{12}(4), 674--687.
#'   \url{http://www.stata-journal.com/sjpdf.html?articlenum=st0275}
#'
#' @examples
#'   #-------------- Simpler examples
#'
#'   # Generate times from a Weibull model including a binary
#'   # treatment variable, with log(hazard ratio) = -0.5, and censoring
#'   # after 5 years:
#'   covs <- data.frame(id = 1:100, trt = stats::rbinom(100, 1L, 0.5))
#'   s1 <- simsurv(lambdas = 0.1, gammas = 1.5, betas = c(trt = -0.5),
#'                 x = covs, maxt = 5)
#'   head(s1)
#'
#'   # Generate times from a Gompertz model:
#'   s2 <- simsurv(dist = "gompertz", lambdas = 0.1, gammas = 0.05, x = covs)
#'
#'   # Generate times from a 2-component mixture Weibull model:
#'   s3 <- simsurv(lambdas = c(0.1, 0.05), gammas = c(1, 1.5),
#'                 mixture = TRUE, pmix = 0.5, x = covs, maxt = 5)
#'
#'   # Generate times from user-defined log hazard function:
#'   fn <- function(t, x, betas, ...)
#'     (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3)
#'   s4 <- simsurv(loghazard = fn, x = covs, maxt = 1.5)
#'
#'   # Generate times from a Weibull model with diminishing treatment effect:
#'   s5 <- simsurv(lambdas = 0.1, gammas = 1.5, betas = c(trt = -0.5),
#'                 x = covs, tde = c(trt = 0.05), tdefunction = "log")
#'
#'   #-------------- Complex examples
#'
#'   # Here we present an example of simulating survival times
#'   # based on a joint longitudinal and survival model
#'
#'   # First we define the hazard function to pass to simsurv
#'   # (NB this is a Weibull proportional hazards regression submodel
#'   # from a joint longitudinal and survival model with a "current
#'   # value" association structure).
#'   haz <- function(t, x, betas, ...) {
#'       betas[["shape"]] * (t ^ (betas[["shape"]] - 1)) * exp(
#'         betas[["betaEvent_intercept"]] +
#'         betas[["betaEvent_binary"]] * x[["Z1"]] +
#'         betas[["betaEvent_continuous"]] * x[["Z2"]] +
#'         betas[["betaEvent_assoc"]] * (
#'           betas[["betaLong_intercept"]] +
#'           betas[["betaLong_slope"]] * t +
#'           betas[["betaLong_binary"]] * x[["Z1"]] +
#'           betas[["betaLong_continuous"]] * x[["Z2"]]
#'         )
#'       )
#'   }
#'
#'   # Then we construct data frames with the true parameter
#'   # values and the covariate data for each individual
#'   set.seed(5454) # set seed before simulating data
#'   N <- 20        # number of individuals
#'
#'   # Population (fixed effect) parameters
#'   betas <- data.frame(
#'     shape                = rep(2,    N),
#'     betaEvent_intercept  = rep(-11.9,N),
#'     betaEvent_binary     = rep(0.6,  N),
#'     betaEvent_continuous = rep(0.08, N),
#'     betaEvent_assoc      = rep(0.03, N),
#'     betaLong_binary      = rep(-1.5, N),
#'     betaLong_continuous  = rep(1,    N),
#'     betaLong_intercept   = rep(90,   N),
#'     betaLong_slope       = rep(2.5,  N)
#'   )
#'
#'   # Individual-specific (random effect) parameters
#'   b_corrmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#'   b_sds     <- c(20, 3)
#'   b_means   <- rep(0, 2)
#'   b_z       <- MASS::mvrnorm(n = N, mu = b_means, Sigma = b_corrmat)
#'   b         <- sapply(1:length(b_sds), function(x) b_sds[x] * b_z[,x])
#'   betas$betaLong_intercept <- betas$betaLong_intercept + b[,1]
#'   betas$betaLong_slope     <- betas$betaLong_slope     + b[,2]
#'
#'   # Covariate data
#'   covdat <- data.frame(
#'     Z1 = stats::rbinom(N, 1, 0.45), # a binary covariate
#'     Z2 = stats::rnorm(N, 44, 8.5)   # a continuous covariate
#'   )
#'
#'   # Then we simulate the survival times based on the
#'   # hazard function, covariates, and true parameter values
#'   times <- simsurv(hazard = haz, x = covdat, betas = betas, maxt = 10)
#'   head(times)
#'
simsurv <- function(dist = c("weibull", "exponential", "gompertz"),
                    lambdas, gammas, x, betas, tde, tdefunction = NULL,
                    mixture = FALSE, pmix = 0.5, hazard, loghazard,
                    cumhazard, logcumhazard,
                    idvar = NULL, ids = NULL, nodes = 15,
                    maxt = NULL, interval = c(1E-8, 500),
                    seed = sample.int(.Machine$integer.max, 1), ...) {
  set.seed(seed)
  dist <- match.arg(dist)
  if (!is.numeric(interval) || !length(interval) == 2)
    stop("'interval' should a length 2 numeric vector.")
  if (!all(interval > 0))
    stop("Both 'interval' limits must be positive.")
  if (!is.null(maxt) && (interval[2] <= maxt))
    stop("The upper limit of 'interval' must be greater than 'maxt'.")
  if (missing(lambdas))
    lambdas <- NULL
  if (missing(gammas))
    gammas <- NULL
  if (missing(betas))
    betas <- NULL
  if (missing(tde))
    tde <- NULL
  if (missing(hazard))
    hazard <- NULL
  if (missing(loghazard))
    loghazard <- NULL
  if (missing(cumhazard))
    cumhazard <- NULL
  if (missing(logcumhazard))
    logcumhazard <- NULL
  user_defined <- length(c(hazard, loghazard, cumhazard, logcumhazard))
  if (user_defined > 1L)
    stop("Only one of 'hazard', 'loghazard', 'cumhazard' ",
         "or 'logcumhazard' can be specified.")
  ok_args <- c("t", "x", "betas")
  if (!is.null(hazard)) {
    if (!is.function(hazard))
      stop("'hazard' should be a function.")
    nm_args <- names(formals(hazard))
    if (!all(ok_args %in% nm_args))
      stop("'hazard' function should have the following named arguments: ",
           paste(ok_args, collapse = ", "))
  }
  if (!is.null(loghazard)) {
    if (!is.function(loghazard))
      stop("'loghazard' should be a function.")
    nm_args <- names(formals(loghazard))
    if (!all(ok_args %in% nm_args))
      stop("'loghazard' function should have the following named arguments: ",
           paste(ok_args, collapse = ", "))
  }
  if (!is.null(cumhazard)) {
    if (!is.function(cumhazard))
      stop("'cumhazard' should be a function.")
    nm_args <- names(formals(cumhazard))
    if (!all(ok_args %in% nm_args))
      stop("'cumhazard' function should have the following named arguments: ",
           paste(ok_args, collapse = ", "))
  }
  if (!is.null(logcumhazard)) {
    if (!is.function(logcumhazard))
      stop("'logcumhazard' should be a function.")
    nm_args <- names(formals(logcumhazard))
    if (!all(ok_args %in% nm_args))
      stop("'logcumhazard' function should have the following named arguments: ",
           paste(ok_args, collapse = ", "))
  }
  x <- validate_x(x)
  if (!is.null(betas)) {
    betas <- validate_betas(betas)
    if (!user_defined && !is(betas, "list") && is.vector(betas) &&
      (!all(names(betas) %in% colnames(x))))
        stop("The named elements in 'betas' should correspond to named ",
             "columns in the data frame specified in 'x'.")
  }
  if (!is.null(tde)) {
    tde <- validate_tde(tde)
    if (!all(names(tde) %in% colnames(x)))
      stop("The named elements in 'tde' should correspond to named ",
           "columns in the data frame specified in 'x'.")
  }

  if (!is.null(ids) == is.null(idvar))
    stop("Both 'idvar' and 'ids' must be supplied together.")
  if (!is.null(ids)) {
    N <- length(ids) # number of individuals
    if (any(duplicated(ids)))
      stop("The 'ids' vector must specify unique ID values.")
  } else {
    N <- nrow(x) # number of individuals
    ids <- seq(N)
  }
  if (is.null(hazard) &&
      is.null(loghazard) &&
      is.null(cumhazard) &&
      is.null(logcumhazard) &&
      is.null(tde)) { # standard parametric or mixture distribution
    if (!mixture) { # non-mixture dist; use analytic form for inverted survival
      inverted_surv <- get_inverted_surv(dist = dist,
                                         lambdas = lambdas, gammas = gammas)
      tt <- sapply(ids, function(i) {
        x_i <- subset_df(x, i, idvar = idvar)
        betas_i <- subset_df(betas, i, idvar = idvar)
        u_i <- stats::runif(1)
        t_i <- inverted_surv(u = u_i, x = x_i, betas = betas_i)
        return(t_i)
      })
    } else { # mixture dist; use root finding to invert survival
      survival <- get_survival(dist = dist, lambdas = lambdas, gammas = gammas,
                               mixture = mixture, pmix = pmix)
      tt <- sapply(ids, function(i) {
        x_i <- subset_df(x, i, idvar = idvar)
        betas_i <- subset_df(betas, i, idvar = idvar)
        u_i <- stats::runif(1)
        # check whether S(t) is still greater than random uniform variable u_i at the
        # upper limit of uniroot's interval (otherwise uniroot will return an error)
        at_limit <- rootfn_surv(interval[2], survival = survival, x = x_i,
                                betas = betas_i, u = u_i, ...)
        if (is.nan(at_limit)) {
          STOP_nan_at_limit()
        } else if (at_limit > 0) {
          if (is.null(maxt)) { # no censoring time
            STOP_increase_limit()
          } else { # individual will be censored anyway, so just return interval[2]
            return(interval[2])
          }
        } else {
          t_i <- stats::uniroot(
            rootfn_surv, survival = survival, x = x_i, betas = betas_i,
            u = u_i, ..., interval = interval)$root
        }
        return(t_i)
      })
    }
  } else if (is.null(hazard) &&
             is.null(loghazard) &&
             is.null(cumhazard) &&
             is.null(logcumhazard)) { # # standard parametric with tde
    if (is.null(tdefunction)) {
      tdefunction <- function(x) x # just returns input value (ie. time)
    } else if (!is.null(tdefunction)) {
      tdefunction <- tryCatch(match.fun(tdefunction), error = function(e)
        stop("'tdefunction' should be a function.", call. = FALSE))
    }
    hazard <- get_hazard(dist = dist, lambdas = lambdas, gammas = gammas,
                         mixture = mixture, pmix = pmix, tde = tde,
                         tdefunction = tdefunction)
    qq <- get_quadpoints(nodes)
    tt <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas, i, idvar = idvar)
      u_i <- stats::runif(1)
      # check whether S(t) is still greater than random uniform variable u_i at the
      # upper limit of uniroot's interval (otherwise uniroot will return an error)
      at_limit <- rootfn_hazard(interval[2], hazard = hazard, x = x_i,
                                betas = betas_i, u = u_i, qq = qq,
                                tde = tde, tdefunction = tdefunction)
      if (is.nan(at_limit)) {
        STOP_nan_at_limit()
      } else if (at_limit > 0) {
        if (is.null(maxt)) { # no censoring time
          STOP_increase_limit()
        } else { # individual will be censored anyway, so just return interval[2]
          return(interval[2])
        }
      } else {
        t_i <- stats::uniroot(
          rootfn_hazard, hazard = hazard, x = x_i, betas = betas_i,
          u = u_i, qq = qq, tde = tde, tdefunction = tdefunction,
          interval = interval)$root
      }
      return(t_i)
    })
  } else if (!is.null(cumhazard) ||
             !is.null(logcumhazard)) { # user-defined cum hazard or log cum hazard
    if (!is.null(tde))
      stop("'tde' cannot be specified with a user-defined [log] cumulative hazard ",
           "function; please just incorporate the time dependent effects ",
           "into the [log] cumulative hazard function you are defining.")
    if (!is.null(logcumhazard)) {
      cumhazard <- function(t, x, betas, ...) exp(logcumhazard(t, x, betas, ...))
    }
    tt <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas, i, idvar = idvar)
      u_i <- stats::runif(1)
      # check whether S(t) is still greater than random uniform variable u_i at the
      # upper limit of uniroot's interval (otherwise uniroot will return an error)
      at_limit <- rootfn_cumhazard(interval[2], cumhazard = cumhazard, x = x_i,
                                   betas = betas_i, u = u_i, ...)
      if (is.nan(at_limit)) {
        STOP_nan_at_limit()
      } else if (at_limit > 0) {
        if (is.null(maxt)) { # no censoring time
          STOP_increase_limit()
        } else { # individual will be censored anyway, so just return interval[2]
          return(interval[2])
        }
      } else {
        t_i <- stats::uniroot(
          rootfn_cumhazard, cumhazard = cumhazard, x = x_i, betas = betas_i,
          u = u_i, ..., interval = interval)$root
      }
      return(t_i)
    })
  } else { # user-defined hazard or log hazard
    if (!is.null(tde))
      stop("'tde' cannot be specified with a user-defined [log] hazard ",
           "function; please just incorporate the time dependent effects ",
           "into the [log] hazard function you are defining.")
    if (!is.null(loghazard)) {
      hazard <- function(t, x, betas, ...) exp(loghazard(t, x, betas, ...))
    }
    qq <- get_quadpoints(nodes)
    tt <- sapply(ids, function(i) {
      x_i <- subset_df(x, i, idvar = idvar)
      betas_i <- subset_df(betas, i, idvar = idvar)
      u_i <- stats::runif(1)
      # check whether S(t) is still greater than random uniform variable u_i at the
      # upper limit of uniroot's interval (otherwise uniroot will return an error)
      at_limit <- rootfn_hazard(interval[2], hazard = hazard, x = x_i,
                                betas = betas_i, u = u_i, qq = qq, ...)
      if (is.nan(at_limit)) {
        STOP_nan_at_limit()
      } else if (at_limit > 0) {
        if (is.null(maxt)) { # no censoring time
          STOP_increase_limit()
        } else { # individual will be censored anyway, so just return interval[2]
          return(interval[2])
        }
      } else {
        t_i <- stats::uniroot(
          rootfn_hazard, hazard = hazard, x = x_i, betas = betas_i,
          u = u_i, qq = qq, ..., interval = interval)$root
      }
      return(t_i)
    })
  }
  if (!is.null(maxt)) {
    if (maxt <= 0)
      stop("'maxt' must be positive.")
    d <- as.integer(tt < maxt)
    tt <- tt * d + maxt * (1 - d)
  } else {
    d <- rep(1, N)
  }
  ret <- data.frame(id = if (!is.null(ids)) ids else seq(N),
                    eventtime = tt, status = d, row.names = NULL)
  if (!is.null(idvar))
    colnames(ret)[[1L]] <- idvar
  ret <- structure(ret, seed = seed)
  return(ret)
}


#----------- internal

# Return the inverted survival function, rearranged to solve for t
#
# @param dist The name of the distribution for the baseline hazard
# @param lambdas The scale parameter(s) for the baseline hazard
# @param gammas The shape parameter(s) for the baseline hazard
get_inverted_surv <- function(dist = c("weibull", "exponential", "gompertz"),
                              lambdas = NULL, gammas = NULL) {
  validate_lambdas(lambdas = lambdas, dist = dist, mixture = FALSE)
  validate_gammas(gammas = gammas, dist = dist, mixture = FALSE)
  if (dist == "weibull") { # weibull inverted survival
    inv_surv <- function(u, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      t <- (-log(u) / (lambdas[1L] * exp(eta))) ^ (1 / gammas[1L])
      return(t)
    }
  } else if (dist == "exponential") { # exponential inverted survival
    inv_surv <- function(u, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      t <- (-log(u) / (lambdas[1L] * exp(eta)))
      return(t)
    }
  } else if (dist == "gompertz") { # gompertz inverted survival
    inv_surv <- function(u, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      t <- (1 / gammas[1L]) * log(((-gammas[1L]*log(u)) / (lambdas[1L]*exp(eta))) + 1)
      return(t)
    }
  }
  return(inv_surv)
}

# Return the survival function for the specified distribution
#
# @param dist The name of the distribution for the baseline hazard
# @param lambdas The scale parameter(s) for the baseline hazard
# @param gammas The shape parameter(s) for the baseline hazard
# @param mixture Logical specifying whether to use a 2-component mixture model
#   for the baseline hazard distribution
# @param pmix Scalar specifying the mixture parameter, must be between 0 and 1
get_survival <- function(dist = c("weibull", "exponential", "gompertz"),
                         lambdas = NULL, gammas = NULL, mixture = FALSE, pmix = 0.5) {
  validate_lambdas(lambdas = lambdas, dist = dist, mixture = mixture)
  validate_gammas(gammas = gammas, dist = dist, mixture = mixture)
  if (pmix < 0 || pmix > 1)
    stop("'pmix' must be between 0 and 1.", call. = FALSE)
  if (dist == "weibull" && !mixture) { # weibull survival
    surv <- function(t, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      basesurv <- exp(-lambdas[1L] * (t ^ gammas[1L]))
      return(basesurv ^ exp(eta))
    }
  } else if (dist == "weibull") { # weibull-weibull survival
    surv <- function(t, x, betas) {
      eta <-  if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      basesurv1 <- exp(-lambdas[1L] * (t ^ gammas[1L]))
      basesurv2 <- exp(-lambdas[2L] * (t ^ gammas[2L]))
      return((pmix * basesurv1 + (1 - pmix) * basesurv2) ^ exp(eta))
    }
  } else if (dist == "exponential" && !mixture) { # exponential survival
    surv <- function(t, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      basesurv <- exp(-lambdas[1L] * t)
      return(basesurv ^ exp(eta))
    }
  } else if (dist == "exponential") { # exponential-exponential survival
    surv <- function(t, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      basesurv1 <- exp(-lambdas[1L] * t)
      basesurv2 <- exp(-lambdas[2L] * t)
      return((pmix * basesurv1 + (1 - pmix) * basesurv2) ^ exp(eta))
    }
  } else if (dist == "gompertz" && !mixture) { # gompertz survival
    surv <- function(t, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      basesurv <- exp(-lambdas[1L]/gammas[1L] * (exp(gammas[1L] * t) - 1))
      return(basesurv ^ exp(eta))
    }
  } else if (dist == "gompertz") { # gompertz-gompertz survival
    surv <- function(t, x, betas) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      basesurv1 <- exp(-lambdas[1L]/gammas[1L] * (exp(gammas[1L] * t) - 1))
      basesurv2 <- exp(-lambdas[2L]/gammas[2L] * (exp(gammas[2L] * t) - 1))
      return((pmix * basesurv1 + (1 - pmix) * basesurv2) ^ exp(eta))
    }
  }
  return(surv)
}

# Return the hazard function for the specified distribution
#
# @param dist The name of the distribution for the baseline hazard
# @param lambdas The scale parameter(s) for the baseline hazard
# @param gammas The shape parameter(s) for the baseline hazard
# @param mixture Logical specifying whether to use a 2-component mixture model
#   for the baseline hazard distribution
# @param pmix Scalar specifying the mixture parameter, must be between 0 and 1
# @param tde Vector of log hazard ratios for the time dependent effects
# @param tdefunction Optionally, a function to use for the interaction with time
#   when creating time dependent effects
get_hazard <- function(dist = c("weibull", "exponential", "gompertz"),
                       lambdas = NULL, gammas = NULL, mixture = FALSE,
                       pmix = 0.5, tde, tdefunction = NULL) {
  validate_lambdas(lambdas = lambdas, dist = dist, mixture = mixture)
  validate_gammas(gammas = gammas, dist = dist, mixture = mixture)
  if (pmix < 0 || pmix > 1)
    stop("'pmix' must be between 0 and 1.", call. = FALSE)
  if (dist == "weibull" && !mixture) { # weibull hazard
    haz <- function(t, x, betas, tde, tdefunction) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      if (!is.null(tde))
        eta <- eta +
          sum(sapply(names(tde), function(i) tde[[i]] * tdefunction(t) * x[[i]]))
      basehaz <- gammas[1L] * lambdas[1L] * (t ^ (gammas[1L] - 1))
      return(basehaz * exp(eta))
    }
  } else if (dist == "weibull") { # weibull-weibull hazard
    haz <- function(t, x, betas, tde, tdefunction) {
      eta <-  if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      if (!is.null(tde))
        eta <- eta +
          sum(sapply(names(tde), function(i) tde[[i]] * tdefunction(t) * x[[i]]))
      basesurv1 <- exp(-lambdas[1L] * (t ^ gammas[1L]))
      basesurv2 <- exp(-lambdas[2L] * (t ^ gammas[2L]))
      basehaz1 <- gammas[1L] * lambdas[1L] * (t ^ (gammas[1L] - 1))
      basehaz2 <- gammas[2L] * lambdas[2L] * (t ^ (gammas[2L] - 1))
      denom <- pmix * basesurv1 + (1 - pmix) * basesurv2
      numer <- pmix * basesurv1 * basehaz1 + (1 - pmix) * basesurv2 * basehaz2
      basehaz <- numer / denom
      return(basehaz * exp(eta))
    }
  } else if (dist == "exponential" && !mixture) { # exponential hazard
    haz <- function(t, x, betas, tde, tdefunction) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      if (!is.null(tde))
        eta <- eta +
          sum(sapply(names(tde), function(i) tde[[i]] * tdefunction(t) * x[[i]]))
      basehaz <- lambdas[1L]
      return(basehaz * exp(eta))
    }
  } else if (dist == "exponential") { # exponential-exponential hazard
    haz <- function(t, x, betas, tde, tdefunction) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      if (!is.null(tde))
        eta <- eta +
          sum(sapply(names(tde), function(i) tde[[i]] * tdefunction(t) * x[[i]]))
      basesurv1 <- exp(-lambdas[1L] * t)
      basesurv2 <- exp(-lambdas[2L] * t)
      denom <- pmix * basesurv1 + (1 - pmix) * basesurv2
      numer <- pmix * basesurv1 * lambdas[1L] + (1 - pmix) * basesurv2 * lambdas[2L]
      basehaz <- numer / denom
      return(basehaz * exp(eta))
    }
  } else if (dist == "gompertz" && !mixture) { # gompertz hazard
    haz <- function(t, x, betas, tde, tdefunction) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      if (!is.null(tde))
        eta <- eta +
          sum(sapply(names(tde), function(i) tde[[i]] * tdefunction(t) * x[[i]]))
      basehaz <- lambdas[1L] * exp(gammas[1L] * t)
      return(basehaz * exp(eta))
    }
  } else if (dist == "gompertz") { # gompertz-gompertz hazard
    stop("Gompertz 2-component mixture with tde not yet implemented.", call. = FALSE)
    haz <- function(t, x, betas, tde, tdefunction) {
      eta <- if (!is.null(betas))
        sum(sapply(names(betas), function(i) betas[[i]] * x[[i]])) else 0L
      if (!is.null(tde))
        eta <- eta +
          sum(sapply(names(tde), function(i) tde[[i]] * tdefunction(t) * x[[i]]))
      return(NULL)
    }
  }
  return(haz)
}

# Ensure the vector of lambdas/gammas is positive and the correct length
#
# @param lambdas,gammas The user input to the lambdas/gammas argument
# @param dist Character string specifying the survival distribution
# @param mixture Logical specifying whether to use a mixture distribution
validate_lambdas <- function(lambdas = NULL, dist, mixture) {
  if (is.null(lambdas)) {
    stop("'lambdas' must be specified.", call. = FALSE)
  } else if (mixture && (!length(lambdas) == 2L)) {
    stop("'lambdas' should be length 2.", call. = FALSE)
  } else if (!mixture && (!length(lambdas) == 1L)) {
    stop("'lambdas' should be length 1.", call. = FALSE)
  }
  if (any(lambdas < 0))
    stop("'lambdas' should be positive.", call. = FALSE)
}
validate_gammas <- function(gammas = NULL, dist, mixture) {
  if (dist == "exponential") { # exponential
    if (!is.null(gammas))
      stop("'gammas' should not be specified for exponential models.", call. = FALSE)
  } else { # weibull or gompertz
    if (is.null(gammas)) {
      stop("'gammas' must be specified for weibull and gompertz models.", call. = FALSE)
    } else if (mixture && (!length(gammas) == 2L)) {
      stop("'gammas' should be length 2.", call. = FALSE)
    } else if (!mixture && (!length(gammas) == 1L)) {
      stop("'gammas' should be length 1.", call. = FALSE)
    }
    if (any(gammas < 0))
      stop("'gammas' should be positive.", call. = FALSE)
  }
}

# Calculate the cumulative hazard at time t using Gauss-Kronrod quadrature,
# and then returning the survival probability (i.e. exp(-cumhaz))
# minus a random uniform.
# Setting the returned value of this function to zero, and solving for t,
# should provide the simulated survival time for one individual.
#
# @param t The event time, unknown but solution to be found using \code{uniroot}
# @param hazard The user-defined hazard function, with named arguments t, x,
#   betas, ...
# @param x Vector of covariate data to be supplied to hazard.
# @param betas Vector of parameter values to be supplied to hazard.
# @param qq The standardised quadpoints and quadweights returned by a call to
#   get_quadpoints.
# @param ... Further arguments passed to hazard.
rootfn_hazard <- function(t, hazard, x = NULL, betas = NULL,
                   u = stats::runif(1), qq = get_quadpoints(nodes = 15), ...) {
  qpts <- unstandardise_quadpoints(qq$points, 0, t)
  qwts <- unstandardise_quadweights(qq$weights, 0, t)
  cumhaz <- sum(unlist(lapply(1:length(qpts), function(q) {
    qwts[[q]] * hazard(t = qpts[[q]], x = x, betas = betas, ...)
  })))
  surv <- exp(-cumhaz)
  return(surv - u)
}

# Function for calculating the survival probability at time t minus a
# random uniform.
# Setting the returned value of this function to zero, and solving for t,
# should provide the simulated survival time for one individual.
#
# @param t The event time, unknown but solution to be found using \code{uniroot}
# @param cumhazard The user-defined cumulative hazard function, with named
#   arguments x, betas, aux
# @param x Vector of covariate data to be supplied to survival.
# @param betas Vector of parameter values to be supplied to survival.
# @param ... Further arguments passed to survival.
rootfn_cumhazard <- function(t, cumhazard, x = NULL, betas = NULL,
                             u = stats::runif(1), ...) {
  cumhaz <- cumhazard(t = t, x = x, betas = betas, ...)
  surv <- exp(-cumhaz)
  return(surv - u)
}

# Function for calculating the survival probability at time t minus a
# random uniform.
# Setting the returned value of this function to zero, and solving for t,
# should provide the simulated survival time for one individual.
#
# @param t The event time, unknown but solution to be found using \code{uniroot}
# @param hazard The survival function, with named arguments x, betas, aux
# @param x Vector of covariate data to be supplied to survival.
# @param betas Vector of parameter values to be supplied to survival.
# @param ... Further arguments passed to survival.
rootfn_surv <- function(t, survival, x = NULL, betas = NULL,
                        u = stats::runif(1), ...) {
  surv <- survival(t = t, x = x, betas = betas, ...)
  return(surv - u)
}

# Check that x is either NULL or a data frame
#
# @param x Object to check
validate_x <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x) && !is.data.frame(x))
    stop("'", nm, "' should be a data frame.")
  x
}

# Check that x is either NULL, a named vector, a data frame, or
# a list of data frames
#
# @param x Object to check
validate_betas <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x) && !is.data.frame(x) && !is(x, "list")) {
    if (!is.vector(x)) {
      stop("'", nm, "' should be a vector, a data frame or a list of data frames.")
    } else if (is.null(names(x))) {
      stop("'", nm, "' should be a named vector.")
    }
  }
  if (!is.null(x) && is(x, "list")) {
    checks <- sapply(x, is.data.frame)
    if (!all(checks))
      stop("'", nm, "' should be a vector, a data frame or a list of data frames.")
  }
  x
}

# Check that x is either NULL, or a named vector
#
# @param x Object to check
validate_tde <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.vector(x) || is.null(names(x))) {
    stop("'", nm, "' should be a named vector.")
  }
  x
}

# Extract rows of data corresponding to row i (if idvar is NULL) or
# individual i (if idvar is not NULL)
#
# @param x A data frame or list of data frames
# @param i The row index or the id value
# @param idvar The name of the ID variable
subset_df <- function(x, i, idvar = NULL) {
  if (is.null(x)) {
    return(x)
  } else if (is.data.frame(x) && is.null(idvar)) {
    return(x[i, , drop = FALSE])
  } else if (is.data.frame(x)) {
    check_for_idvar_and_id(x, idvar = idvar, id = i)
    return(x[x[[idvar]] == i, , drop = FALSE])
  } else if (is(x, "list") && is.null(idvar)) {
    return(lapply(x, function(tmp) tmp[i, , drop=FALSE]))
  } else if (is(x, "list")) {
    return(lapply(x, function(tmp) {
      check_for_idvar_and_id(tmp, idvar = idvar, id = i)
      tmp[tmp[[idvar]] == i, , drop=FALSE]}))
  } else if (is.vector(x)) {
    return(x)
  } else {
    stop("'x' should be NULL, a vector, a data frame, or a list of data frames.")
  }
}

# Check that the individual's ID appears in the data frame
#
# @param df A data frame (being either 'x' or 'betas' in the simsurv call)
# @param idvar The name of the ID variable
# @param id The current individual's ID value
check_for_idvar_and_id <- function(df, idvar, id) {
  if (!idvar %in% colnames(df))
    stop("The variable '", idvar, "' does not appear in all data frames.", call. = FALSE)
  if (!id %in% df[[idvar]])
    stop("The individual '", id, "' does not appear in all data frames.", call. = FALSE)
}

# Consistent error messages
STOP_nan_at_limit <- function() {
  stop("Could not evaluate the survival probability for some individuals at ",
       "the upper limit of 'interval' (likely because the hazard was infinite). ",
       "Try specifying a smaller upper limit using the 'interval' argument.",
       call. = FALSE)
}
STOP_increase_limit <- function() {
  stop("Could not find the simulated survival time for some individuals within ",
       "the specified interval. Consider increasing the upper limit of the ",
       "interval using the 'interval' argument.", call. = FALSE)
}


# Convert a standardised quadrature node to an unstandardised value based on
# the specified integral limits
#
# @param t A vector of standardised quadrature nodes
# @param a The lower limit of the integral
# @param b The upper limit of the integral
unstandardise_quadpoints <- function(t, a, b) {
  ((b - a) / 2) * t + ((b + a) / 2)
}

# Convert a standardised quadrature weight to an unstandardised value based on
# the specified integral limits
#
# @param t A vector of standardised quadrature weights
# @param a The lower limit of the integral
# @param b The upper limit of the integral
unstandardise_quadweights <- function(t, a, b) {
  ((b - a) / 2) * t
}

# Function to return standardised GK quadrature points and weights
#
# @param nodes The required number of quadrature nodes
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
  if (!is.numeric(nodes) || (length(nodes) > 1L)) {
    stop("'nodes' should be a numeric vector of length 1.")
  } else if (nodes == 15) {
    list(
      points = c(
        -0.991455371120812639207,
        -0.949107912342758524526,
        -0.86486442335976907279,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
        0,
        0.2077849550078984676007,
        0.405845151377397166907,
        0.5860872354676911302941,
        0.741531185599394439864,
        0.86486442335976907279,
        0.9491079123427585245262,
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))
  } else if (nodes == 11) {
    list(
      points = c(
        -0.984085360094842464496,
        -0.906179845938663992798,
        -0.754166726570849220441,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.754166726570849220441,
        0.906179845938663992798,
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))
  } else stop("'nodes' must be either 7, 11 or 15.")
}
