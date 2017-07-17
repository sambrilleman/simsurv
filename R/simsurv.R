#' Simulate survival times from a user-specified hazard function
#'
#' Simulate survival times from a user-specified hazard function.The hazard
#' function is allowed to be time-dependent, and the resulting cumulative
#' hazard function does not need to have a closed-form solution.
#' The package is based on the approach described in Crowther and Lambert (2013),
#' whereby the cumulative hazard is evaluated using numerical quadrature and
#' survival times are generated using an iterative algorithm which nests the
#' quadrature-based evaluation of the cumulative hazard inside Brent's (1973)
#' univariate root finder (for the latter
#' the \code{\link[stats]{uniroot}} function is used). Not requiring a
#' closed form solution to the cumulative hazard function has the benefit that
#' survival times can be generated for complex models such as joint
#' longitudinal and survival models; the \strong{Examples} section provides
#' an example of this.
#'
#' @export
#'
#' @param hazfn The user specified hazard function, with named arguments
#'   \code{t}, \code{x} and \code{pars}. See the \strong{Details} section for
#'   a description of these arguments.
#' @param x A data frame, or possibly a list of data frames, containing the
#'   covariates to be supplied to \code{hazfn}. If \code{idvar = NULL} then
#'   each row of the data frame should supply covariate data for one individual.
#' @param pars A data frame, or possibly a list of data frames, containing the
#'   parameter values to be supplied to \code{hazfn}. If \code{idvar = NULL}
#'   then each row of the data frame should supply parameter values for one
#'   individual.
#' @param idvar The name of the ID variable identifying individuals. This is
#'   only required when \code{x} and/or \code{pars} contain multiple rows per
#'   individual. Otherwise, if \code{idvar = NULL} then each row of \code{x}
#'   and \code{pars} is assumed to correspond to a different individual.
#' @param ids A vector containing the unique values of \code{idvar} (i.e. the
#'   unique individual IDs). This is only required when \code{x} and/or
#'   \code{pars} contain multiple rows per individual. Otherwise, if
#'   \code{ids = NULL} then each row of \code{x} and \code{pars} is assumed to
#'   correspond to a different individual.
#' @param maxt The maximum event time. For simulated event times greater than
#'   \code{maxt}, the event time (\code{"eventtime"}) returned in the data frame
#'   will be truncated at \code{maxt} and the event indicator (\code{"status"})
#'   will be set to zero indicating that the individual was right-censored.
#' @param qnodes Integer specifying the number of quadrature nodes to use for
#'   the Gauss-Kronrod quadrature.
#' @param interval The interval over which to search for the
#'   \code{\link[stats]{uniroot}} corresponding to each simulated event time.
#' @param ... Other arguments passed to \code{hazfn}.
#'
#' @return A data frame a row for each individual, and the following three
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
#'   4118–-4134. \url{doi:10.1002/sim.5823}
#'
#'   Bender R, Augustin T, and Blettner M. (2005) Generating survival times to
#'   simulate Cox proportional hazards models. \emph{Statistics in Medicine}
#'   \strong{24}(11), 1713--1723.
#'
#'   Brent R. (1973) \emph{Algorithms for Minimization without Derivatives}.
#'   Englewood Cliffs, NJ: Prentice-Hall.
#'
#' @examples
#'   #######
#'   # Here we present an example of simulating survival times
#'   # based on a joint longitudinal and survival model
#'
#'   # First we define the hazard function to pass to simsurv
#'   # (NB this is a Weibull proportional hazards regression submodel
#'   # from a joint longitudinal and survival model with a "current
#'   # value" association structure).
#'   weibull_ph_hazfn <- function(t, x, pars) {
#'       pars[["shape"]] * (t ^ (pars[["shape"]] - 1)) * exp(
#'         pars[["betaEvent_intercept"]] +
#'         pars[["betaEvent_binary"]] * x[["Z1"]] +
#'         pars[["betaEvent_continuous"]] * x[["Z2"]] +
#'         pars[["betaEvent_assoc"]] * (
#'           pars[["betaLong_intercept"]] +
#'           pars[["betaLong_slope"]] * t +
#'           pars[["betaLong_binary"]] * x[["Z1"]] +
#'           pars[["betaLong_continuous"]] * x[["Z2"]]
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
#'   covs <- data.frame(
#'     Z1 = rbinom(N, 1, 0.45), # a binary covariate
#'     Z2 = rnorm(N, 44, 8.5)   # a continuous covariate
#'   )
#'
#'   # Then we simulate the survival times based on the
#'   # hazard function, covariates, and true parameter values
#'   s1 <- simsurv(hazfn = weibull_ph_hazfn, x = covs, pars = betas, maxt = 10)
#'   head(s1)
#'
simsurv <- function(hazfn, x = NULL, pars = NULL, idvar = NULL, ids = NULL,
                    maxt = NULL, qnodes = 15, interval = c(0, 500),
                    seed = NULL, ...) {
  if (!is.null(seed))
    set.seed(seed)
  if (!is.function(hazfn))
    stop("'hazfn' should be a function")
  ok_args <- c("t", "x", "pars")
  nm_args <- names(formals(hazfn))
  if (!all(ok_args %in% nm_args))
    stop("'hazfn' function should have the following named arguments: ",
         paste(ok_args, collapse = ", "))
  x <- validate_df(x)
  pars <- validate_df(pars)
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
  tt <- sapply(ids, function(i) {
    x_i <- subset_df(x, i, idvar = idvar)
    pars_i <- subset_df(pars, i, idvar = idvar)
    u_i <- runif(1)
    # check whether S(t) is still greater than random uniform variable u_i at the
    # upper limit of uniroot's interval (otherwise uniroot will return an error)
    at_limit <- rootfn(interval[2], hazfn = hazfn, x = x_i, pars = pars_i,
                       u = u_i, qnodes = qnodes, ...)
    t_i <- if (at_limit > 0) interval[2] else
      stats::uniroot(rootfn, hazfn = hazfn, x = x_i, pars = pars_i,
                     u = u_i, qnodes = qnodes, ..., interval = interval)$root
    return(t_i)
  })
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
  return(ret)
}


#----------- internal

# Function for calculating the cumulative hazard at time t using Gauss-Kronrod
# quadrature, and then returning the survival probability (i.e. exp(-cumhaz))
# minus a random uniform.
# Setting the returned value of this function to zero, and solving for t,
# should provide the simulated survival time for one individual.
#
# @param t The event time, unknown but solution to be found using \code{uniroot}
# @param hazfn The user specified hazard function, with named arguments x,
#   beta, aux
# @param x Vector of covariate data to be supplied to hazfn
# @param pars Vector of parameter values to be supplied to hazfn
# @param qnodes Integer specifying the number of quadrature nodes
# @param ... Further arguments passed to hazfn.
rootfn <- function(t, hazfn, x = NULL, pars = NULL,
                   u = runif(1), qnodes = 15, ...) {
  qq <- get_quadpoints(qnodes)
  qpts <- lapply(qq$points,  unstandardise_quadpoints,  0, t)
  qwts <- lapply(qq$weights, unstandardise_quadweights, 0, t)
  cumhaz <- sum(sapply(seq(qnodes), function(q) {
    qwts[[q]] * hazfn(t = qpts[[q]], x = x, pars = pars, ...)
  }))
  return(exp(-cumhaz) - u)
}

# Check that x is either NULL, a data frame, or a list of data frames
#
# @param x Object to check
validate_df <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x) && !is.data.frame(x) && !is(x, "list"))
    stop("'", nm, "' should be a data frame or a list of data frames.")
  if (!is.null(x) && is(x, "list")) {
    checks <- sapply(x, is.data.frame)
    if (!all(checks))
      stop("'", nm, "' should be a data frame or a list of data frames.")
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
  } else {
    stop("'x' should be NULL, a data frame, or a list of data frames.")
  }
}

# Check that the individual's ID appears in the data frame
#
# @param df A data frame (being either 'x' or 'pars' in the simsurv call)
# @param idvar The name of the ID variable
# @param id The current individual's ID value
check_for_idvar_and_id <- function(df, idvar, id) {
  if (!idvar %in% colnames(df))
    stop("The variable '", idvar, "' does not appear in all data frames.", call. = FALSE)
  if (!id %in% df[[idvar]])
    stop("The individual '", id, "' does not appear in all data frames.", call. = FALSE)
}

# Convert a standardised quadrature node to an unstandardised value based on
# the specified integral limits
#
# @param t An unstandardised quadrature node
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_quadpoints <- function(t, a, b) {
  if (!identical(length(t), 1L) || !is.numeric(t))
    stop("'t' should be a single numeric value.", call. = FALSE)
  if (!all(is.numeric(a), is.numeric(b)))
    stop("'a' and 'b' should be numeric.", call. = FALSE)
  if (!length(a) %in% c(1L, length(b)))
    stop("'a' and 'b' should be vectors of length 1, or, be the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
  ((b - a) / 2) * t + ((b + a) / 2)
}

# Convert a standardised quadrature weight to an unstandardised value based on
# the specified integral limits
#
# @param t An unstandardised quadrature weight
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_quadweights <- function(t, a, b) {
  if (!identical(length(t), 1L) || !is.numeric(t))
    stop("'t' should be a single numeric value.", call. = FALSE)
  if (!all(is.numeric(a), is.numeric(b)))
    stop("'a' and 'b' should be numeric.", call. = FALSE)
  if (!length(a) %in% c(1L, length(b)))
    stop("'a' and 'b' should be vectors of length 1, or, be the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
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
    stop("'quadnodes' should be a numeric vector of length 1.")
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
  } else stop("'quadnodes' must be either 7, 11 or 15.")
}
