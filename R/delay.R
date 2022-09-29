# mkuhn, 2021-04-06
# delay distribution functions


#' Delayed Exponential Distribution
#'
#' @description
#' Density, distribution function, quantile function, random generation and restricted mean survival time function for the delayed exponential distribution.
#' There is an initial delay phase (parameter `delay1`) where no events occur. After that, `rate1` applies.
#' Optionally, a second phase is possible where the hazard rate might change (parameters `delay2` and `rate2`).
#'
#' @details
#' Additional arguments are forwarded via `...` to the underlying functions of the exponential distribution in the `stats`-package.
#' If only a single initial delay phase is there, the numerical arguments other than `n` are recycled to the length of the result (as with the exponential distribution in `stats`).
#' With two phases, the arguments are **not** recycled. Only the first element of delays and rates are used as it otherwise becomes ambiguous which delay and rate parameter apply for observations in different phases.
#' Generally, only the first elements of the logical arguments are used.
#'
#' @param x A numeric vector of values for which to get the density.
#' @param q A numeric vector of quantile values.
#' @param t A numeric vector of times that restrict the mean survival. Default is `+Inf`, i.e., the unrestricted mean survival time.
#' @param p A numeric vector of probabilities.
#' @param n integer. Number of random observations requested.
#' @param delay1 numeric. The first delay, must be non-negative.
#' @param rate1 numeric. The event rate, must be non-negative.
#' @param delay numeric. Alias for first delay.
#' @param rate numeric. Alias for first rate.
#' @param delay2 numeric. The second delay, must be non-negative.
#' @param rate2 numeric. The second event rate, must be non-negative.
#' @param ... further arguments are passed on to the underlying non-delayed function, e.g., `lower.tail=` to [stats::pexp()]
#' @return Functions pertaining to the delayed exponential distribution:
#' * `dexp_delayed` gives the density
#' * `pexp_delayed` gives the distribution function
#' * `qexp_delayed` gives the quantile function
#' * `rexp_delayed` generates a pseudo-random sample
#' * `mexp_delayed` gives the restricted mean survival time
#'
#' The length of the result is determined by `n` for `rexp_delayed`, and is the maximum of the lengths of the numerical arguments for the other functions,
#' R's recycling rules apply when only single initial delay phase is used.
#' @seealso stats::Exponential
#' @keywords distribution
#' @name DelayedExponential
NULL

#' @rdname DelayedExponential
#' @export
dexp_delayed <- function(x, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL, log = FALSE) {
  stopifnot( length(log) >= 1L, is.logical(log) )
  log <- log[[1L]] # only first value of log is used
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(rate1)) )

  # check for easy case: only a single phase
  if ( is.null(delay2) ) {
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(stats::dexp(x = x - delay1, rate = rate1, log = log))
  }

  # two phases
  # take only first value of parameter arguments!
  if ( length(delay1) > 1L || length(rate1) > 1L || length(delay2) > 1L || length(rate2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }
  # first phase
  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  dvals <- stats::dexp(x = x - delay1, rate = rate1, log = log)
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # second phase
  # update density for observations that lie in the 2nd phase
  phase2Ind <- which(x >= delay2)
  if (length(phase2Ind)) {
    dvals[phase2Ind] <- stats::dexp(x = x[phase2Ind] - delay2, rate = rate2, log = log)
    dvals[phase2Ind] <- if (log) dvals[phase2Ind] - rate1 * (delay2 - delay1) else dvals[phase2Ind] * exp(-rate1 * (delay2 - delay1))
  }

  dvals
}

#' @rdname DelayedExponential
#' @export
pexp_delayed <- function(q, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(rate1)) )

  # check for easy case: only a single delay
  if ( is.null(delay2) ) {
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(stats::pexp(q = q - delay1, rate = rate1, ...))
  }

  # two phases
  if ( length(delay1) > 1L || length(rate1) > 1L || length(delay2) > 1L || length(rate2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }
  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # first phase
  pvals <- stats::pexp(q = q - delay1, rate = rate1, ...)
  # check if we have observations in 2nd phase
  phase2Ind <- which(q > delay2)
  if (length(phase2Ind)) pvals[phase2Ind] <- stats::pexp(q = q[phase2Ind] - delay2 + rate1/rate2 * (delay2 - delay1), rate = rate2, ...)

  pvals
}

#' @rdname DelayedExponential
#' @export
qexp_delayed <- function(p, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL, lower.tail = TRUE, log.p = FALSE) {
  #lower.tail = TRUE, log.p = FALSE
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(rate1)) )

  # check for easy case: only a single delay phase
  if ( is.null(delay2) ){
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(delay1 + stats::qexp(p = p, rate = rate1, lower.tail = lower.tail, log.p = log.p)) #lower.tail = lower.tail, log.p = log.p
  }

  # two phases
  if ( length(delay1) > 1L || length(rate1) > 1L || length(delay2) > 1L || length(rate2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }

  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  # we need both delay2 AND rate2
  stopifnot( is.finite(delay2), is.finite(rate2) )
  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # first phase
  qvals <- delay1 + stats::qexp(p = p, rate = rate1, lower.tail = lower.tail, log.p = log.p)

  # second phase
  # check if we have observations in 2nd phase
  # transform p-values to canonical meaning (lower.tail=T, log.p=F)
  if (!lower.tail) p <- 1L - p
  if (log.p) p <- exp(p)
  phase2Ind <- which(p >= pexp_delayed(q = delay2, delay1 = delay1, rate1 = rate1, lower.tail = TRUE, log.p = FALSE))
  if (length(phase2Ind)) qvals[phase2Ind] <- delay2 - rate1/rate2 * (delay2 - delay1) - log(1-p[phase2Ind])/rate2 #stats::qexp(p = p[phase2Ind], rate = rate2, lower.tail = TRUE, log.p = FALSE)

  qvals
}

#' @rdname DelayedExponential
#' @export
rexp_delayed <- function(n, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(rate1)) )

  # single phase
  # check for easy case: only a single delay
  if ( is.null(delay2) ){
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(delay1 + stats::rexp(n = n, rate = rate1))
  }

  # two phases
  if ( length(delay1) > 1L || length(rate1) > 1L || length(delay2) > 1L || length(rate2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }
  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # use inverse CDF-method
  qexp_delayed(p = stats::runif(n = n), delay1 = delay1, rate1 = rate1, delay2 = delay2, rate2 = rate2)
}


#' @rdname DelayedExponential
#' @export
mexp_delayed <- function(t=+Inf, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(rate1)) )

  # single phase
  # calculate for single phase delayed exponential
  if ( is.null(delay2) ){
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(pmin.int(t, delay1) + pexp_delayed(q = t, delay1 = delay1, rate1 = rate1) / rate1)
  }

  # two phases
  if ( length(delay1) > 1L || length(rate1) > 1L || length(delay2) > 1L || length(rate2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }
  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  pmin.int(t, delay1) + pexp_delayed(q = pmin.int(t, delay2), delay1 = delay1, rate1 = rate1) / rate1 +
    pexp_delayed(q=delay2, delay1 = delay1, rate1 = rate1, lower.tail = FALSE) * pexp_delayed(q=t, delay1 = delay2, rate1 = rate2) / rate2
}



#' Delayed Weibull Distribution
#'
#' @description
#' Density, distribution function, quantile function and random generation for the delayed Weibull distribution.
#' Besides the additional parameter `delay`, the other two Weibull-parameters are in principle retained as in R's stats-package:
#' * `shape`
#' * `scale` (as inverse of rate)
#'
#' @details
#' Additional arguments are forwarded via `...` to the underlying functions of the exponential distribution in the stats-package.
#'
#' The numerical arguments other than `n` are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @param x A numeric vector of values for which to get the density.
#' @param q A numeric vector of quantile values.
#' @param t A numeric vector of times that restrict the mean survival. Default is `+Inf`, i.e., the unrestricted mean survival time.
#' @param p A numeric vector of probabilities.
#' @param n integer. Number of random observations requested.
#' @param delay1 numeric. The first delay, must be non-negative.
#' @param shape1 numeric. First shape parameter, must be positive.
#' @param scale1 numeric. First scale parameter (inverse of rate), must be positive.
#' @param delay numeric. Alias for first delay.
#' @param shape numeric. Alias for first shape.
#' @param scale numeric. Alias for first scale.
#' @param ... further arguments are passed on to the underlying non-delayed function, e.g., [stats::dweibull()]
#' @return Functions pertaining to the delayed Weibull distribution:
#' * `dweib_delayed` gives the density
#' * `pweib_delayed` gives the distribution function
#' * `qweib_delayed` gives the quantile function
#' * `rweib_delayed` generates a pseudo-random sample
#' * `mweib_delayed` gives the restricted mean survival time
#'
#' The length of the result is determined by `n` for `rweib_delayed`, and is the maximum of the lengths of the numerical arguments for the other functions, R's recycling rules apply.
#' @keywords distribution
#' @name DelayedWeibull
NULL

#' @rdname DelayedWeibull
#' @export
dweib_delayed <- function(x, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1,
                          delay2 = NULL, shape2 = NULL, scale2 = 1, log = FALSE) {
  stopifnot( length(log) >= 1L, is.logical(log) )
  log <- log[[1L]] # only first value of log is used

  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(shape1), is.finite(scale1)) )


  # check for easy case: only a single phase
  if ( is.null(delay2) ) {
    if (!is.null(shape2) || ! missing(scale2)) warning("Arguments shape2= and/or scale2= are ignored, as argument delay2= is not set.", call. = FALSE)
    return(stats::dweibull(x = x - delay1, shape = shape1, scale = scale1, log = log))
  }


  # two phases
  # take only first value of parameter arguments!
  if ( length(delay1) > 1L || length(shape1) > 1L || length(scale1) > 1L || length(delay2) > 1L || length(shape2) > 1L || length(scale2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }

  delay1 <- delay1[[1L]]
  shape1 <- shape1[[1L]]
  scale1 <- scale1[[1L]]
  # we need both delay2 AND shape2 AND scale2
  delay2 <- delay2[[1L]]
  shape2 <- shape2[[1L]]
  scale2 <- scale2[[1L]]
  stopifnot( is.finite(delay2), is.finite(shape2), is.finite(scale2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # first phase densities
  dvals <- stats::dweibull(x = x - delay1, shape = shape1, scale = scale1, log = log)
  # second phase
  # update density for observations that lie in the 2nd phase
  phase2Ind <- which(x >= delay2)
  if (length(phase2Ind)) {
    dvals[phase2Ind] <- stats::dweibull(x = x[phase2Ind] - delay2, shape = shape2, scale = scale2, log = log)
    dvals[phase2Ind] <- if (log) dvals[phase2Ind] - ((delay2 - delay1)/scale1)^shape1 else dvals[phase2Ind] * exp(-((delay2 - delay1)/scale1)^shape1)
  }

  dvals
}

#' @rdname DelayedWeibull
#' @export
pweib_delayed <- function(q, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1,
                          delay2 = NULL, shape2 = NULL, scale2 = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(shape1), is.finite(scale1)) )

  # check for easy case: only a single delay
  if ( is.null(delay2) ) {
    if (!is.null(shape2) || ! missing(scale2)) warning("Arguments shape2= and/or scale2= are ignored, as argument delay2= is not set.", call. = FALSE)
    return(stats::pweibull(q = q - delay1, shape = shape1, scale = scale1, lower.tail = lower.tail, log.p = log.p))
  }

  # two phases
  if ( length(delay1) > 1L || length(shape1) > 1L || length(scale1) > 1L || length(delay2) > 1L || length(shape2) > 1L || length(scale2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }

  delay1 <- delay1[[1L]]
  shape1 <- shape1[[1L]]
  scale1 <- scale1[[1L]]
  # we need both delay2 AND shape2 AND scale2
  delay2 <- delay2[[1L]]
  shape2 <- shape2[[1L]]
  scale2 <- scale2[[1L]]
  stopifnot( is.finite(delay2), is.finite(shape2), is.finite(scale2) )

  # first phase
  pvals <- stats::pweibull(q = q - delay1, shape = shape1, scale = scale1, lower.tail = lower.tail, log.p = log.p)

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # check if we have observations in 2nd phase
  phase2Ind <- which(q > delay2)
  if (length(phase2Ind)) {
    # probability values to exceed time q
    pvals[phase2Ind] <- exp(-((delay2 - delay1)/scale1)^shape1 - ((q[phase2Ind] - delay2)/scale2)^shape2)
    if (lower.tail) pvals[phase2Ind] <- 1L - pvals[phase2Ind]
    if (log.p) pvals[phase2Ind] <- log(pvals[phase2Ind])
  }

  pvals
}

#' @rdname DelayedWeibull
#' @export
qweib_delayed <- function(p, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1,
                          delay2 = NULL, shape2 = NULL, scale2 = 1, lower.tail = TRUE, log.p = FALSE) {
  stopifnot( is.logical(lower.tail), length(lower.tail) >= 1L, is.logical(log.p), length(log.p) >= 1L)
  lower.tail <- isTRUE(lower.tail[[1L]])
  log.p <- isTRUE(log.p[[1L]])

  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(shape1), is.finite(scale1)) )

  # check for easy case: only a single delay phase
  if ( is.null(delay2) ){
    if (!is.null(shape2) || ! missing(scale2)) warning("Arguments shape2= and/or scale2= are ignored, as argument delay2= is not set.", call. = FALSE)
    return(delay1 + stats::qweibull(p = p, shape = shape1, scale = scale1, lower.tail = lower.tail, log.p = log.p))
  }

  # two phases
  if ( length(delay1) > 1L || length(shape1) > 1L || length(scale1) > 1L || length(delay2) > 1L || length(shape2) > 1L || length(scale2) > 1L ){
    warning("In two-phase setting we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }

  delay1 <- delay1[[1L]]
  shape1 <- shape1[[1L]]
  scale1 <- scale1[[1L]]
  # we need both delay2 AND shape2 AND scale2
  delay2 <- delay2[[1L]]
  shape2 <- shape2[[1L]]
  scale2 <- scale2[[1L]]
  stopifnot( is.finite(delay2), is.finite(shape2), is.finite(scale2) )

  # first phase
  qvals <- delay1 + stats::qweibull(p = p, shape = shape1, scale = scale1, lower.tail = lower.tail, log.p = log.p)

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # second phase
  # check if we have observations in 2nd phase
  # transform p-values to canonical meaning (lower.tail=T, log.p=F)
  if (!lower.tail) p <- 1L - p
  if (log.p) p <- exp(p)
  phase2Ind <- which(p >= pweib_delayed(q = delay2, delay1 = delay1, shape1 = shape1, scale1 = scale1, lower.tail = TRUE, log.p = FALSE))
  if (length(phase2Ind)) qvals[phase2Ind] <- delay2 + scale2 * (-log(1L-p[phase2Ind]) - ((delay2 - delay1)/scale1)^shape1)^(1/shape2)

  qvals
}

#' @rdname DelayedWeibull
#' @export
rweib_delayed <- function(n, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1,
                          delay2 = NULL, shape2 = NULL, scale2 = 1){
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stopifnot( all(is.finite(delay1), is.finite(shape1), is.finite(scale1)) )

  # single phase
  # check for easy case: only a single delay phase
  if ( is.null(delay2) ){
    if (!is.null(shape2) || ! missing(scale2)) warning("Arguments shape2= and/or scale2= are ignored, as argument delay2= is not set.", call. = FALSE)
    return(delay1 + stats::rweibull(n = n, shape = shape1, scale = scale1))
  }

  # two phases
  if ( length(delay1) > 1L || length(shape1) > 1L || length(scale1) > 1L || length(delay2) > 1L || length(shape2) > 1L || length(scale2) > 1L ){
    warning("In two-phase setting, we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
  }

  delay1 <- delay1[[1L]]
  shape1 <- shape1[[1L]]
  scale1 <- scale1[[1L]]
  # we need both delay2 AND shape2 AND scale2
  delay2 <- delay2[[1L]]
  shape2 <- shape2[[1L]]
  scale2 <- scale2[[1L]]
  stopifnot( is.finite(delay2), is.finite(shape2), is.finite(scale2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # use inverse CDF-method
  qweib_delayed(p = stats::runif(n = n, min = 0L, max = 1L),
                delay1 = delay1, shape1 = shape1, scale1 = scale1,
                delay2 = delay2, shape2 = shape2, scale2 = scale2)
}

#' @rdname DelayedWeibull
#' @export
mweib_delayed <- function(t=+Inf, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1,
                          delay2 = NULL, shape2 = NULL, scale2 = 1) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stopifnot( is.numeric(t) )
  stopifnot( all(is.finite(delay1), is.finite(shape1), is.finite(scale1)) )

  # prepare return value
  mvals <- t
  is.na(mvals) <- is.na(t) # propagate NAs

  # check for easy case: only a single delay
  if ( is.null(delay2) ){
    # single phase
    if (!is.null(shape2) || ! missing(scale2)) warning("Arguments shape2= and/or scale2= are ignored, as argument delay2= is not set.", call. = FALSE)

    afterInd <- which(t > delay1)
    # make use of lower incomplete gamma function which is calculated as gamma * pgamma
    mvals[afterInd] <- delay1 + scale1 / shape1 * gamma(1/shape1) * stats::pgamma(q = ((t[afterInd]-delay1)/scale1)^shape1, shape = 1/shape1)

  } else {
    # two phases
    if ( length(delay1) > 1L || length(shape1) > 1L || length(scale1) > 1L || length(delay2) > 1L || length(shape2) > 1L || length(scale2) > 1L ){
      warning("In two-phase setting, we do not recycle parameters. Only the 1st value of the parameter arguments is used!", call. = FALSE)
    }

    delay1 <- delay1[[1L]]
    shape1 <- shape1[[1L]]
    scale1 <- scale1[[1L]]
    # we need both delay2 AND shape2 AND scale2
    delay2 <- delay2[[1L]]
    shape2 <- shape2[[1L]]
    scale2 <- scale2[[1L]]
    stopifnot( is.finite(delay2), is.finite(shape2), is.finite(scale2) )

    # check delay constraint
    if (delay1 >= delay2) {
      stop("First delay phase must antedate the second delay phase!", call. = FALSE)
    }

    phase1Ind <- which(delay1 < t & t <= delay2)
    phase2Ind <- which(t > delay2)
    if (length(phase1Ind)) {
      mvals[phase1Ind] <- delay1 + scale1 / shape1 * gamma(1/shape1) * stats::pgamma(q = ((t[phase1Ind]-delay1)/scale1)^shape1, shape = 1/shape1)
    }
    if (length(phase2Ind)){
      mvals[phase2Ind] <- delay1 + scale1 / shape1 * gamma(1/shape1) * stats::pgamma(q = ((delay2-delay1)/scale1)^shape1, shape = 1/shape1) +
        pweib_delayed(q = delay2, delay1 = delay1, shape1 = shape1, scale1 = scale1, lower.tail = FALSE) * scale2 / shape2 * gamma(1/shape2) * stats::pgamma(q = ((t[phase2Ind]-delay2)/scale2)^shape2, shape = 1/shape2)
    }
  }# 2phase

  mvals
}


#' Get delay distribution function
#' @param distribution character(1). delay distribution.
#' @param type character(1). type of function, cdf: cumulative distribution function, density or random function
#' @param twoPhase logical(1). For `type='param'`, do we model two phases?
#' @param twoGroup logical(1). For type='param', do we have two groups?
#' @param bind character. For type='param', names of parameters that are bind between the two groups.
#' @param transformed logical(1). For type='param', do we need parameter names transformed (as used inside the optimization function?)
#' @return selected distribution function or parameter names
#' @include delay_estimation.R
getDist <- function(distribution = c("exponential", "weibull"), type = c("cdf", "prob", "density", "random", "param"),
                    twoPhase = FALSE, twoGroup = FALSE, bind = NULL, transformed = FALSE) {
  distribution <- match.arg(distribution)
  type <- match.arg(type)


  switch(type,
         # cumulative distribution function
         prob = ,
         cdf     = c(pexp_delayed, pweib_delayed),
         # density function
         density = c(dexp_delayed, dweib_delayed),
         # random draw function
         random  = c(rexp_delayed, rweib_delayed),
         param   = {

           par_list <- list(exponential = c("delay1", "rate1", "delay2", "rate2")[seq_len(2L*(1L + twoPhase))],
                            weibull = c("delay1", "shape1", "scale1", "delay2", "shape2", "scale2")[seq_len(3L*(1L + twoPhase))])

           if (transformed) {
             par_list <- purrr::map(par_list, ~ paste0(.x, "_tr"))
             if (! is.null(bind) && any(nzchar(bind))) bind <- paste0(bind, "_tr")
           }

           if (twoGroup) {
             # bind only parameters from chosen distribution
             myPars <- par_list[[1L + (distribution == 'weibull')]]
             bind <- intersect(myPars, bind) #intersect: enforces original order from myPars
             par_gr <- purrr::map(par_list, ~ setdiff(.x, bind))
             # bind parameters first
             purrr::map(par_gr, ~ c(bind,
                                    paste(rep.int(.x, times = 2L), rep(c("x", "y"), each = length(.x)),
                                          sep = ".")))
           } else par_list

         },
         stop(glue("Unknown attribute of distribution {distribution}."), call. = FALSE)
  )[[1L + (distribution == 'weibull')]]
}

