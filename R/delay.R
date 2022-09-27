# mkuhn, 2021-04-06
# delay distribution functions


#' Delayed Exponential Distribution
#'
#' @description
#' Density, distribution function, quantile function, random generation and restricted mean survival time function for the delayed exponential distribution.
#' There is an initial delay phase (parameter `delay1`) where no events occur.
#' Optionally, a second phase is possible where the hazard rate might change (parameters `delay2` and `rate2`).
#'
#' @details
#' Additional arguments are forwarded via `...` to the underlying functions of the exponential distribution in the `stats`-package.
#' Unlike the distribution functions from `stats` the arguments are **not** recycled. The delay and rate related parameters must have length 1 as it otherwise becomes ambiguous which delay and rate parameter to apply, depending on the observations' phase.
#' Only the first elements of the logical arguments are used.
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
#' The length of the result is determined by `n` for `rexp_delayed`, and is the maximum of the lengths of the numerical arguments for the other functions, R's recycling rules apply.
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

  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  stopifnot( is.numeric(delay1), is.numeric(rate1), is.finite(delay1), is.finite(rate1) )

  # first phase
  dvals <- stats::dexp(x = x - delay1, rate = rate1, log = log)

  # check for easy case: only a single phase
  if ( is.null(delay2) ) {
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(dvals)
  }

  # two phases
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.numeric(delay2), is.numeric(rate2), is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # second phase
  # update density for observations that lie in the 2nd phase
  phase2Ind <- which(x >= delay2)
  if (length(phase2Ind)) {
    dvals[phase2Ind] <- stats::dexp(x = x - delay2, rate = rate2, log = log)[phase2Ind]
    dvals[phase2Ind] <- if (log) dvals[phase2Ind] - rate1 * (delay2 - delay1) else dvals[phase2Ind] * exp(-rate1 * (delay2 - delay1))
  }

  #if (log) dvals - log(K) else dvals/K
  dvals

}

#' @rdname DelayedExponential
#' @export
pexp_delayed <- function(q, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  stopifnot( is.numeric(delay1), is.numeric(rate1), is.finite(delay1), is.finite(rate1) )

  # first phase
  pvals <- stats::pexp(q = q - delay1, rate = rate1, ...)

  # check for easy case: only a single delay
  if ( is.null(delay2) ) {
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(pvals)
  }

  # second phase
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.numeric(delay2), is.numeric(rate2), is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # check if we have observations in 2nd phase
  phase2Ind <- which(q >= delay2)
  if (length(phase2Ind)) pvals[phase2Ind] <- stats::pexp(q = q - delay2 + rate1/rate2 * (delay2 - delay1), rate = rate2, ...)[phase2Ind]

  pvals
}

#' @rdname DelayedExponential
#' @export
qexp_delayed <- function(p, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  stopifnot( is.numeric(delay1), is.numeric(rate1), is.finite(delay1), is.finite(rate1) )

  # first phase
  qvals <- delay1 + stats::qexp(p = p, rate = rate1, ...)

  # check for easy case: only a single delay phase
  if ( is.null(delay2) ){
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(qvals)
  }

  # second phase
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.numeric(delay2), is.numeric(rate2), is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  # check if we have observations in 2nd phase
  phase2Ind <- which(p >= pexp_delayed(q = delay2, delay1 = delay1, rate1 = rate1))
  if (length(phase2Ind)) qvals[phase2Ind] <- delay2 - rate1/rate2 * (delay2 - delay1) + stats::qexp(p = p, rate = rate2, ...)[phase2Ind]

  qvals

}

#' @rdname DelayedExponential
#' @export
rexp_delayed <- function(n, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  stopifnot( is.numeric(delay1), is.numeric(rate1), is.finite(delay1), is.finite(rate1) )

  # single phase
  # check for easy case: only a single delay
  if ( is.null(delay2) ){
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(delay1 + stats::rexp(n = n, rate = rate1))
  }

  # two phases
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.numeric(delay2), is.numeric(rate2), is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must always antedate the second delay phase!", call. = FALSE)
  }

  # use inverse CDF-method
  qexp_delayed(p = stats::runif(n = n), delay1 = delay1, rate1 = rate1, delay2 = delay2, rate2 = rate2)

}


#' @rdname DelayedExponential
#' @export
mexp_delayed <- function(t=+Inf, delay1 = 0, rate1 = 1, delay = delay1, rate = rate1, delay2 = NULL, rate2 = NULL) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(rate)) if (missing(rate1)) rate1 <- rate else warning("Argument rate= is ignored as rate1= is given!", call. = FALSE)

  delay1 <- delay1[[1L]]
  rate1 <- rate1[[1L]]
  stopifnot( is.numeric(delay1), is.numeric(rate1), is.finite(delay1), is.finite(rate1) )

  # single phase
  # calculate for single phase delayed exponential
  if ( is.null(delay2) ){
    if (!is.null(rate2)) warning("Argument rate2= is ignored, as argument delay2= is not set.", call. = FALSE)
    return(pmin.int(t, delay1) + pexp_delayed(q = t, delay1 = delay1, rate1 = rate1) / rate1)
  }

  # two phases
  # we need both delay2 AND rate2
  delay2 <- delay2[[1L]]
  rate2 <- rate2[[1L]]
  stopifnot( is.numeric(delay2), is.numeric(rate2), is.finite(delay2), is.finite(rate2) )

  # check delay constraint
  if (delay1 >= delay2) {
    stop("First delay phase must antedate the second delay phase!", call. = FALSE)
  }

  pmin.int(t, delay1) + pexp_delayed(q = pmin.int(t, delay2), delay1 = delay1, rate1 = rate1) / rate1 +
    (1L-pexp_delayed(q=delay2, delay1 = delay1, rate1 = rate1)) * pexp_delayed(q=t, delay1 = delay2, rate1 = rate2) / rate2
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
dweib_delayed <- function(x, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stats::dweibull(x = x - delay1, shape = shape1, scale = scale1, ...)
}

#' @rdname DelayedWeibull
#' @export
pweib_delayed <- function(q, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  stats::pweibull(q = q - delay1, shape = shape1, scale = scale1, ...)
}

#' @rdname DelayedWeibull
#' @export
qweib_delayed <- function(p, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  delay1 + stats::qweibull(p = p, shape = shape1, scale = scale1, ...)
}

#' @rdname DelayedWeibull
#' @export
rweib_delayed <- function(n, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1){
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  delay1 + stats::rweibull(n = n, shape = shape1, scale = scale1)
}

#' @rdname DelayedWeibull
#' @export
mweib_delayed <- function(t=+Inf, delay1, shape1, scale1 = 1, delay = delay1, shape = shape1, scale = scale1, ...) {
  if (!missing(delay)) if (missing(delay1)) delay1 <- delay else warning("Argument delay= is ignored as delay1= is given!", call. = FALSE)
  if (!missing(shape)) if (missing(shape1)) shape1 <- shape else warning("Argument shape= is ignored as shape1= is given!", call. = FALSE)
  if (!missing(scale)) if (missing(scale1)) scale1 <- scale else warning("Argument scale= is ignored as scale1= is given!", call. = FALSE)

  ifelse(test = t <= delay1, yes = t,
         # make use of lower incomplete gamma function which is calculated as gamma * pgamma
         no = delay1 + scale1 / shape1 * gamma(1/shape1) * stats::pgamma(q = ((t-delay1)/scale1)^shape1, shape = 1/shape1))
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

