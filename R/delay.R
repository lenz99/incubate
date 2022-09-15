# mkuhn, 2021-04-06
# delay distribution functions


#' Delayed Exponential Distribution
#'
#' @description
#' Density, distribution function, quantile function, random generation and restricted mean survival time function for the delayed exponential distribution with `rate`-parameter.
#'
#' @details
#' Additional arguments are forwarded via `...` to the underlying functions of the exponential distribution in the `stats`-package.
#' The numerical arguments other than `n` are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @param x A numeric vector of values for which to get the density.
#' @param q A numeric vector of quantile values.
#' @param t A numeric vector of times that restrict the mean survival. Default is `+Inf`, i.e., the unrestricted mean survival time.
#' @param p A numeric vector of probabilities.
#' @param n integer. Number of random observations requested.
#' @param delay1 numeric. The first delay, must be non-negative.
#' @param rate1 numeric. The event rate, must be non-negative.
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
dexp_delayed <- function(x, delay1 = 0, rate1 = 1, delay2 = NULL, rate2 = NULL, log = FALSE) {
  stopifnot( is.numeric(delay1), is.numeric(rate1), all(is.finite(delay1)), all(is.finite(rate1)) )
  log <- isTRUE(log)

  dvals <- stats::dexp(x = x - delay1, rate = rate1, log = log)
  # check for easy case: only a single delay
  if ( is.null(delay2) ) return(dvals)

  # we need both delay2 AND rate2
  stopifnot( is.numeric(delay2), is.numeric(rate2) )

  # manually recycle delay and rate parameters to ensure they have common length, respectively
  ndelay <- max(length(delay1), length(delay2), na.rm = TRUE)
  nrate <- max(length(rate1), length(rate2), na.rm = TRUE)

  delay1 <- rep_len(delay1, length.out = ndelay); delay2 <- rep_len(delay2, length.out = ndelay)
  rate1  <- rep_len(rate1, length.out = nrate); rate2 <- rep_len(rate2, length.out = nrate)

  # check delay constraint
  if (any(delay1 >= delay2)) {
    stop("First delay phase must always antedate the second delay phase!")
  }

  # check if we have observations in 2nd phase
  phase2Ind <- which(x >= delay2)
  if (length(phase2Ind)) dvals[phase2Ind] <- stats::dexp(x = x - delay2, rate = rate2, log = log)[phase2Ind] * exp(-rate1 * (delay2 - delay1))

  #if (log) dvals - log(K) else dvals/K
  dvals

}

#' @rdname DelayedExponential
#' @export
pexp_delayed <- function(q, delay1 = 0, rate1 = 1, delay2 = NULL, rate2 = NULL, ...) {
  stopifnot( is.numeric(delay1), is.numeric(rate1), all(is.finite(delay1)), all(is.finite(rate1)) )

  pvals <- stats::pexp(q = q - delay1, rate = rate1, ...)

  # check for easy case: only a single delay
  if ( is.null(delay2) ) return(pvals)

  # we need both delay2 AND rate2
  stopifnot( is.numeric(delay2), is.numeric(rate2) )

  # manually recycle delay and rate parameters to ensure they have common length, respectively
  ndelay <- max(length(delay1), length(delay2), na.rm = TRUE)
  nrate <- max(length(rate1), length(rate2), na.rm = TRUE)

  delay1 <- rep_len(delay1, length.out = ndelay); delay2 <- rep_len(delay2, length.out = ndelay)
  rate1 <- rep_len(rate1, length.out = nrate); rate2 <- rep_len(rate2, length.out = nrate)

  # check delay constraint
  if (any(delay1 >= delay2)) {
    stop("First delay phase must always antedate the second delay phase!")
  }

  # check if we have observations in 2nd phase
  phase2Ind <- which(q >= delay2)
  if (length(phase2Ind)) pvals[phase2Ind] <- stats::pexp(q = q - delay2 + rate1/rate2 * (delay2 - delay1), rate = rate2, ...)[phase2Ind]

  pvals
}

#' @rdname DelayedExponential
#' @export
qexp_delayed <- function(p, delay1 = 0, rate1 = 1, delay2 = NULL, rate2 = NULL, ...) {
  stopifnot( is.numeric(delay1), is.numeric(rate1), all(is.finite(delay1)), all(is.finite(rate1)) )

  qvals <- delay1 + stats::qexp(p = p, rate = rate1, ...)
  # check for easy case: only a single delay
  if ( is.null(delay2) ) return(qvals)

  # we need both delay2 AND rate2
  stopifnot( is.numeric(delay2), is.numeric(rate2) )

  # manually recycle delay and rate parameters to ensure they have common length, respectively
  ndelay <- max(length(delay1), length(delay2), na.rm = TRUE)
  nrate <- max(length(rate1), length(rate2), na.rm = TRUE)

  delay1 <- rep_len(delay1, length.out = ndelay); delay2 <- rep_len(delay2, length.out = ndelay)
  rate1 <- rep_len(rate1, length.out = nrate); rate2 <- rep_len(rate2, length.out = nrate)

  # check delay constraint
  if (any(delay1 >= delay2)) {
    stop("First delay phase must always antedate the second delay phase!")
  }

  # check if we have observations in 2nd phase
  phase2Ind <- which(p >= pexp_delayed(q = delay2, delay1 = delay1, rate1 = rate1))
  if (length(phase2Ind)) qvals[phase2Ind] <- delay2 - rate1/rate2 * (delay2 - delay1) + stats::qexp(p = p, rate = rate2, ...)[phase2Ind]

  qvals

}

#' @rdname DelayedExponential
#' @export
rexp_delayed <- function(n, delay1 = 0, rate1 = 1, delay2 = NULL, rate2 = NULL) {
  stopifnot( is.numeric(delay1), is.numeric(rate1), all(is.finite(delay1)), all(is.finite(rate1)) )

  # check for easy case: only a single delay
  if ( is.null(delay2) ) return(delay1 + stats::rexp(n = n, rate = rate1))

  # we need both delay2 AND rate2
  stopifnot( is.numeric(delay2), is.numeric(rate2) )

  # manually recycle delay and rate parameters to ensure they have common length, respectively
  ndelay <- max(length(delay1), length(delay2), na.rm = TRUE)
  nrate <- max(length(rate1), length(rate2), na.rm = TRUE)

  delay1 <- rep_len(delay1, length.out = ndelay); delay2 <- rep_len(delay2, length.out = ndelay)
  rate1 <- rep_len(rate1, length.out = nrate); rate2 <- rep_len(rate2, length.out = nrate)

  # check delay constraint
  if (any(delay1 >= delay2)) {
    stop("First delay phase must always antedate the second delay phase!")
  }

  # use inverse CDF-method
  qexp_delayed(p = stats::runif(n = n), delay1 = delay1, rate1 = rate1, delay2 = delay2, rate2 = rate2)

}


#' @rdname DelayedExponential
#' @export
mexp_delayed <- function(t=+Inf, delay1 = 0, rate1 = 1, delay2 = NULL, rate2 = NULL) {
  stopifnot( is.numeric(delay1), is.numeric(rate1), all(is.finite(delay1)), all(is.finite(rate1)) )

  # calculate for single phase delayed exponential
  res_mvals <- pmin.int(t, delay1) + pexp_delayed(q = if (is.null(delay2)) t else pmin.int(t, delay2),
                                                 delay1 = delay1, rate1 = rate1) / rate1

  # directly return when only single phase
  if ( is.null(delay2) ) return(res_mvals)

  # we need both delay2 AND rate2
  stopifnot( is.numeric(delay2), is.numeric(rate2) )

  # manually recycle delay and rate parameters to ensure they have common length, respectively
  ndelay <- max(length(delay1), length(delay2), na.rm = TRUE)
  nrate <- max(length(rate1), length(rate2), na.rm = TRUE)

  delay1 <- rep_len(delay1, length.out = ndelay); delay2 <- rep_len(delay2, length.out = ndelay)
  rate1 <- rep_len(rate1, length.out = nrate); rate2 <- rep_len(rate2, length.out = nrate)

  # check delay constraint
  if (any(delay1 >= delay2)) {
    stop("First delay phase must always antedate the second delay phase!")
  }

  res_mvals + (1L-pexp_delayed(q=delay2, delay1 = delay1, rate1 = rate1)) * pexp_delayed(q=t, delay1 = delay2, rate1 = rate2) / rate2
}


#' Delayed Weibull Distribution
#'
#' @description
#' Density, distribution function, quantile function and random generation for the delayed Weibull distribution with parameters
#' as in the Weibull distribution functions in R's stats-package, namely:
#' * `delay`
#' * `shape`
#' * `scale` (inverse of rate)
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
#' @param delay numeric. The delay, must be non-negative.
#' @param shape numeric. Shape parameter, must be positive.
#' @param scale numeric. Scale parameter (inverse of rate), must be positive.
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
dweib_delayed <- function(x, delay1, shape1, scale1 = 1, ...) stats::dweibull(x = x - delay1, shape = shape1, scale = scale1, ...)
#' @rdname DelayedWeibull
#' @export
pweib_delayed <- function(q, delay1, shape1, scale1 = 1, ...) stats::pweibull(q = q - delay1, shape = shape1, scale = scale1, ...)
#' @rdname DelayedWeibull
#' @export
qweib_delayed <- function(p, delay1, shape1, scale1 = 1, ...) delay1 + stats::qweibull(p = p, shape = shape1, scale = scale1, ...)
#' @rdname DelayedWeibull
#' @export
rweib_delayed <- function(n, delay1, shape1, scale1 = 1) delay1 + stats::rweibull(n = n, shape = shape1, scale = scale1)
#' @rdname DelayedWeibull
#' @export
mweib_delayed <- function(t=+Inf, delay1, shape1, scale1 = 1, ...) ifelse(test = t <= delay1, yes = t,
                                                                       # make use of lower incomplete gamma function which is calculated as gamma * pgamma
                                                                       no = delay1 + scale1 / shape1 * gamma(1/shape1) * stats::pgamma(q = ((t-delay1)/scale1)^shape1, shape = 1/shape1))


#' Get delay distribution function
#' @param distribution character(1). delay distribution.
#' @param type character(1). type of function, cdf: cumulative distribution function, density or random function
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
           # # XXXX quick fix for 2-phase param
           # if (distribution == 'exponential' && ! twoGroup) {
           #   stopifnot( ! twoGroup ) # XXX implication of twoGroup=TRUE and bind=.. are not implemented yet!
           #
           #   #return(c("delay", "rate", "delay2", "rate2"))
           #   return(rownames(if (isTRUE(transformed)) PARAM_TRANSF_EXP_MAT else PARAM_TRANSF_EXP_INV)[seq_len(2L*(1L + twoPhase))])
           # }

           par_list <- list(exponential = c("delay1", "rate1", "delay2", "rate2")[seq_len(2L*(1L + twoPhase))],
                            weibull = c("delay1", "shape1", "scale1", "delay2", "shape2", "scale2")[seq_len(3L*(1L + twoPhase))])

           if (transformed){
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

