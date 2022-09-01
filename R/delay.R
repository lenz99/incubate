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
#' @param delay numeric. The delay, must be non-negative.
#' @param rate numeric. The event rate, must be non-negative.
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
dexp_delayed <- function(x, delay = 0, rate = 1, delay2 = NULL, rate2 = NULL, log = FALSE) {
  stopifnot( is.numeric(delay), is.numeric(rate), all(is.finite(delay)), all(is.finite(rate)) )
  log <- isTRUE(log)

  # check for easy case: only a single delay
  if ( is.null(delay2) ) return(stats::dexp(x = x - delay, rate = rate, log = log))

  # we need both delay2 AND rate2
  stopifnot( is.numeric(delay2), is.numeric(rate2) )

  # recycle delay and rate parameters to have common length, respectively
  ndelayMax <- max(length(delay), length(delay2), na.rm = TRUE)
  nrateMax <- max(length(rate), length(rate2), na.rm = TRUE)

  delay <- rep_len(delay, length.out = ndelayMax)
  delay2 <- rep_len(delay2, length.out = ndelayMax)

  rate <- rep_len(rate, length.out = nrateMax)
  rate2 <- rep_len(rate2, length.out = nrateMax)

  # check delay constraint
  if (any(delay >= delay2)) {
    stop("First delay must always antedate the second delay!")
  }

  # normalizing constant
  K <- 1L - exp(-rate * (delay2 - delay)) + exp(-rate2 * (delay2 - delay))
  stopifnot( all(K > 0L) )


  dvals <- stats::dexp(x = x - delay, rate = rate, log = log)
  # check if we have observations in 2nd phase
  phase2Ind <- which(x >= delay2)
  if (length(phase2Ind)) dvals[phase2Ind] <- stats::dexp(x = x - delay, rate = rate2, log = log)[phase2Ind]

  if (log) dvals - log(K) else dvals/K

  ## XXX think & continue here
  #+ density with lower rate starts lower but has smaller decline and will eventually lie above the density with higher rate
  #+ this means, our density_delayed can jump upwards at delay2 even when rate goes down (rate2 < rate)!!
  #+ Is this what we want?
  # ggplot() + xlim(0, 8) +
  ## geom_function(fun = dexp, args = list(rate = .5)) +
  ##geom_function(fun = dexp, args = list(rate = .2), col = "red")
}

#' @rdname DelayedExponential
#' @export
pexp_delayed <- function(q, delay = 0, rate = 1, ...) stats::pexp(q = q - delay, rate = rate, ...)
#' @rdname DelayedExponential
#' @export
qexp_delayed <- function(p, delay = 0, rate = 1, ...) delay + stats::qexp(p = p, rate = rate, ...)
#' @rdname DelayedExponential
#' @export
rexp_delayed <- function(n, delay = 0, rate = 1) delay + stats::rexp(n = n, rate = rate)
#' @rdname DelayedExponential
#' @export
mexp_delayed <- function(t=+Inf, delay = 0, rate = 1, ...) ifelse(test = t <= delay, yes = t,
                                                              no = delay + 1/rate * pexp_delayed(q=t, delay=delay, rate = rate))


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
dweib_delayed <- function(x, delay, shape, scale = 1, ...) stats::dweibull(x = x - delay, shape = shape, scale = scale, ...)
#' @rdname DelayedWeibull
#' @export
pweib_delayed <- function(q, delay, shape, scale = 1, ...) stats::pweibull(q = q - delay, shape = shape, scale = scale, ...)
#' @rdname DelayedWeibull
#' @export
qweib_delayed <- function(p, delay, shape, scale = 1, ...) delay + stats::qweibull(p = p, shape = shape, scale = scale, ...)
#' @rdname DelayedWeibull
#' @export
rweib_delayed <- function(n, delay, shape, scale = 1) delay + stats::rweibull(n = n, shape = shape, scale = scale)
#' @rdname DelayedWeibull
#' @export
mweib_delayed <- function(t=+Inf, delay, shape, scale = 1, ...) ifelse(test = t <= delay, yes = t,
                                                                       # make use of lower incomplete gamma function which is calculated as gamma * pgamma
                                                                       no = delay + scale / shape * gamma(1/shape) * stats::pgamma(q = ((t-delay)/scale)^shape, shape = 1/shape))


#' Get delay distribution function
#' @param distribution character(1). delay distribution.
#' @param type character(1). type of function, cdf: cumulative distribution function, density or random function
#' @param twoGroup logical(1). Do we have two groups?
#' @param bind character. Names of parameters that are bind between the two groups.
#' @return selected distribution function or parameter names
getDist <- function(distribution = c("exponential", "weibull"), type = c("cdf", "prob", "density", "random", "param"),
                    twoGroup = FALSE, bind = NULL) {
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
           par_list <- list(exponential = c("delay", "rate"),
                            weibull = c("delay", "shape", "scale"))
           if (twoGroup) {
             # bind only parameters from chosen distribution
             myPars <- par_list[[1L + (distribution == 'weibull')]]
             bind <- intersect(myPars, bind) #intersect: enforces original order
             par_gr <- purrr::map(par_list, ~ setdiff(., bind))
             # bind parameters first
             purrr::map(par_gr, ~ c(bind,
                                    paste(rep(., 2L), rep(c("x", "y"), each = length(.)),
                                          sep = ".")))
           } else par_list

         },
         stop(glue("Unknown attribute of distribution {distribution}."), call. = FALSE)
  )[[1L + (distribution == 'weibull')]]
}

