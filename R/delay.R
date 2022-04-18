# mkuhn, 2021-04-06
# delay distribution functions


#' Density of delayed exponential distribution.
#' @param x numeric. Values for which to get the density.
#' @param delay numeric. The delay, must be non-negative.
#' @param rate numeric. The event rate, must be non-negative.
#' @param ... further arguments to `stats::dexp`
#' @export
dexp_delayed <- function(x, delay, rate = 1, ...) dexp(x = x - delay, rate = rate, ...)
#' @export
pexp_delayed <- function(q, delay, rate = 1, ...) pexp(q = q - delay, rate = rate, ...)
#' @export
rexp_delayed <- function(n, delay, rate = 1) rexp(n = n, rate = rate) + delay

#' Density of delayed Weibull distribution
#' @param x numeric. Values for which to get the density.
#' @param delay numeric. The delay, must be non-negative.
#' @param shape numeric. Shape parameter, must be positive.
#' @param scale numeric. Scale parameter (inverse of rate), must be positive.
#' @param ... further arguments to `stats::deweibull`
#' @export
dweib_delayed <- function(x, delay, shape, scale = 1, ...) dweibull(x = x - delay, shape = shape, scale = scale, ...)
#' @export
pweib_delayed <- function(q, delay, shape, scale = 1, ...) pweibull(q = q - delay, shape = shape, scale = scale, ...)
#{cat(sprintf('min: %.3f, delay: %.3f — shape: %.2f — scale: %.2f\n', min(q), delay, shape, scale)); }
#' @export
rweib_delayed <- function(n, delay, shape, scale = 1) rweibull(n = n, shape = shape, scale = scale) + delay


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
         stop(glue::glue("Unknown attribute of distribution {distribution}."))
  )[[1L + (distribution == 'weibull')]]
}

