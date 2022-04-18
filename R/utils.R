# mkuhn, 2021-10-11
# Utility functions used within this package

#' Compare two numeric vectors
#' @param x numeric vector
#' @param y numeric vector
#' @param tol numeric. Numerical tolerance level for comparison.
#' @seealso `dplyr::near`
near <- function (x, y, tol = .Machine$double.eps^0.4){
  abs(x - y) < tol
}

#' Estimate rounding error based on given sample of metric values
#' The idea is to check at which level of rounding the sample values do not change.
#' @param obs numeric. Metric values from a sample to estimate the corresponding rounding error
#' @param roundDigits integer. Which level of rounding to test? Negative numbers round to corresponding powers of 10
#' @param maxObs integer. How many observations to consider at most? If the provided sample has more observations a sub-sample is used.
#' @return estimated rounding error
estimRoundingError <- function(obs, roundDigits = seq.int(-4L, 6L), maxObs = 100L) {
  stopifnot( is.numeric(maxObs), length(maxObs) == 1L )
  maxObs <- trunc(maxObs)
  # drop NA and Inf
  obs <- obs[is.finite(obs)]

  if (maxObs > 1L && length(obs) > maxObs){
    obs <- obs[round(seq.int(from = 1, to = length(obs), length.out = maxObs))]
  }


  # digits to round to
  roundDigits <- unique(trunc(roundDigits))

  #XXX think here: fails for small obs close to 0!
  if (all(abs(obs) < .1)) obs <- 1L + obs

  rDigInd <- purrr::map_lgl(.x = roundDigits,
                            .f = ~all(near(x = obs, y = round(obs, digits = .x),
                                           tol = 2L * 10L**min(-2L, -.x-1L))) )
  10L**-if (!any(rDigInd)) max(roundDigits)+1L else if (all(rDigInd)) min(roundDigits)-1L else roundDigits[which.max(rDigInd)]
}
