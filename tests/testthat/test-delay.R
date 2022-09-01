# mkuhn, 2022-03-16
# tests for the delayed distribution functions


test_that('density of delayed distributions', {
  tPoints <- seq.int(from = -1, to = 11, length.out = 27)

  rateVals <- c(.06, .32, .821, 1.14, 1.78, 2.19, 5.116, 11.2, rlnorm(n=3, meanlog = 1, sdlog = 3))
  shapeVals <- c(.0001, .0052, .14, .87, 1.26, 3.9, 14.8, 21, rpois(n=3, lambda = 8))

  # delayed exponential with delay=0 coincides with stats::dexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(dexp_delayed(x = tPoints, delay = 0L, rate = .x),
                                                    stats::dexp(x = tPoints, rate = .x)))

  # delayed Weibull with delay=0 coincides with stats::dweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
              .f = ~expect_identical(dweib_delayed(x = tPoints, delay = 0L, scale = .x, shape = .y),
                                     stats::dweibull(x = tPoints, scale = .x, shape = .y)))

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  delayTimes <- c(0, 1, 4, 7, 11, 13) + .11
  delayTimes <- unique(c(delayTimes, rpois(n=length(rateVals) - length(delayTimes), lambda = 11)))
  delayTimes <- c(delayTimes, rnorm(n = length(rateVals) - length(delayTimes), mean = 7, sd = .71))
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(dexp_delayed(x = tPoints, delay = .x, rate = .y),
                                  dweib_delayed(x = tPoints, delay = .x, shape = 1L, scale = .y**-1)))

  # delayed exponential density differs by the factor exp(lambda * alpha) from the standard delay density
  #+ on log-scale, which is numerically more robust, this becomes a difference by lambda * alpha
  purrr::walk2(.x  = delayTimes, .y = rateVals,
               .f = ~ expect_equal(dexp_delayed(x = tPoints, delay = .x, rate = .y, log=TRUE),
                                   dplyr::if_else(tPoints < .x, true = -Inf, false = dexp(x = tPoints, rate = .y, log = TRUE) + .x * .y),
                                   tolerance = .01))

})


test_that('distribution function of delayed distributions', {
  qPoints <- seq.int(from = 0, to = 1, length.out = 27)

  rateVals <- c(.1, .52, 1.4, 1.58, 3.9, 11.2)
  shapeVals <- c(.001, .14, .84, 1.6, 5.9, 13.8)
  # CDF of delayed exponential with delay=0 coincides with stats::pexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(pexp_delayed(q = qPoints, delay = 0L, rate = .x),
                                                           stats::pexp(q = qPoints, rate = .x)))

  # CDF of delayed Weibull with delay=0 coincides with stats::pweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
               .f = ~expect_identical(pweib_delayed(q = qPoints, delay = 0L, scale = .x, shape = .y),
                                      stats::pweibull(q = qPoints, scale = .x, shape = .y)))

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  delayTimes <- c(0, 1, 4, 7, 11, 13)
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints, delay = .x, rate = .y),
                                  pweib_delayed(q = qPoints, delay = .x, shape = 1L, scale = .y**-1)))
})


test_that('(restricted) mean survival time of delayed distributions', {
  # restriction within delay period:
  earlyTP <- seq.int(from = 0, to = 5, length.out = 11)
  expect_identical(mexp_delayed(t = earlyTP, delay = max(earlyTP), rate = 1/pi), expected = earlyTP)
  expect_identical(mweib_delayed(t = earlyTP, delay = max(earlyTP), shape = 1/pi, scale = 1/pi), expected = earlyTP)

  # restriction longer than delay period (for exponential)
  expect_identical(mexp_delayed(t = 5 + earlyTP, delay = 5, rate = 1/pi), expected = 5 + pi * pexp_delayed(q = 5 + earlyTP, delay = 5, rate = 1/pi))

  tPoints <- c(seq.int(from = 6, to = 11, length.out = 13L), +Inf)

  # some random delay value
  delayTimes <- runif(n = 7L, min = 0, max = 7)

  settingDF <- tidyr::expand_grid(tPoints = c(seq.int(from = 6, to = 11, length.out = 6), +Inf),
                     delayTimes = runif(n = 7, min = 0, max = 7),
                     rates = c(.1, .5, 1, 2, 5))

  # restricted mean survival of exponential coincides with that from weibull with shape =1
  purrr::pwalk(.l = settingDF,
               .f = ~expect_equal(mexp_delayed(t = ..1, delay = ..2, rate = ..3),
                                  mweib_delayed(t = ..1, delay = ..2, shape = 1L, scale = ..3**-1)))

  # simulated data
  set.seed(1234)
  simObs1 <- rweib_delayed(n = 123456, delay = 5, shape = 0.5, scale = 2)
  simObs2 <- rweib_delayed(n = 123456, delay = 5, shape = 2.1, scale = 2)

  expect_equal(mean(pmin.int(simObs1, 6.5)), expected = mweib_delayed(t = 6.5, delay = 5, shape = 0.5, scale = 2), tolerance = .003)
  expect_equal(mean(pmin.int(simObs2, 6.5)), expected = mweib_delayed(t = 6.5, delay = 5, shape = 2.1, scale = 2), tolerance = .003)
})

