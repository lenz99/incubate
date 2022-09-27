# mkuhn, 2022-03-16
# tests for the delayed distribution functions


test_that('density of delayed distributions', {
  tPoints <- seq.int(from = -1, to = 11, length.out = 27)

  rateVals <- c(.06, .32, .821, 1.14, 1.78, 2.19, 5.116, 11.2, rlnorm(n=3, meanlog = 1, sdlog = 3))
  shapeVals <- c(.0001, .0052, .14, .87, 1.26, 3.9, 14.8, 21, rpois(n=3, lambda = 8))


  # delayed exponential with delay=0 coincides with stats::dexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(dexp_delayed(x = tPoints, delay1 = 0L, rate1 = .x),
                                                    stats::dexp(x = tPoints, rate = .x)))

  # delayed Weibull with delay=0 coincides with stats::dweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
              .f = ~expect_identical(dweib_delayed(x = tPoints, delay1 = 0L, scale1 = .x, shape1 = .y),
                                     stats::dweibull(x = tPoints, scale = .x, shape = .y)))

  delayTimes <- c(0, 1, 4, 7, 11, 13) + .11
  delayTimes <- unique(c(delayTimes, rpois(n=length(rateVals) - length(delayTimes), lambda = 11)))
  delayTimes <- c(delayTimes, rnorm(n = length(rateVals) - length(delayTimes), mean = 7, sd = .71))

  # either argument delay=/rate= or delay1=/rate1= can be used
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_identical(dexp_delayed(x = tPoints, delay = .x, rate = .y),
                                      dexp_delayed(x = tPoints, delay1 = .x, rate1 = .y)))
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_identical(dexp_delayed(x = tPoints, delay = .x, rate = .y),
                                      dexp_delayed(x = tPoints, delay1 = .x, rate = .y)))
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_identical(dexp_delayed(x = tPoints, delay = .x, rate1 = .y),
                                      dexp_delayed(x = tPoints, delay1 = .x, rate = .y)))
  # using both delay= and delay1= results in a warning
  expect_warning(dexp_delayed(x = tPoints, delay = 5, delay1 = 8), regexp = "is ignored")
  # using both rate= and rate1= results in a warning
  expect_warning(dexp_delayed(x = tPoints, delay = 3, rate = 1.1, rate1 = .8), regexp = "is ignored")

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(dexp_delayed(x = tPoints,  delay1 = .x, rate1 = .y),
                                  dweib_delayed(x = tPoints, delay1 = .x, shape1 = 1L, scale1 = .y**-1)))

  # using both delay= and delay1= results in a warning
  expect_warning(dweib_delayed(x = tPoints, delay = 5, delay1 = 8, shape = .2), regexp = "is ignored")
  # using both scale= and scale1= results in a warning
  expect_warning(dweib_delayed(x = tPoints, delay = 3, shape = 1.8, scale = 1.1, scale1 = .8), regexp = "is ignored")
  # using both shape= and shape1= results in a warning
  expect_warning(dweib_delayed(x = tPoints, delay = 3, shape = 1.8, shape1 = 1.3, scale = .8), regexp = "is ignored")

  # delayed exponential density differs by the factor exp(lambda * alpha) from the standard delay density
  #+ on log-scale, which is numerically more robust, this becomes a difference by lambda * alpha
  purrr::walk2(.x  = delayTimes, .y = rateVals,
               .f = ~ expect_equal(dexp_delayed(x = tPoints, delay1 = .x, rate1 = .y, log = TRUE),
                                   dplyr::if_else(tPoints < .x, true = -Inf, false = dexp(x = tPoints, rate = .y, log = TRUE) + .x * .y),
                                   tolerance = .01))


  # examples for 1-phase and 2-phase Weibull density with scale-parameters
  expect_equal(dweib_delayed(x = 2.5, delay1 = 1, shape1 = .7, scale1 = 2),
               expected = .7 / 2 * ((2.5-1)/2)^(.7-1) * exp(-((2.5-1)/2)**.7))
  expect_equal(dweib_delayed(x = 5.5, delay1 = 1, shape1 = .7, scale1 = 2, delay2 = 4, shape2 = 1.4, scale2 = 1.7),
               expected = exp(- ((4 - 1)/2)^.7) * 1.4/1.7 * ((5.5 - 4)/1.7)^(1.4-1) * exp(-((5.5 - 4)/1.7)^1.4))
})


test_that('distribution function of delayed distributions', {
  qPoints <- seq.int(from = 0, to = 1, length.out = 27)

  rateVals <- c(.1, .52, 1.4, 1.58, 3.9, 11.2)
  shapeVals <- c(.001, .14, .84, 1.6, 5.9, 13.8)
  # CDF of delayed exponential with delay=0 coincides with stats::pexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(pexp_delayed(q = qPoints, delay1 = 0L, rate1 = .x),
                                                           stats::pexp(q = qPoints, rate = .x)))

  # CDF of delayed Weibull with delay=0 coincides with stats::pweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
               .f = ~expect_identical(pweib_delayed(q = qPoints, delay1 = 0L, scale1 = .x, shape1 = .y),
                                      stats::pweibull(q = qPoints, scale = .x, shape = .y)))

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  delayTimes <- c(0, 1, 4, 7, 11, 13)
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints,  delay1 = .x, rate1 = .y),
                                  pweib_delayed(q = qPoints, delay1 = .x, shape1 = 1L, scale1 = .y**-1)))
})


test_that('(restricted) mean survival time of delayed distributions', {
  # restriction within delay period:
  earlyTP <- seq.int(from = 0, to = 5, length.out = 11)
  expect_identical(mexp_delayed(t = earlyTP, delay1 = max(earlyTP), rate1 = 1/pi),
                   expected = earlyTP)
  # two-phase exponential
  expect_identical(mexp_delayed(t = earlyTP, delay1 = max(earlyTP), rate1 = 1/pi, delay2 = max(earlyTP)+1, rate2 = 2/pi),
                   expected = earlyTP)
  expect_identical(mweib_delayed(t = earlyTP, delay1 = max(earlyTP), shape1 = 1/pi, scale1 = 1/pi),
                   expected = earlyTP)

  # restriction longer than delay period (for exponential)
  expect_equal(mexp_delayed(t = 5 + earlyTP, delay1 = 5, rate1 = 1/pi),
               expected = 5 + pi * pexp_delayed(q = 5 + earlyTP, delay1 = 5, rate1 = 1/pi))

  # 2-phase
  expect_equal(mexp_delayed(t = 5 + earlyTP, delay1 = 5, rate1 = 1/pi, delay2 = 7.5, rate2 = 2/pi),
               expected = 5 + pi * pexp_delayed(q = pmin.int(5 + earlyTP, 7.5), delay1 = 5, rate1 = 1/pi) +
                 pi/2 * pexp_delayed(q = 5 + earlyTP, delay1 = 7.5, rate1 = 2/pi) * exp(-1/pi * (7.5 - 5)))


  settingDF <- tidyr::expand_grid(tPoints = c(seq.int(from = 6, to = 11, length.out = 7), +Inf),
                                  # some random delay value
                                  delayTimes = runif(n = 7, min = 0, max = 7),
                                  rates = c(.1, .5, 1, 2, 5))

  # restricted mean survival of exponential coincides with that from Weibull with shape =1
  purrr::pwalk(.l = settingDF,
               .f = ~expect_equal(mexp_delayed(t = ..1,  delay1 = ..2, rate1 = ..3),
                                  mweib_delayed(t = ..1, delay1 = ..2, shape = 1L, scale1 = ..3**-1)))

  # simulated Weibull data
  expect_equal(mweib_delayed(t = 6.5, delay1 = 5, shape1 = 0.5, scale1 = 2),
               expected = mean(pmin.int(rweib_delayed(n = 123456, delay1 = 5, shape1 = 0.5, scale1 = 2), 6.5)),
               tolerance = .003)
  expect_equal(mweib_delayed(t = 6.5, delay1 = 5, shape1 = 2.1, scale1 = 2),
               expected = mean(pmin.int(rweib_delayed(n = 123456, delay = 5, shape1 = 2.1, scale1 = 2), 6.5)),
               tolerance = .003)
})

