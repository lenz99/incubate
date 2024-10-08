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
  qPoints <- seq.int(from = 0, to = 8, length.out = 27)

  rateVals <- c(.1, .52, 1.4, 1.58, 3.9, 11.2)
  shapeVals <- c(.001, .14, .84, 1.6, 5.9, 13.8)
  # CDF of delayed exponential with delay=0 coincides with stats::pexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(pexp_delayed(q = qPoints, delay1 = 0L, rate1 = .x),
                                                    expected = stats::pexp(q = qPoints, rate = .x)))

  # CDF of delayed Weibull with delay=0 coincides with stats::pweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
               .f = ~expect_identical(pweib_delayed(q = qPoints, delay1 = 0L, scale1 = .x, shape1 = .y),
                                      expected = stats::pweibull(q = qPoints, scale = .x, shape = .y)))

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  delayTimes <- c(0, 1, 4, 7, 11, 13)
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints,  delay1 = .x, rate1 = .y),
                                  pweib_delayed(q = qPoints, delay1 = .x, shape1 = 1L, scale1 = .y**-1)))

  # lower.tail=T/F adds up to 1
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, lower.tail = TRUE),
                                      expected = 1L - pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, lower.tail = FALSE)))

  # 2-phase
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, delay2 = .x + 1, rate2 = sqrt(.y), lower.tail = TRUE),
                                      expected = 1L - pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, delay2 = .x + 1, rate2 = sqrt(.y), lower.tail = FALSE)))

  # log.p=
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, log.p = TRUE),
                                  expected = log(pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, log.p = FALSE))))

  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, delay2 = .x + 1, rate2 = sqrt(.y), log.p = TRUE),
                                  expected = log(pexp_delayed(q = qPoints, delay1 = .x, rate1 = .y, delay2 = .x + 1, rate2 = sqrt(.y), log.p = FALSE))))
})


test_that('quantile function (percentage point function) of delayed distributions', {
  del1 <- 1
  pVals <- seq.int(from = 0, to = 1, length.out = 19)
  expect_equal(qexp_delayed(p = pVals, delay1 = del1, rate1 = .67),
               expected = del1 - log(1-pVals)/.67)
  expect_equal(qweib_delayed(p = pVals, delay1 = 1, shape1 = 2, scale1 = 1.5),
               expected = del1 + sqrt(-log(1-pVals)) * 1.5)

  # lower.tail=
  purrr::walk(.x = c(.003, .1, .35, .51, 1, 4, 11),
              .f = ~{
                expect_equal(qexp_delayed(p = pVals, delay1 = del1, rate1 = .x, lower.tail = TRUE),
                             expected = qexp_delayed(p = 1-pVals, delay1 = del1, rate1 = .x, lower.tail = FALSE))
                expect_equal(qexp_delayed(p = pVals, delay1 = del1, rate1 = .x, delay2 = del1+.89, rate2 = .34, lower.tail = TRUE),
                             expected = qexp_delayed(p = 1-pVals, delay1 = del1, rate1 = .x, delay2 = del1+.89, rate2 = .34, lower.tail = FALSE))
              })

  purrr::walk2(.x = c(.003, .1, .35, .51, 1, 4, 11), .y = sample(c(.0001, .11, .5, 1, 2.1, pi, exp(2))),
              .f = ~ {
                expect_equal(qweib_delayed(p = pVals, delay1 = del1, shape1 = .x, scale1 = .y, lower.tail = TRUE),
                             expected = qweib_delayed(p = 1-pVals, delay1 = del1, shape1 = .x, scale1 = .y, lower.tail = FALSE))
                expect_equal(qweib_delayed(p = pVals, delay1 = del1, shape1 = .x, scale1 = .y, delay2 = del1+.93, shape2 = .9, scale2 = 1.1, lower.tail = TRUE),
                             expected = qweib_delayed(p = 1-pVals, delay1 = del1, shape1 = .x, scale1 = .y, delay2 = del1+.93, shape2 = .9, scale2 = 1.1, lower.tail = FALSE))
                })

  # log.p= works
  purrr::walk(.x = seq.int(from = 0, to = 1, length.out = 19),
              .f = ~expect_equal(qexp_delayed(p = log(.x), delay1 = del1, rate1 = .67, log.p = TRUE),
                                     expected = qexp_delayed(p = .x, delay1 = del1, rate1 = .67, log.p = FALSE)))
  purrr::walk(.x = seq.int(from = 0, to = 1, length.out = 19),
              .f = ~expect_equal(qweib_delayed(p = log(.x), delay1 = del1, shape1 = .67, log.p = TRUE),
                                     expected = qweib_delayed(p = .x, delay1 = del1, shape1 = .67, log.p = FALSE)))
  # 2-phase
  purrr::walk(.x = seq.int(from = 0, to = 1, length.out = 19),
              .f = ~expect_equal(qexp_delayed(p = log(.x), delay1 = del1, rate1 = .67, delay2 = del1 + 1, rate2 = .77, log.p = TRUE),
                                     expected = qexp_delayed(p = .x, delay1 = del1, rate1 = .67, delay2 = del1 + 1, rate2 = .77, log.p = FALSE)))
  purrr::walk(.x = seq.int(from = 0, to = 1, length.out = 19),
              .f = ~expect_equal(qweib_delayed(p = log(.x), delay1 = del1, shape1 = .67, delay2 = del1 + 1, shape2 = 1.1, scale2 = .8, log.p = TRUE),
                                     expected = qweib_delayed(p = .x, delay1 = del1, shape1 = .67, delay2 = del1 + 1, shape2 = 1.1, scale2 = .8, log.p = FALSE)))


  # qweib inverts pweib
  qVals <- seq.int(0, 5)
  expect_equal(qweib_delayed(p = pweib_delayed(q = qVals, delay1 = del1, shape1 = 1.5, scale1 = .8),
                             delay1 = del1, shape1 = 1.5, scale1 = .8),
               expected = pmax.int(qVals, del1)) # qweib_delayed can't go below delay1
  # two-phase
  expect_equal(qweib_delayed(p = pweib_delayed(q = qVals, delay1 = del1, shape1 = 1.5, scale1 = .8,
                                               delay2 = del1+.5, shape2 = .89, scale2 = .61),
                             delay1 = del1, shape1 = 1.5, scale1 = .8,
                             delay2 = del1+.5, shape2 = .89, scale2 = .61),
               expected = pmax.int(qVals, del1)) # can't go below delay1
  # pweib inverts qweib
  pVals <- seq.int(0, 1, length.out = 13L)
  expect_equal(pweib_delayed(q = qweib_delayed(p = pVals, delay1 = del1, shape1 = 1.5, scale1 = .8),
                             delay1 = del1, shape1 = 1.5, scale1 = .8),
               expected = pVals)
  # two-phase
  expect_equal(pweib_delayed(q = qweib_delayed(p = pVals, delay1 = del1, shape1 = 1.5, scale1 = .8,
                                               delay2 = del1+.5, shape2 = .89, scale2 = .61),
                             delay1 = del1, shape1 = 1.5, scale1 = .8,
                             delay2 = del1+.5, shape2 = .89, scale2 = .61),
               expected = pVals)
})


test_that('(restricted) mean survival time of delayed distributions', {
  withr::local_package("tidyr")
  withr::local_package("dplyr")

  # restriction within delay period:
  earlyTP <- seq.int(from = 0, to = 5, length.out = 11L)

  expect_identical(mexp_delayed(t = earlyTP, delay1 = max(earlyTP), rate1 = 1/pi),
                   expected = earlyTP)
  expect_identical(mexp_delayed(t = earlyTP, delay1 = max(earlyTP), rate1 = 1/pi, delay2 = max(earlyTP)+1, rate2 = 2/pi),
                   expected = earlyTP) # 2-phase exponential
  expect_identical(mweib_delayed(t = earlyTP, delay1 = max(earlyTP), shape1 = 1/pi, scale1 = sqrt(pi)),
                   expected = earlyTP)
  expect_identical(mweib_delayed(t = earlyTP, delay1 = max(earlyTP), shape1 = 1/pi, scale1 = sqrt(pi), delay2 = max(earlyTP) + 1, shape2 = exp(-1), scale2 = .37),
                   expected = earlyTP) # 2-phase Weibull

  # restriction beyond the delay period (for exponential)
  expect_equal(mexp_delayed(t = 5 + earlyTP, delay1 = 5, rate1 = 1/pi),
               expected = 5 + pi * pexp_delayed(q = 5 + earlyTP, delay1 = 5, rate1 = 1/pi))

  # 2-phase
  expect_equal(mexp_delayed(t = 5 + earlyTP, delay1 = 5, rate1 = 1/pi, delay2 = 7.5, rate2 = 2/pi),
               expected = 5 + pi * pexp_delayed(q = pmin.int(5 + earlyTP, 7.5), delay1 = 5, rate1 = 1/pi) +
                 pi/2 * pexp_delayed(q = 5 + earlyTP, delay1 = 7.5, rate1 = 2/pi) * exp(-1/pi * (7.5 - 5)))


  settingDF <- tidyr::expand_grid(# some random delay value
                                  delay1 = runif(n = 3, min = 0, max = 7),
                                  rates1 = c(.1, .6, 1, 2, 5),
                                  delay2 = delay1 + runif(n = length(delay1), min = 0.5, max = 5),
                                  rates2 = c(.11, .31, .91, 1.51, 3.21))
  settingDF <- settingDF[settingDF$delay1 < settingDF$delay2,]

  tPoints <-  c(seq.int(from = 6, to = 11, length.out = 7L), +Inf)

  # restricted mean survival can never exceed the restricted time
  #+as survival curve never goes beyond 1, so integral is less than length of time-axis
  purrr::pwalk(.l = tidyr::expand_grid(tPoints, unique(settingDF[,c("delay1", "rates1")])),
               .f = ~expect_lte(mexp_delayed(t = ..1,  delay1 = ..2, rate1 = ..3),
                                expected = ..1))
  purrr::pwalk(.l = tidyr::expand_grid(tPoints, settingDF %>% dplyr::slice_sample(n = 11)),
               .f = ~expect_lte(mexp_delayed(t = ..1,  delay1 = ..2, rate1 = ..3, delay2 = ..4, rate2 = ..5),
                                expected = ..1))

  # restricted mean survival of exponential coincides with that from Weibull with shape =1
  purrr::pwalk(.l = unique(settingDF[,c("delay1", "rates1")]),
               .f = ~expect_equal(mexp_delayed(t = tPoints,  delay1 = ..1, rate1 = ..2),
                                  mweib_delayed(t = tPoints, delay1 = ..1, shape1 = 1L, scale1 = ..2**-1)))
  # 2-phase
  purrr::pwalk(.l = settingDF,
               .f = ~expect_equal(mexp_delayed(t = tPoints,  delay1 = ..1, rate1 = ..2, delay2 = ..3, rate2 = ..4),
                                  mweib_delayed(t = tPoints, delay1 = ..1, shape1 = 1, scale1 = ..2**-1, delay2 = ..3, shape2 = 1, scale2 = ..4**-1)))


  # simulated data
  expect_equal(mexp_delayed(t = 6.5, delay1 = 5, rate1 = 0.5),
               expected = mean(pmin.int(rexp_delayed(n = 123456, delay1 = 5, rate1 = 0.5), 6.5)),
               tolerance = .003)
  expect_equal(mexp_delayed(t = 6.5, delay1 = 5, rate1 = 1.5),
               expected = mean(pmin.int(rexp_delayed(n = 123456, delay1 = 5, rate1 = 1.5), 6.5)),
               tolerance = .003)
  expect_equal(mweib_delayed(t = 6.5, delay1 = 5, shape1 = 0.5, scale1 = 2),
               expected = mean(pmin.int(rweib_delayed(n = 123456, delay1 = 5, shape1 = 0.5, scale1 = 2), 6.5)),
               tolerance = .003)
  expect_equal(mweib_delayed(t = 6.5, delay1 = 5, shape1 = 2.1, scale1 = 2),
               expected = mean(pmin.int(rweib_delayed(n = 123456, delay = 5, shape1 = 2.1, scale1 = 2), 6.5)),
               tolerance = .003)
  # 2-phase
  expect_equal(mexp_delayed(t = 8.5, delay1 = 5, rate1 = 0.5, delay2 = 7, rate2 = .1),
               expected = mean(pmin.int(rexp_delayed(n = 123456, delay1 = 5, rate1 = 0.5, delay2 = 7, rate2 = .1), 8.5)),
               tolerance = .003)
  expect_equal(mexp_delayed(t = 11.5, delay1 = 5, rate1 = 0.35, delay2 = 7, rate2 = .1),
               expected = mean(pmin.int(rexp_delayed(n = 123456, delay1 = 5, rate1 = 0.35, delay2 = 7, rate2 = .1), 11.5)),
               tolerance = .003)
  expect_equal(mweib_delayed(t = 6.5, delay1 = 5, shape1 = 0.5, scale1 = 2, delay2 = 7, shape2 = .91, scale2 = 1.1),
               expected = mean(pmin.int(rweib_delayed(n=123456, delay1 = 5, shape1 = 0.5, scale1 = 2, delay2 = 7, shape2 = .91, scale2 = 1.1), 6.5)),
               tolerance = .003)
  expect_equal(mweib_delayed(t = 8, delay1 = 5, shape1 = 0.5, scale1 = 2, delay2 = 7, shape2 = .91, scale2 = 1.1),
               expected = mean(pmin.int(rweib_delayed(n=123456, delay1 = 5, shape1 = 0.5, scale1 = 2, delay2 = 7, shape2 = .91, scale2 = 1.1), 8)),
               tolerance = .003)
})

