# mkuhn, 2021-04-07
# examples for MPS-fitting

library('purrr', warn.conflicts = FALSE)
library('tidyr', warn.conflicts = FALSE)
library('dplyr', warn.conflicts = FALSE)
library('future')
library('future.apply')
library('future.callr')

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
  delayTimes <- c(0, 1, 4, 7, 11, 13, rpois(n=length(rateVals) - 6L, lambda = 9))
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(dexp_delayed(x = tPoints, delay = .x, rate = .y),
                                  dweib_delayed(x = tPoints, delay = .x, shape = 1L, scale = .y**-1)))
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
               .f = ~expect_equal(dexp_delayed(x = tPoints, delay = .x, rate = .y),
                                  dweib_delayed(x = tPoints, delay = .x, shape = 1L, scale = .y**-1)))
})


test_that('GOF-test on single-group exponentials', {
  future::plan(future.callr::callr, workers = parallelly::availableCores(omit = 1L))

  # GOF-tests on true exponential data with varying sample size, delay and rate
  # results in a matrix of dimension #scenarios x #replications
  fitting_expos <- future.apply::future_replicate(n = 467L, simplify = FALSE, expr = {
    scenarios <- tidyr::expand_grid(n = c(10, 25, 50), delay = c(0, 5, 15), rate = c(.01, .2, .4, 1, 1.5, 4))
    # fit exponential models with varying n, delay and rate
    purrr::pmap(.l = scenarios,
                .f = ~ delay_model(x = rexp_delayed(n = ..1, delay = ..2, rate = ..3),
                                   distribution = 'exponential'))
  }) %>% purrr::transpose() # get a list of scenarios, each containing its models of replicated data

  # list: for each scenario, the vector of Moran's GOF-test p-value
  GOF_pvals <- list(
    moran = purrr::map(.x = fitting_expos, .f = ~ purrr::map_dbl(., \(fit) test_GOF(fit, method = 'moran')$p.value)),
    pearson = purrr::map(.x = fitting_expos, .f = ~ purrr::map_dbl(., \(fit) test_GOF(fit, method = 'pearson')$p.value)),
    ad = purrr::map(.x = fitting_expos, .f = ~ purrr::map_dbl(., \(fit) test_GOF(fit, method = 'ad')$p.value))
  )

  # to visualize the GOF P-values:
  # as_tibble(GOF_pvals[['moran']], .name_repair = 'unique') %>%
  #   pivot_longer(cols = everything()) %>%
  #   ggplot(aes(x=value, col = name)) +
  #   geom_freqpoly(binwidth = .12) + xlim(0, 1)

  # expect uniform P-values for GOF under valid H0
  # go over the three types of GOF-tests, within test type, check each scenario
  # use purrr::flatten(GOF_pvals) as data argument to walk to test *all* in one go
  purrr::walk(GOF_pvals[['pearson']], .f = ~expect_equal(object = mean(.), expected = 0.5,
                                                       tolerance = .25, info = 'pearson'))
  purrr::walk(GOF_pvals[['ad']], .f = ~expect_equal(object = mean(.), expected = 0.5,
                                                         tolerance = .25, info = 'AD'))

  # Moran's test:
  purrr::walk(GOF_pvals[['moran']], .f = ~expect_equal(object = mean(.), expected = 0.5,
                                                            tolerance = .15, info = 'moran'))
  GOF_moran_pvals_KSpval <- purrr::map_dbl(GOF_pvals[['moran']], .f = ~suppressWarnings(ks.test(., y = 'punif')$p.value))
  # not too many small P-values
  expect_lte(length(which(GOF_moran_pvals_KSpval < .05)) / length(GOF_moran_pvals_KSpval),
             expected = .3)
  # some high P-values
  expect_gte(length(which(GOF_moran_pvals_KSpval > .6)) / length(GOF_moran_pvals_KSpval),
             expected = .15)


  future::plan(sequential)
})

test_that("Test difference in delay for two exponential fits", {

  future::plan(future.callr::callr, workers = parallelly::availableCores(omit = 1L))

  set.seed(12345)
  x <- rexp_delayed(13L, delay = 11, rate = .05)
  y <- rexp_delayed(17L, delay = 11, rate = .08)

  # increasing effect
  te_diff_delays <- purrr::map(purrr::set_names(c(0, 9, 19)),
                               ~ test_diff(x = x + .x, y = y, param = "delay", R = 399))

  # null model (no effect) has a high P-value
  expect_gt(purrr::chuck(te_diff_delays, "0", "P", "boot"), expected = .1)

  # the bigger the effect (=difference in delay) the smaller the P-value
  expect_true(all(diff(purrr::map_dbl(te_diff_delays, ~ purrr::chuck(., "P", "boot"))) < 0L))
  # negative correlation betw effect size and P-values
  expect_lt(cor(x = as.integer(names(te_diff_delays)),
                y = purrr::map_dbl(te_diff_delays, ~ purrr::chuck(., "P", "boot"))),
            expected = -.66)

  #test effect of sample size: increase n and power goes up.
  set.seed(123456)
  # data with difference in delay by 2.5 time units
  #+but different sample sizes
  xs <- purrr::map(purrr::set_names(c(9, 20, 32, 37)),
      ~ 6.5 + rexp(., rate = .07))

  ys <- purrr::map(purrr::set_names(c(10, 19, 30, 38)),
            ~ 9 + rexp(., rate = .07))

  te_diff_delays_n <- purrr::map2(.x = xs, .y = ys,
       .f = ~ suppressWarnings(test_diff(x = .x, y = .y, param = "delay", R = 399)))

  expect_lt(cor(x = as.integer(names(te_diff_delays_n)),
                y = purrr::map_dbl(te_diff_delays_n, ~ purrr::chuck(., "P", "boot"))),
            expected = -.66)

  future::plan(sequential)
})



test_that("Test difference in delay when H0 is true (no difference in delay)", {

  plan(future.callr::callr, workers = parallelly::availableCores(omit = 1L))

  set.seed(20210506)

  testres_P_H0 <- future.apply::future_vapply(X = seq(12), FUN.VALUE = double(1L),
                                              FUN = function(dummy) {
    x <- rexp_delayed(13, delay = 4, rate = .07)
    y <- rexp_delayed(11, delay = 4, rate = .07)

    # return P-value of bootstrap test
    Pval <- NA_real_
    try(Pval <- purrr::chuck(test_diff(x = x, y = y, param = "delay", distribution = "exp", R = 201), "P", "boot"),
        silent = TRUE)

    Pval
  }, future.seed = TRUE)

  # drop failed tests (i.e. P=NA) because of trouble in `optim`
  testres_P_H0 <- testres_P_H0[is.finite(testres_P_H0)]

  # KS-test does not reject H0: uniform distribution
  expect_gt(suppressWarnings(stats::ks.test(x = testres_P_H0, "punif")$p.value), expected = .2)

  future::plan(sequential)

})



