# mkuhn, 2021-04-07
# examples for MPS-fitting

library('purrr', warn.conflicts = FALSE)
library('tidyr', warn.conflicts = FALSE)
library('dplyr', warn.conflicts = FALSE)
library('future')
library('future.apply')
library('future.callr')

test_that("Fit delayed Exponentials", {
  # from a call to 9 + rexp(13, rate = 0.5)
  xx <- c(9.37584220431745, 9.43687826953828, 9.44079428166151, 9.63324003852904,
    9.76594032067806, 9.80526794679463, 9.90732720028609, 10.3573373407125,
    10.596041733315, 10.6229753642434, 11.1074129025543, 11.5750403608287,
    16.3800762999327)

  fd_exp <- delay_model(xx, distribution = "expon")
  coef_exp <- coef(fd_exp)

  # optim does converge properly now for this data vector!
  expect_identical(fd_exp$convergence, expected = 0L) # used to be 52L

  expect_equal(coef_exp[1L], expected = c(delay = 9), tolerance = .03)
  expect_equal(coef_exp[2L], expected = c(rate = 0.5), tolerance = .3)

  set.seed(20210429)
  yy <- rexp_delayed(21L, delay = 10.5, rate = .9)

  fd_exp2 <- delay_model(x = xx, y = yy, distribution = "exp")
  coef_exp2 <- coef(fd_exp2)

  expect_identical(fd_exp2[["convergence"]], expected = 0L)
  expect_equal(as.numeric(coef_exp2[1:2]), expected = as.numeric(coef_exp), tolerance = .002)


  # bind delay
  fd_exp2b <- delay_model(x = xx, y = yy, distribution = "exp", bind = "delay")
  coef_exp2b <- coef(fd_exp2b)

  expect_identical(fd_exp2b[["convergence"]], expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(as.numeric(coef_exp2b[1L]),
               expected = min(coef_exp2[grepl(pattern = "delay", names(coef_exp2), fixed = TRUE)]),
               tolerance = .01)

  # bind delay + rate
  fd_exp2c <- delay_model(x = xx, y = yy, distribution = "exp", bind = c("delay", "rate"))
  coef_exp2c <- coef(fd_exp2c)

  fd_exp3 <- delay_model(x = c(xx, yy), distribution = "exp")
  coef_exp3 <- coef(fd_exp3)

  expect_identical(length(coef_exp3), length(coef_exp2c))
  # similar coefficients (but the criterion is not identical, weighted mean for both groups vs overall mean)
  expect_equal(coef_exp2c[1L], coef_exp3[1L], tolerance = .01)
  expect_equal(coef_exp2c[2L], coef_exp3[2L], tolerance = .05)
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

  testres_P_H0 <- future.apply::future_vapply(X = seq(12), FUN.VALUE = double(1),
                                              FUN = function(dummy) {
    x <- 4 + rexp(13, rate = .07)
    y <- 4 + rexp(11, rate = .07)

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


test_that("Fit delayed Weibull", {

  # from Dumonceaux and Antle (1973)
  # as cited by Cheng (1982)
  maxFloodLvl <- c(
    0.654, 0.613, 0.315, 0.449, 0.297,
    0.402, 0.379, 0.423, 0.379, 0.3235,
    0.269, 0.740, 0.418, 0.412, 0.494,
    0.416, 0.338, 0.392, 0.484, 0.265
  )
  fd_maxFl <- incubate::delay_model(maxFloodLvl, distribution = "weib")
  coef_maxFl <- coef(fd_maxFl)

  expect_identical(fd_maxFl$convergence, expected = 0L)
  expect_equal(coef_maxFl[1L], expected = c(delay=0.244), tolerance = .01)
  expect_equal(coef_maxFl[2L], expected = c(shape=1.310), tolerance = .01)
  expect_equal(coef_maxFl[3L], expected = c(scale=.202),  tolerance = .01)


  # beach pollution
  # from Steen & Stickler (1976)
  # as cited by Cheng (1982)
  pollution <- c(
    1364, 2154, 2236, 2518, 2527,
    2600, 3009, 3045, 4109, 5500,
    5800, 7200, 8400, 8400, 8900,
    11500, 12700, 15300, 18300, 20400
  )

  fd_poll <- incubate::delay_model(pollution, distribution = "wei")
  objFun_poll <- fd_poll$objFun
  coef_poll <- coef(fd_poll)

  expect_identical(fd_poll$convergence, expected = 0L)
  expect_identical(length(coef_poll), expected = 3L)
  expect_lt(objFun_poll(coef_poll), expected = 3.75)
  expect_equal(coef_poll, expected = c(delay=1085, shape=0.95, scale=6562), tolerance = .001)

  # fit in two groups
  fd_wb2 <- incubate::delay_model(x = maxFloodLvl, y = pollution, distribution = "weib")
  coef_wb2 <- coef(fd_wb2)

  expect_identical(fd_wb2[["convergence"]], expected = 0L)
  expect_identical(length(coef_wb2), expected = 2L*3L)
  expect_equal(as.numeric(coef_wb2[1L]), expected = as.numeric(coef_maxFl[1L]), tolerance = .001)
  expect_equal(as.numeric(coef_wb2[2L]), expected = as.numeric(coef_maxFl[2L]), tolerance = .001)
  expect_equal(as.numeric(coef_wb2[3L]), expected = as.numeric(coef_maxFl[3L]), tolerance = .001)
  expect_equal(as.numeric(coef_wb2[4L]), expected = as.numeric(coef_poll[1L]), tolerance = .001)
  expect_equal(as.numeric(coef_wb2[5L]), expected = as.numeric(coef_poll[2L]), tolerance = .001)
  expect_equal(as.numeric(coef_wb2[6L]), expected = as.numeric(coef_poll[3L]), tolerance = .001)

  # fit in two groups with binding
  set.seed(20210430)
  rweib <- getDist("weib", type = "ran")
  fd_wb2b <- incubate::delay_model(x = rweib(n=37, delay = 7, shape = .8, scale = 3),
                               y = rweib(n=51, delay = 5, shape = 1.2, scale = 1.5),
                               distribution = "weib", bind = "delay")

  coef_wb2b <- coef(fd_wb2b)

  expect_identical(fd_wb2b[["convergence"]], expected = 0L)
  expect_identical(length(coef_wb2b), expected = 2L*3L-1L) #delay is bound
  # expect a delay close to the minimum of the two true delay parameters
  expect_equal(coef_wb2b[1L], expected = c(delay = 5), tolerance = .01)
})

