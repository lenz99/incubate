# mkuhn, 2021-04-07
# testing the delay estimation,
# in particular parameter estimates, convergence etc. from the model fit object

test_that("Fit delayed Exponentials", {

  # single group ------------------------------------------------------------

  # rexp_delayed(16, delay = 9, rate = 0.5)
  xx <- c(9.3758422, 9.436878, 9.440794, 9.63324, 9.6421,
          9.76594, 9.80527, 9.907327, 10.357, 10.371,
          10.596, 10.623, 11.1074, 11.575, 11.849, 16.38)

  fd_exp <- delay_model(xx, distribution = "expon")
  coef_exp <- coef(fd_exp)

  expect_type(fd_exp$data, type = 'double')
  expect_identical(length(fd_exp$data), expected = 16L)
  expect_named(coef(fd_exp), expected = c('delay1', 'rate1'))
  expect_named(fd_exp$optimizer, expected = c('parOpt', 'convergence', 'message', 'counts', 'optim_args'))
  # optim converges properly for this data vector!
  # * convergence=51 is warning from L-BFGS-B
  # * convergence=52 is error from L-BFGS-B
  expect_identical(purrr::chuck(fd_exp, 'optimizer', 'convergence'), expected = 0L)
  expect_type(purrr::chuck(fd_exp, 'optimizer', 'optim_args'), type = 'list')

  expect_equal(coef_exp[["delay1"]], expected = 9, tolerance = .04)
  expect_equal(coef_exp[["rate1"]], expected = 0.5, tolerance = .37)

  # update does not change structure
  fd_exp_oa <- purrr::pluck(fd_exp, 'optimizer', 'optim_args')
  expect_identical(update(fd_exp, optim_args = fd_exp_oa), expected = fd_exp)
  # effect of worse start values
  fd_exp_updW <- update(fd_exp, optim_args = purrr::assign_in(fd_exp_oa, 'par', c(1, .1)))
  expect_gt(min(fd_exp_updW$optimizer$counts), expected = min(fd_exp$optimizer$counts))
  # objective function does not reach a much better (=lower) value when starting with worse start value
  # (add a little safety margin as buffer)
  expect_gte(fd_exp_updW$val + 1e-05, expected = fd_exp$val)

  # MLE0 fit ----------------------------------------------------------------
  fd_exp_MLE0 <- delay_model(xx, distribution = 'expon', method = 'MLEn')

  expect_type(fd_exp_MLE0$data, type = 'double')
  expect_identical(length(fd_exp_MLE0$data), expected = 16L)
  expect_identical(length(coef(fd_exp_MLE0)), expected = 2L)

  expect_identical(coef(fd_exp_MLE0)[['delay1']], expected = xx[[1L]])
  # MLE's later delay is compensated for by higher estimated rate
  expect_gt(coef(fd_exp_MLE0)[['rate1']], expected = coef_exp[['rate1']])

  expect_named(fd_exp_MLE0$optimizer, expected = c("parOpt", "convergence", "message", "counts"))
  expect_identical(purrr::chuck(fd_exp_MLE0, 'optimizer', 'convergence'), expected = 0L)
  fd_exp_MLE0_opt <- attr(fd_exp_MLE0$objFun, which = 'opt', exact = TRUE)
  expect_type(fd_exp_MLE0_opt, type = 'list')
  expect_named(fd_exp_MLE0_opt, expected = c('par', 'value', 'convergence', 'message', 'counts'))


  # 2nd group ---------------------------------------------------------------

  set.seed(20220429)
  yy <- rexp_delayed(27L, delay1 = 10.2, rate1 = .9)

  fd_exp2 <- delay_model(x = xx, y = yy, distribution = "expon")
  coef_exp2 <- coef(fd_exp2)

  expect_type(fd_exp2$data, type = 'list')
  # data gets sorted
  expect_identical(fd_exp2$data$y, sort(yy))
  expect_named(fd_exp2$optimizer, expected = c("parOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2, 'optimizer', 'convergence'), expected = 0L)
  expect_type(purrr::chuck(fd_exp2, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef_exp2[1:2]), expected = as.numeric(coef_exp), tolerance = .01)


  # bind delay
  fd_exp2b <- delay_model(x = xx, y = yy, distribution = "expon", bind = "delay1")
  coef_exp2b <- coef(fd_exp2b)

  expect_named(fd_exp2b$optimizer, expected = c("parOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b[[1L]],
               expected = min(coef_exp2[grepl(pattern = "delay1", names(coef_exp2), fixed = TRUE)]),
               tolerance = .01)

  # bind delay + rate
  fd_exp2c <- delay_model(x = xx, y = yy, distribution = "expon", bind = c("delay1", "rate1"))
  coef_exp2c <- coef(fd_exp2c)

  fd_exp3 <- delay_model(x = c(xx, yy), distribution = "expon")
  coef_exp3 <- coef(fd_exp3)

  expect_identical(length(coef_exp3), length(coef_exp2c))
  # similar coefficients (but the criterion is not identical, weighted mean for both groups vs overall mean)
  expect_equal(coef_exp2c[1L], coef_exp3[1L], tolerance = .01)
  expect_equal(coef_exp2c[2L], coef_exp3[2L], tolerance = .05)
})



test_that("Fit delayed Weibull with MPSE", {

  fd_maxFl <- delay_model(susquehanna, distribution = "weib")
  coef_maxFl <- coef(fd_maxFl)

  expect_identical(purrr::chuck(fd_maxFl, 'optimizer', 'convergence'), expected = 0L)
  expect_equal(coef_maxFl, expected = c(delay1=0.244, shape1=1.310, scale1=.202), tolerance = .005)

  fd_poll <- delay_model(pollution, distribution = "weib")
  objFun_poll <- fd_poll$objFun
  coef_poll <- coef(fd_poll)

  expect_identical(purrr::chuck(fd_poll, 'optimizer', 'convergence'), expected = 0L)
  expect_identical(length(coef_poll), expected = 3L)
  expect_lt(objFun_poll(incubate:::extractPars(parV = coef_poll, distribution = "weibull", transform = TRUE)), expected = 3.75)
  expect_equal(coef_poll, expected = c(delay1=1085, shape1=0.95, scale1=6562), tolerance = .001)

  # fit in two groups
  fd_wb2 <- incubate::delay_model(x = susquehanna, y = pollution, distribution = "weib")
  coef_wb2 <- coef(fd_wb2)

  expect_identical(purrr::chuck(fd_wb2, 'optimizer', 'convergence'), expected = 0L)
  expect_identical(length(coef_wb2), expected = 2L*3L)
  expect_equal(as.numeric(coef_wb2[1:3]), expected = as.numeric(coef_maxFl), tolerance = .001)
  # there might occur quite some differences in the coefficients
  expect_equal(as.numeric(coef_wb2[4:6]), expected = as.numeric(coef_poll), tolerance = .25)

  # fit model for two groups with binding
  set.seed(20210430)
  fd_wb2b <- delay_model(x = rweib_delayed(n=37, delay1 = 7, shape1 = .8, scale1 = 3),
                         y = rweib_delayed(n=51, delay1 = 5, shape1 = 1.2, scale1 = 1.5),
                         distribution = "weib", bind = "delay1")

  coef_wb2b <- coef(fd_wb2b)

  expect_identical(purrr::chuck(fd_wb2b, 'optimizer', 'convergence'), expected = 0L)
  expect_identical(length(coef_wb2b), expected = 2L*3L-1L) #delay is bound
  # expect a delay close to the minimum of the two true delay parameters
  expect_equal(coef_wb2b[1L], expected = c(delay1 = 5), tolerance = .02)
})


test_that('Confidence intervals', {
  testthat::skip_on_cran()
  set.seed(1234)


  obs_sim <- rexp_delayed(n = 29L, delay = 5, rate = .3)

  fm <- delay_model(x = obs_sim, distribution = 'expon')

  # use default smd_factor
  bs_data <- incubate:::bsDataStep(fm, R = 1199, useBoot = FALSE)
  # set up identical bs_data from boot-package
  bs_data_bt <- incubate:::bsDataStep(fm, R = 3, useBoot = TRUE)
  bs_data_bt$R <- NCOL(bs_data)
  bs_data_bt$t <- t(bs_data)

  # agreement of own implementation and boot-implementation
  purrr::walk(.x = c('normal', 'lognormal', 'quantile', 'logquantile', 'quantile0'),
              .f = ~ {
                ci_own <- confint(fm, bs_data = bs_data, bs_infer = .x)
                ci_boot <- confint(fm, bs_data = bs_data_bt, bs_infer = .x, useBoot = TRUE)
                expect_equal(ci_own, ci_boot, tolerance = 1e-3)})

})
