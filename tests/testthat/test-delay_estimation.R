# mkuhn, 2021-04-07
# testing the delay estimation,
# in particular parameter estimates, convergence etc. from the model fit object

test_that('Parameter extraction and transformation', {
  par_exp1 <- c(delay1 = 3, rate1 = .8)

  expect_identical(incubate:::extractPars(parV = par_exp1), par_exp1)
  expect_identical(incubate:::extractPars(parV = par_exp1, group = "x"), par_exp1)
  expect_identical(incubate:::extractPars(parV = par_exp1, group = "y"), par_exp1) # group is ignored here
  expect_identical(incubate:::extractPars(parV = par_exp1, transform = TRUE), c(delay1_tr = par_exp1[[1L]], rate1_tr = log(par_exp1[[2L]])))

  par_exp2 <- c(delay1.x = 2.8, rate1.x = .81, delay1.y = 5.1, rate1.y = 1.1)
  expect_identical(incubate:::extractPars(parV = par_exp2), par_exp2)
  expect_identical(incubate:::extractPars(parV = par_exp2, group = "z"), par_exp2) # strange groups are ignored
  expect_identical(incubate:::extractPars(parV = par_exp2, group = "x"), setNames(par_exp2[1:2], c("delay1", "rate1")))
  expect_identical(incubate:::extractPars(parV = par_exp2, group = "y"), setNames(par_exp2[3:4], c("delay1", "rate1")))

  par_exp2b <- c(rate1 = .23, delay1.x = 2.8, delay1.y = 5.1)
  expect_identical(incubate:::extractPars(parV = par_exp2b), par_exp2b)
  expect_identical(incubate:::extractPars(parV = par_exp2b, group = "x"), setNames(par_exp2b[c(2, 1)], c("delay1", "rate1")))
  expect_identical(incubate:::extractPars(parV = par_exp2b, group = "y"), setNames(par_exp2b[c(3, 1)], c("delay1", "rate1")))
  expect_identical(incubate:::extractPars(parV = par_exp2b, transform = TRUE),
                   c(rate1_tr = log(par_exp2b[["rate1"]]), delay1_tr.x = par_exp2b[[2]], delay1_tr.y = par_exp2b[[3]]))
  expect_identical(incubate:::extractPars(parV = par_exp2b, group = "x", transform = TRUE),
                   c(delay1_tr = par_exp2b[[2]], rate1_tr = log(par_exp2b[[1]])))
  # idem-potent (round-trip)
  expect_identical(incubate:::extractPars(incubate:::extractPars(parV = par_exp2b, transform = TRUE), transform = TRUE), par_exp2b)


  par_weib1 <- c(delay = 5, shape = .8, scale = 1.2)
  expect_identical(incubate:::extractPars(parV = par_weib1, distribution = "weibull"), par_weib1)

  par_weib1s <- c(delay1 = 5, shape1 = .8)
  expect_identical(incubate:::extractPars(parV = par_weib1s, distribution = "weibull"), par_weib1s)
  expect_identical(incubate:::extractPars(parV = par_weib1s, distribution = "weibull", transform = TRUE),
                   c(delay1_tr = par_weib1s[["delay1"]], shape1_tr = log(par_weib1s[["shape1"]])))

  par_weib2 <- c(delay1.x = 1, shape1.x = 1.8, scale1.x = 34, delay1.y = 4.2, shape1.y = 3.4, scale1.y = 12)
  par_weib2.y <- setNames(par_weib2[-c(1:3)], nm = c("delay1", "shape1", "scale1"))
  expect_identical(incubate:::extractPars(parV = par_weib2, distribution = "weibull"), par_weib2)
  expect_identical(incubate:::extractPars(par_weib2, distribution = "weibull", group = "y"),
                   par_weib2.y)
  expect_identical(incubate:::extractPars(par_weib2, distribution = "weibull", group = "y", transform = TRUE),
                   setNames(c(par_weib2.y[1], log(par_weib2.y[-1])), nm = paste0(names(par_weib2.y), "_tr")))

  par_weib2s <- par_weib2[-c(3, 6)] # drop scale
  expect_identical(incubate:::extractPars(par_weib2s, distribution = "weibull", group = "y"),
                   setNames(par_weib2s[c(3,4)], nm = c("delay1", "shape1")))

  expect_identical(incubate:::extractPars(par_weib2s, distribution = "weib", group = "y", transform = TRUE),
                   setNames(c(par_weib2s[3], log(par_weib2s[4])), nm = c("delay1_tr", "shape1_tr")))

})


test_that("Fit delayed Exponentials", {

  # single group ------------------------------------------------------------

  # from a call to rexp_delayed(16, delay = 9, rate = 0.5)
  xx <- c(9.3758422, 9.436878, 9.440794, 9.63324, 9.6421,
          9.76594, 9.80527, 9.907327, 10.357, 10.371,
          10.596, 10.623, 11.1074, 11.575, 11.849, 16.38)

  fd_exp <- delay_model(xx, distribution = "expon")
  coef_exp <- coef(fd_exp)

  expect_type(fd_exp$data, type = 'double')
  expect_identical(length(fd_exp$data), expected = 16L)
  expect_named(coef(fd_exp), expected = c('delay1', 'rate1'))
  expect_named(fd_exp$optimizer, expected = c('parOpt', "valOpt", 'profiled', 'convergence', 'message', 'counts', 'optim_args'))
  # optim converges properly for this data vector!
  # * convergence=51 is warning from L-BFGS-B
  # * convergence=52 is error from L-BFGS-B
  expect_identical(purrr::chuck(fd_exp, "optimizer", "valOpt"), fd_exp[["criterion"]]) # same for MSPE
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
  expect_gte(fd_exp_updW$criterion + 1e-05, expected = fd_exp$criterion)

  # MLE fits -----------------------------------------------------------
  fd_exp_MLEn <- delay_model(xx, distribution = 'expon', method = 'MLEn')

  expect_type(fd_exp_MLEn$data, type = 'double')
  expect_identical(length(fd_exp_MLEn$data), expected = 16L)
  expect_identical(length(coef(fd_exp_MLEn)), expected = 2L)

  # MLEn uses analytical solution when in single group mode
  expect_named(fd_exp_MLEn$optimizer, expected = c("parOpt", "valOpt", 'profiled', "convergence", "message", "counts"))
  expect_identical(purrr::chuck(fd_exp_MLEn, 'optimizer', 'convergence'), expected = 0L)
  fd_exp_MLEn_opt <- attr(fd_exp_MLEn$objFun, which = 'opt', exact = TRUE)
  expect_type(fd_exp_MLEn_opt, type = 'list')
  expect_named(fd_exp_MLEn_opt, expected = c('par', 'value', 'convergence', 'message', 'counts'))

  expect_identical(coef(fd_exp_MLEn)[['delay1']], expected = xx[[1L]])
  # MLE's later delay is compensated for by higher estimated rate
  purrr::walk(seq_along(coef_exp), .f = ~ expect_gt(coef(fd_exp_MLEn)[[.x]], expected = coef_exp[[.x]]))
  expect_equal(coef(fd_exp_MLEn)[['rate1']], expected = (mean(xx) - min(xx))**-1)


  fd_exp_MLEc <- delay_model(xx, distribution = 'expon', method = 'MLEc')

  expect_type(fd_exp_MLEc$data, type = 'double')
  expect_identical(length(fd_exp_MLEc$data), expected = 16L)
  expect_identical(length(coef(fd_exp_MLEc)), expected = 2L)

  expect_named(fd_exp_MLEc$optimizer, expected = c("parOpt", "valOpt", "profiled", "convergence", "message", "counts", "optim_args"))
  expect_identical(purrr::chuck(fd_exp_MLEc, 'optimizer', 'convergence'), expected = 0L)
  # no exact solution for MLEc
  expect_null(attr(fd_exp_MLEc$objFun, which = 'opt', exact = TRUE))


  # check objective function
  # objective function is for transformed parameters:
  purrr::walk2(.x = runif(n=7, min=0.1, max=9),   #delay1
               .y = runif(n=7, min=0.001, max=2), #rate1
               .f = ~ expect_equal(- length(xx) * (log(.y) - .y * (mean(xx) - .x)),
                                   expected = fd_exp_MLEn$objFun(pars = incubate:::extractPars(parV = c(delay1=.x, rate1=.y), distribution = "expon", transform = TRUE)))
  )


  # 2nd group ---------------------------------------------------------------

  set.seed(20220429)
  yy <- rexp_delayed(27L, delay1 = 10.2, rate1 = .9)

  fd_exp2 <- delay_model(x = xx, y = yy, distribution = "expon")
  coef_exp2 <- coef(fd_exp2)

  expect_type(fd_exp2$data, type = 'list')
  # data gets sorted
  expect_identical(fd_exp2$data$y, sort(yy))
  expect_named(fd_exp2$optimizer, expected = c("parOpt", "valOpt", "profiled", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2, 'optimizer', 'convergence'), expected = 0L)
  expect_type(purrr::chuck(fd_exp2, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef_exp2[1:2]), expected = as.numeric(coef_exp), tolerance = .01)

  # no delay in 2nd group
  set.seed(20221124)
  zz <- rexp_delayed(13L, delay1=0, rate1=1.1)
  fd_exp2nd <- delay_model(x = xx, y = zz, distribution = "expon") #nd = no delay
  coef_exp2nd <- coef(fd_exp2nd)

  expect_identical(purrr::chuck(fd_exp2nd, "optimizer", "convergence"), expected = 0L)
  expect_equal(coef_exp2nd[[3]], 0) # no delay in 3rd group

  # MLE fits -----
  fd_exp2_MLEn <- delay_model(x = xx, y = yy, distribution = "expon", method = "MLEn")
  expect_type(fd_exp2_MLEn$data, type = "list")
  expect_identical(fd_exp2_MLEn$data$y, sort(yy))
  expect_named(fd_exp2_MLEn$optimizer, expected = c("parOpt", "valOpt", "profiled", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2_MLEn, 'optimizer', 'convergence'), expected = 0L) # converged
  expect_type(purrr::chuck(fd_exp2_MLEn, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef(fd_exp2_MLEn)[1:2]), expected = as.numeric(coef(fd_exp_MLEn)), tolerance = .01)
  # MLEn fits have larger (=later) delay and (to recover a higher rate)
  purrr::walk(.x = seq_along(coef(fd_exp2)), .f = ~ expect_gte(coef(fd_exp2_MLEn)[.x], coef(fd_exp2)[.x]))

  fd_exp2_MLEc <- delay_model(x = xx, y = yy, distribution = "expon", method = "MLEc")
  expect_type(fd_exp2_MLEc$data, type = "list")
  expect_identical(fd_exp2_MLEc$data$y, sort(yy))
  expect_named(fd_exp2_MLEc$optimizer, expected = c("parOpt", "valOpt", "profiled", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2_MLEc, 'optimizer', 'convergence'), expected = 0L) # converged
  expect_type(purrr::chuck(fd_exp2_MLEc, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef(fd_exp2_MLEc)[1:2]), expected = as.numeric(coef(fd_exp_MLEc)), tolerance = .01)

  # bind delay
  fd_exp2b <- delay_model(x = xx, y = yy, distribution = "expon", bind = "delay1")
  coef_exp2b <- coef(fd_exp2b)

  expect_named(coef_exp2b, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b$optimizer, expected = c("parOpt", "valOpt", "profiled", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b[[1L]],
               expected = min(coef_exp2[grepl(pattern = "delay1", names(coef_exp2), fixed = TRUE)]),
               tolerance = .01)

  # bind delay with MLEn
  fd_exp2b_MLEn <- delay_model(x = xx, y = yy, distribution = "expon", bind = "delay1", method = "MLEn")
  coef_exp2b_MLEn <- coef(fd_exp2b_MLEn)

  expect_named(coef_exp2b_MLEn, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_MLEn$optimizer, expected = c("parOpt", "valOpt", "profiled", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_MLEn, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b_MLEn[[1L]],
               expected = min(coef(fd_exp2_MLEn)[grepl(pattern = "delay1", names(coef(fd_exp2_MLEn)), fixed = TRUE)]),
               tolerance = .01)

  # bind delay with MLEc
  fd_exp2b_MLEc <- delay_model(x = xx, y = yy, distribution = "expon", bind = "delay1", method = "MLEc")
  coef_exp2b_MLEc <- coef(fd_exp2b_MLEc)

  expect_named(coef_exp2b_MLEc, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_MLEc$optimizer, expected = c("parOpt", "valOpt", "profiled", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_MLEc, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b_MLEc[[1L]],
               expected = min(coef(fd_exp2_MLEc)[grepl(pattern = "delay1", names(coef(fd_exp2_MLEc)), fixed = TRUE)]),
               tolerance = .01)


  # bind delay + rate
  fd_exp2c <- delay_model(x = xx, y = yy, distribution = "expon", bind = c("delay1", "rate1"))
  coef_exp2c <- coef(fd_exp2c)

  # data combinded
  fd_exp2comb <- delay_model(x = c(xx, yy), distribution = "expon")
  coef_exp2comb <- coef(fd_exp2comb)

  expect_named(coef_exp2c, expected = c("delay1", "rate1"))
  expect_identical(length(coef_exp2comb), length(coef_exp2c))
  # similar coefficients (but the criterion is not identical, weighted mean for both groups vs overall mean)
  expect_equal(coef_exp2c[[1L]], coef_exp2comb[[1L]], tolerance = .01)
  expect_equal(coef_exp2c[[2L]], coef_exp2comb[[2L]], tolerance = .05)

  # bind delay + rate for MLE
  fd_exp2c_MLEn <- delay_model(x = xx, y = yy, distribution = "expon", bind = c("delay1", "rate1"), method = "MLEn")
  coef_exp2c_MLEn <- coef(fd_exp2c_MLEn)
  fd_exp2comb_MLEn <- delay_model(x = c(xx, yy), distribution = "expon", method = "MLEn")
  coef_exp2comb_MLEn <- coef(fd_exp2comb_MLEn)

  expect_named(coef_exp2c_MLEn, expected = c("delay1","rate1"))
  expect_identical(length(coef_exp2comb_MLEn), length(coef_exp2c_MLEn))
  purrr::walk(seq_along(coef_exp2c_MLEn), ~expect_equal(coef_exp2c_MLEn[[.x]], coef_exp2comb_MLEn[[.x]], tolerance = .05))

  # MLEc
  fd_exp2c_MLEc <- delay_model(x = xx, y = yy, distribution = "expon", bind = c("delay1", "rate1"), method = "MLEc")
  coef_exp2c_MLEc <- coef(fd_exp2c_MLEc)
  fd_exp2comb_MLEc <- delay_model(x = c(xx, yy), distribution = "expon", method = "MLEc")
  coef_exp2comb_MLEc <- coef(fd_exp2comb_MLEc)

  expect_named(coef_exp2c_MLEc, expected = c("delay1","rate1"))
  expect_identical(length(coef_exp2comb_MLEc), length(coef_exp2c_MLEc))
  purrr::walk(seq_along(coef_exp2c_MLEc), ~expect_equal(coef_exp2c_MLEc[[.x]], coef_exp2comb_MLEc[[.x]], tolerance = .05))

})



test_that("Fit delayed Weibull", {

  # susquehanna is an example dataset within incubate
  fd_maxFl <- delay_model(susquehanna, distribution = "weib")
  coef_maxFl <- coef(fd_maxFl)

  expect_identical(purrr::chuck(fd_maxFl, 'optimizer', 'convergence'), expected = 0L)
  expect_equal(coef_maxFl, expected = c(delay1=0.244, shape1=1.310, scale1=.202), tolerance = .005)

  # MLEn fit to susquehanna
  fd_maxFl_MLEn <- delay_model(susquehanna, distribution = "weib", method = "MLEn")
  coef_maxFl_MLEn <- coef(fd_maxFl_MLEn)
  expect_identical(purrr::chuck(fd_maxFl_MLEn, "optimizer", "convergence"), expected = 0L)
  expect_gte(coef_maxFl_MLEn[[1]], coef_maxFl[[1]])

  # check the objective function
  purrr::pwalk(.l = list(delay1=runif(7, max=0.25), shape1=runif(7, max=3), scale1=runif(7, max=5)),
              .f = ~ {
                xc <- susquehanna - ..1
                expect_equal(-length(susquehanna) * (log(..2) - ..2 * log(..3) + (..2 - 1L) * mean(log(xc)) - mean(xc**..2)/..3**..2),
                             expected = fd_maxFl_MLEn$objFun(incubate:::extractPars(c(delay1=..1, shape1=..2, scale1=..3), distribution = 'weibull', transform = TRUE)))

                })

  fd_maxFl_MLEnp <- delay_model(susquehanna, distribution = "weibu", method = "MLEn", profile = TRUE)

  # MLEn profiling by hand:
  llFun <- function(theta){
    stopifnot( is.numeric(theta), length(theta) == 2L )
    a <- theta[[1L]]
    k <- theta[[2L]]
    (1/k + mean(log(susquehanna-a)) - sum(log(susquehanna - a) * (susquehanna - a)**k)/sum((susquehanna-a)**k))^2 +
      # first factor inverse of harmonic mean of (susquehanna - a)
      (mean(1/(susquehanna-a)) * sum((susquehanna-a)**k)/sum((susquehanna-a)**(k-1)) - k/(k-1))^2
  }

  # the objective function of model coincides with the manual objective function!
  expect_identical(fd_maxFl_MLEnp$objFun(fd_maxFl_MLEnp$optimizer$parOpt), llFun(coef(fd_maxFl_MLEnp)[1:2]))

  # manual optimization
  opt_maxFl_MLEnp_man <- optim(par = c(a=0.255, k=1.8), #c(a=0.165, k=exp(1.3847)), #c(a=0.25, k=1.5),
                               fn = llFun, method = "L-BFGS-B",
                      lower = c(0, 1 + 1.49e-8), upper = c(min(susquehanna)-1.49e-8, +Inf))
  coef_maxFl_MLEnp_man <- purrr::set_names(c(opt_maxFl_MLEnp_man$par, mean((susquehanna-opt_maxFl_MLEnp_man$par["a"])^opt_maxFl_MLEnp_man$par["k"])^(1/opt_maxFl_MLEnp_man$par["k"])),
                                  nm = c("delay1", "shape1", "scale1"))

  #llFun(coef_maxFl_MLEnp_man[1:2])
  #llFun(coef(fd_maxFl_MLEnp)[1:2])

  # # visualize the optimization function landscape
  # objFunLS_maxFl_MLEnp <- tidyr::expand_grid(a = seq.int(0, .26, length.out = 17),
  #                                            k = seq.int(1.3, 3, length.out = 17)) %>%
  #   #head %>%
  #   rowwise() %>%
  #   dplyr::mutate(ll = llFun(theta = c(a,k))) %>%
  #   ungroup()
  #
  # ggplot(objFunLS_maxFl_MLEnp, aes(x = a, y = k, z = log(ll+.1), colour = after_stat(level))) +
  #   geom_contour() +
  #   scale_colour_distiller(palette = "YlGn", direction = -1)

  # estimates are roughly equal
  expect_equal(coef_maxFl_MLEnp_man, coef_maxFl_MLEn, tolerance = .25)


  # pollution is an example dataset within incubate
  fd_poll <- delay_model(pollution, distribution = "weib")
  objFun_poll <- fd_poll$objFun
  coef_poll <- coef(fd_poll)

  expect_identical(purrr::chuck(fd_poll, 'optimizer', 'convergence'), expected = 0L)
  expect_identical(length(coef_poll), expected = 3L)
  expect_lt(objFun_poll(incubate:::extractPars(parV = coef_poll, distribution = "weibull", transform = TRUE)), expected = 3.75)
  expect_equal(coef_poll, expected = c(delay1=1085, shape1=0.95, scale1=6562), tolerance = .001)


  # Cousineau's numerical example, taken from Weibull with delay = 300, shape k = 2 and scale = 100
  x <- c(310, 342, 353, 365, 383, 393, 403, 412, 451, 456)
  densFun_wb <- getDist(distribution = "weib", type = "density")

  fd_wbc_mlen <- delay_model(x = x, distribution = "weib", method = "MLEn")
  expect_equal(coef(fd_wbc_mlen), expected = c(delay1 = 274.8, shape1 = 2.80, scale1 = 126.0), tolerance = .001)
  fd_wbc_mlenp <- delay_model(x = x, distribution = "weib", method = "MLEn", profile = TRUE)
  expect_equal(coef(fd_wbc_mlenp), expected = c(delay1=280.9, shape1=2.62, scale1=119.0), tolerance = .05)
  # our implementation finds a smaller objective value (which is to be minimized)
  expect_lte(fd_wbc_mlenp$optimizer$valOpt, fd_wbc_mlenp$objFun(c(delay1=280.9, shape1=log(2.62))))
  # but log-likelihood is in fact quite similar (actually, Cousineau's solution is slightly better)
  expect_equal(fd_wbc_mlenp$criterion, -sum(densFun_wb(x, delay1 = 280.9, shape1=2.62, scale1=119.0, log = TRUE)), tolerance = .01)
  fd_wbc_mlewp <- delay_model(x = x, distribution = "weib", method = "MLEw", profile = TRUE)
  expect_true(fd_wbc_mlewp$optimizer$profiled)
  expect_equal(coef(fd_wbc_mlewp), c(delay1=283.7, shape1=2.29, scale1=116.0), tolerance = .04)
  # our implementation finds a smaller objective value.
  expect_lte(fd_wbc_mlewp$optimizer$valOpt, fd_wbc_mlewp$objFun(c(delay1=283.7, shape1=log(2.29))))
  # but log-likelihood is in fact quite similar (actually, Cousineau's solution is slightly better)
  expect_equal(fd_wbc_mlewp$criterion, -sum(densFun_wb(x, delay1=283.7, shape1=2.29, scale1=116.0, log = TRUE)), tolerance = .01)

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
