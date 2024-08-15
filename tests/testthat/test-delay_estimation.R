# mkuhn, 2021-04-07
# testing the delay estimation,
# in particular parameter estimates, convergence etc. from the model fit object

test_that('Parameter extraction and transformation', {

  # test getDist for parameter names:
  # on transformed scale (for optimization)
  expect_identical(incubate:::getDist(distribution = "weibull", type = "param", twoPhase = FALSE, twoGroup = TRUE,
                                      bind = "shape1", profiled = FALSE, transformed = TRUE),
                   expected = c("shape1_tr", "delay1_tr.x", "scale1_tr.x", "delay1_tr.y", "scale1_tr.y"))
  expect_identical(incubate:::getDist(distribution = "weibull", type = "param", twoPhase = FALSE, twoGroup = TRUE,
                                      bind = "shape1", profiled = TRUE, transformed = TRUE),
                   expected = c("shape1_tr", "delay1_tr.x", "delay1_tr.y"))
  # original (=non-transformed) scale
  expect_identical(incubate:::getDist(distribution = "weibull", type = "param", twoPhase = FALSE, twoGroup = TRUE,
                                      bind = "shape1", profiled = FALSE, transformed = FALSE),
                   expected = c("shape1", "delay1.x", "scale1.x", "delay1.y", "scale1.y"))
  # getDist takes out profiled parameters even when for original (=non-transformed) scale
  #+ this is needed in delay_estimation (at extractParOptInd)
  expect_identical(incubate:::getDist(distribution = "weibull", type = "param", twoPhase = FALSE, twoGroup = TRUE,
                                      bind = "shape1", profiled = TRUE, transformed = FALSE),
                   expected = c("shape1", "delay1.x", "delay1.y"))

  #' Slow extract-routine based on R's name-matching capabilities as gold standard function.
  #'
  #' Extract certain parameters from a given named parameter vector. External, slow variant that parses parameter names.
  #'
  #' This routine handles parameter vectors for one or two groups. To this end, it parses the parameter names.
  #' Bound parameters (cf. `bind=`) are deduced from the naming: they lack the .x or .y suffix.
  #' it can also transform the given parameters from the standard form (as used in the delayed distribution functions)
  #' to the internal transformed parametrization as used by the optimization function, and vice versa.
  #'
  #' When parameters for a single group are requested the parameters are always returned in the canonical order of the distribution.
  #' Otherwise, the order of parameters is unchanged.
  #' It parses the parameter names to find out if it is twoPhase, twoGroup and, when `transform=TRUE`, which direction to transform.
  #'
  #' @param parV numeric named parameter vector (for single group or two groups)
  #' @param distribution character name. Exponential or Weibull
  #' @param twoPhase logical. Do we model two phases?
  #' @param group character. Extract only parameters for specified group. Recognizes 'x' and 'y'. Other values are ignored.
  #' @param isTransformed logical. Are the parameters transformed? If `NULL` (default) it is deduced from the parameter names.
  #' @param transform logical. Transform the parameters!?
  #' @return requested parameters as named numeric vector
  extractParsTest <- function(parV, distribution = c('exponential', 'weibull'), twoPhase = NULL, group = NULL, isTransformed = NULL, transform = FALSE) {

    stopifnot( is.numeric(parV), all(nzchar(names(parV))) )
    distribution <- match.arg(distribution)

    # contract: parV is canonically ordered (see below as well)
    # extractPars will not check if it needs first to reorder the parameters given by checking the names.
    pNames <- names(parV)

    # isTransformed when all parameter names contain '_tr'
    if (is.null(isTransformed)) {
      trNames <- grepl(pattern = "_tr", pNames, fixed = TRUE)
      stopifnot( sum(trNames) %in% c(0L, length(parV)) ) #all or nothing
      isTransformed <- all(trNames)
    }
    if (is.null(twoPhase)) {
      twoPhase <- any(grepl(pattern = "delay2", pNames, fixed = TRUE))
    }
    # parameter names of single group where we start from!
    oNames <- incubate:::getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, transformed = isTransformed)

    # index of parameters that are *not* group-specific
    idx.nongrp <- which(!grepl(pattern = paste0("[.][xy]$"), pNames))
    #Is there any .x or .y parameter name?
    #+when two groups but all parameters bound then isTwoGroup is FALSE
    isTwoGroup <- length(idx.nongrp) < length(parV)


    # transform parameters if requested
    # if no transformation required, simply extract the relevant parameters
    if ( transform ){

      # transform parameter vector for a single group. Also transforms parameter names.
      # @param parV1: parameter vector for a single group
      # @return transformed parameter vector
      transformPars1Test <- function(parV1){
        # parameter transformation matrices
        PARAM_TRANSF_M <- switch(distribution,
                                 exponential = matrix(c( 1, 0, 0, 0,
                                                         0, 1, 0, 0,
                                                         -1, 0, 1, 0,
                                                         0, 0, 0, 1), nrow = 4L, byrow = TRUE,
                                                      dimnames = list(c("delay1_tr", "rate1_tr", "delay2_tr", "rate2_tr"))),
                                 weibull = matrix(c( 1, 0, 0, 0, 0, 0,
                                                     0, 1, 0, 0, 0, 0,
                                                     0, 0, 1, 0, 0, 0,
                                                     -1, 0, 0, 1, 0, 0,
                                                     0, 0, 0, 0, 1, 0,
                                                     0, 0, 0, 0, 0, 1), nrow = 6L, byrow = TRUE,
                                                  dimnames = list(paste0(c("delay1", "shape1", "scale1", "delay2", "shape2", "scale2"), "_tr"))),
                                 stop("Unknown distribution!", call. = FALSE)
        )
        PARAM_TRANSF_Minv <- switch(distribution,
                                    exponential = matrix(c(1, 0, 0, 0,
                                                           0, 1, 0, 0,
                                                           1, 0, 1, 0,
                                                           0, 0, 0, 1), nrow = 4L, byrow = TRUE,
                                                         dimnames = list(c("delay1", "rate1", "delay2", "rate2"))),
                                    weibull = matrix(c(1, 0, 0, 0, 0, 0,
                                                       0, 1, 0, 0, 0, 0,
                                                       0, 0, 1, 0, 0, 0,
                                                       1, 0, 0, 1, 0, 0,
                                                       0, 0, 0, 0, 1, 0,
                                                       0, 0, 0, 0, 0, 1), nrow = 6L, byrow = TRUE,
                                                     dimnames = list(c("delay1", "shape1", "scale1", "delay2", "shape2", "scale2"))),
                                    stop("Unknown distribution", call. = FALSE)
        )


        PARAM_TRANSF_F <- list(exponential = c(identity, log, log, log),
                               weibull = c(identity, log, #log1p, #identity, #=shape1
                                           log, log, log, log))[[distribution]]
        PARAM_TRANSF_Finv <- list(exponential = c(identity, exp, exp, exp),
                                  weibull = c(identity, exp, #expm1, #identity, #=shape1
                                              exp, exp, exp, exp))[[distribution]]

        stopifnot( length(parV1) <= NCOL(PARAM_TRANSF_M), NCOL(PARAM_TRANSF_M) == length(PARAM_TRANSF_F) )


        if (isTransformed) {
          # b = Ainv %*% Finv(b')
          rlang::set_names(
            as.numeric(PARAM_TRANSF_Minv[seq_along(parV1), seq_along(parV1)] %*%
                         as.numeric(.mapply(FUN = function(f, x) f(x),
                                            dots = list(PARAM_TRANSF_Finv[seq_along(parV1)], parV1),
                                            MoreArgs = NULL))),
            nm = rownames(PARAM_TRANSF_Minv)[seq_along(parV1)])
        } else {
          # b' = F(A %*% b)
          rlang::set_names(
            as.numeric(.mapply(FUN = function(f, x) f(x),
                               dots = list(PARAM_TRANSF_F[seq_along(parV1)],
                                           as.numeric(PARAM_TRANSF_M[seq_along(parV1), seq_along(parV1)] %*% parV1)),
                               MoreArgs = NULL)),
            nm = rownames(PARAM_TRANSF_M)[seq_along(parV1)])
        }
      }# fn transformPars1Test


      # update parameter vector according to transformation
      parV <- if (! isTwoGroup) {
        transformPars1Test(parV1 = parV)
      } else {
        # **twoGroup** case (effectively)
        # target parameter names
        pNames_trgt <- if (isTransformed) {
          sub(pattern = "_tr", replacement = "", pNames, fixed = TRUE) # remove '_tr' string
        } else {
          purrr::map_chr(.x = strsplit(pNames, split = ".", fixed = TRUE),
                         .f = function(s) if (length(s) == 1L) paste0(s, "_tr") else paste0(s[[1L]], "_tr.", s[[2L]]))
        }

        # list of transformed parameters, for both groups x and y
        h <- list(
          x = { idx.grp <- which(endsWith(x = pNames, suffix = ".x"))
          # drop group naming: last two letters
          pNms_grp <- substr(pNames[idx.grp], start = 1L, stop = nchar(pNames[idx.grp], type = "chars") - 2L)
          # apply transform on 1-group parameter vector in canonical order!
          transformPars1Test(parV1 = c(parV[idx.nongrp], rlang::set_names(parV[idx.grp], nm = pNms_grp))[intersect(oNames, c(pNames, pNms_grp))]) },
          y = { idx.grp <- which(endsWith(x = pNames, suffix = ".y"))
          # drop group naming: last two letters
          pNms_grp <- substr(pNames[idx.grp], start = 1L, stop = nchar(pNames[idx.grp], type = "chars") - 2L)
          # apply transform on 1-group parameter vector in canonical order!
          transformPars1Test(parV1 = c(parV[idx.nongrp], rlang::set_names(parV[idx.grp], nm = pNms_grp))[intersect(oNames, c(pNames, pNms_grp))]) }
        )

        # parV transformed for effective twoGroup-setting
        rlang::set_names(c(h[[1L]][pNames_trgt[idx.nongrp]],
                           h[[1L]][! names(h[[1L]]) %in% pNames_trgt[idx.nongrp]],
                           h[[2L]][! names(h[[2L]]) %in% pNames_trgt[idx.nongrp]]),
                         nm = pNames_trgt)
      }

      # reflect transformation in meta-data
      isTransformed <- ! isTransformed
      pNames <- names(parV)
      oNames <- incubate:::getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, transformed = isTransformed)
    }# transform



    # return parameter vector
    if ( ! isTwoGroup || is.null(group) || length(group) != 1L || ! group %in% c('x', 'y') ) {
      parV
    } else {
      # single group extraction
      stopifnot( is.character(group), nzchar(group) )

      # extract all group parameters (=contain .x or .y at the end)
      #+and restore original name (=remove ".x" or ".y")
      # contract: parV is canonically ordered
      parV.gr <- parV[grepl(pattern = paste0(".", group), x = pNames, fixed = TRUE)]
      names(parV.gr) <- substr(names(parV.gr), start = 1L, stop = nchar(names(parV.gr))-2L)
      #parV.gr <- rlang::set_names(parV.gr, nm = setdiff(oNames, pNames[idx.nongrp]))

      # restore original order
      # contract: bind is intersected (only valid names, canonical order)
      c(parV.gr, parV[pNames[idx.nongrp]])[intersect(oNames, c(pNames, names(parV.gr)))]
    }
  }

  extPars_exp1 <- incubate:::objFunFactory(x = rexp_delayed(n=11, delay1 = 5, rate1 = .2), distribution = "expon", twoPhase = FALSE, profiled = FALSE) %>%
    rlang::fn_env() %>% rlang::env_get("extractPars")

  extPars_exp1P <- incubate:::objFunFactory(x = rexp_delayed(n=11, delay1 = 5, rate1 = .2), distribution = "expon", twoPhase = FALSE, profiled = TRUE) %>%
    rlang::fn_env() %>% rlang::env_get("extractPars")

  par_exp1 <- c(delay1 = 3, rate1 = .8)

  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, named = TRUE), par_exp1)
  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, named = FALSE), as.vector(par_exp1))
  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, group = "x", named = TRUE), par_exp1)
  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, group = "y", named = TRUE), par_exp1) # group= is ignored for single group!
  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, group = "z", named = TRUE), par_exp1) # group= is ignored for single group!
  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, transform = TRUE, named = TRUE), c(delay1_tr = par_exp1[[1L]], rate1_tr = log(par_exp1[[2L]])))
  expect_identical(extPars_exp1(par_exp1, isOpt = FALSE, transform = TRUE, named = FALSE), c(par_exp1[[1L]], log(par_exp1[[2L]])))

  # objective function with profiled=TRUE is different
  expect_identical(extPars_exp1P(par_exp1, isOpt = FALSE, named = TRUE), par_exp1)
  expect_identical(extPars_exp1P(par_exp1, isOpt = FALSE, transform = TRUE, named = TRUE), c(delay1_tr = par_exp1[[1L]]))
  expect_identical(extPars_exp1P(par_exp1, isOpt = FALSE, transform = TRUE, named = FALSE), par_exp1[[1L]])

  # my original (but slow) implementation
  expect_identical(extractParsTest(parV = par_exp1), par_exp1)
  expect_identical(extractParsTest(parV = par_exp1, group = "x"), par_exp1)
  expect_identical(extractParsTest(parV = par_exp1, group = "y"), par_exp1) # group= is ignored here
  expect_identical(extractParsTest(parV = par_exp1, transform = TRUE), c(delay1_tr = par_exp1[[1L]], rate1_tr = log(par_exp1[[2L]])))

  # exponential, two groups, unbound
  objFun_exp2 <- incubate:::objFunFactory(x = rexp_delayed(n=3, delay1 = 5, rate1 = .2), y = rexp_delayed(n=3, delay1 = 3, rate1 = .1), distribution = "expon", twoPhase = FALSE)
  extPars_exp2 <- rlang::env_get(rlang::fn_env(objFun_exp2), "extractPars")

  par_exp2 <- c(delay1.x = 2.8, rate1.x = .81, delay1.y = 5.1, rate1.y = 1.1)

  expect_identical(extPars_exp2(par_exp2, isOpt = FALSE, named = TRUE), par_exp2)
  expect_identical(extPars_exp2(par_exp2, isOpt = FALSE, named = TRUE, group = "z"), NULL) # strange groups gives NULL in two group setting
  expect_identical(extPars_exp2(par_exp2, isOpt = FALSE, named = TRUE, group = "x"), setNames(par_exp2[1:2], c("delay1", "rate1")))
  expect_identical(extPars_exp2(par_exp2, isOpt = FALSE, named = TRUE, group = "y"), setNames(par_exp2[3:4], c("delay1", "rate1")))

  # my original (but slow) implementation
  expect_identical(extractParsTest(parV = par_exp2), par_exp2)
  expect_identical(extractParsTest(parV = par_exp2, group = "z"), par_exp2) # strange groups are ignored with extractParsTest!
  expect_identical(extractParsTest(parV = par_exp2, group = "x"), setNames(par_exp2[1:2], c("delay1", "rate1")))
  expect_identical(extractParsTest(parV = par_exp2, group = "y"), setNames(par_exp2[3:4], c("delay1", "rate1")))

  # exponential, two groups, bound
  objFun_exp2b <- incubate:::objFunFactory(x = rexp_delayed(n=3, delay1 = 5, rate1 = .2), y = rexp_delayed(n=3, delay1 = 3, rate1 = .1),
                                           distribution = "expon", twoPhase = FALSE, bind = "rate1")
  extPars_exp2b <- rlang::env_get(rlang::fn_env(objFun_exp2b), "extractPars")

  par_exp2b <- c(rate1 = .23, delay1.x = 2.8, delay1.y = 5.1)
  expect_identical(extPars_exp2b(parV = par_exp2b, isOpt = FALSE, named = TRUE), par_exp2b)
  expect_identical(extPars_exp2b(parV = par_exp2b, isOpt = FALSE, named = TRUE, group = "x"), setNames(par_exp2b[c(2, 1)], c("delay1", "rate1")))
  expect_identical(extPars_exp2b(parV = par_exp2b, isOpt = FALSE, named = TRUE, group = "y"), setNames(par_exp2b[c(3, 1)], c("delay1", "rate1")))
  expect_identical(extPars_exp2b(parV = par_exp2b, isOpt = FALSE, named = TRUE, transform = TRUE),
                   c(rate1_tr = log(par_exp2b[["rate1"]]), delay1_tr.x = par_exp2b[[2]], delay1_tr.y = par_exp2b[[3]]))
  expect_identical(extPars_exp2b(parV = par_exp2b, isOpt = FALSE, named = TRUE, group = "x", transform = TRUE),
                   c(delay1_tr = par_exp2b[[2]], rate1_tr = log(par_exp2b[[1]])))
  # idem-potent (round-trip)
  expect_identical(extPars_exp2b(extPars_exp2b(par_exp2b, isOpt = FALSE, transform = TRUE),
                                 isOpt = TRUE, named = TRUE, transform = TRUE), par_exp2b)

  # my original (but slow) implementation
  expect_identical(extractParsTest(parV = par_exp2b), par_exp2b)
  expect_identical(extractParsTest(parV = par_exp2b, group = "x"), setNames(par_exp2b[c(2, 1)], c("delay1", "rate1")))
  expect_identical(extractParsTest(parV = par_exp2b, group = "y"), setNames(par_exp2b[c(3, 1)], c("delay1", "rate1")))
  expect_identical(extractParsTest(parV = par_exp2b, transform = TRUE),
                   c(rate1_tr = log(par_exp2b[["rate1"]]), delay1_tr.x = par_exp2b[[2]], delay1_tr.y = par_exp2b[[3]]))
  expect_identical(extractParsTest(parV = par_exp2b, group = "x", transform = TRUE),
                   c(delay1_tr = par_exp2b[[2]], rate1_tr = log(par_exp2b[[1]])))
  # idem-potent (round-trip)
  expect_identical(extractParsTest(extractParsTest(parV = par_exp2b, transform = TRUE), transform = TRUE), par_exp2b)


  # weibull
  objFun_weib1 <- incubate:::objFunFactory(x = rweib_delayed(n=3, delay1 = 5, shape1 = 1.2), distribution = "weib", twoPhase = FALSE)
  extPars_weib1 <- rlang::env_get(rlang::fn_env(objFun_weib1), "extractPars")

  par_weib1 <- c(delay1 = 5, shape1 = .8, scale1 = 1.2)
  expect_identical(extPars_weib1(par_weib1, isOpt = FALSE, named = TRUE), par_weib1)
  expect_identical(extractParsTest(par_weib1, distribution = "weibull"), par_weib1)

  # extractParsTest specific: incomplete parameter vector
  par_weib1s <- c(delay1 = 5, shape1 = .8)
  expect_identical(extPars_weib1(par_weib1s, isOpt = FALSE, named = TRUE), c(par_weib1s, scale1=NA)) #missing parameter gets named NA
  expect_identical(extractParsTest(par_weib1s, distribution = "weibull"), par_weib1s)
  expect_identical(extractParsTest(par_weib1s, distribution = "weibull", transform = TRUE),
                   c(delay1_tr = par_weib1s[["delay1"]], shape1_tr = log(par_weib1s[["shape1"]])))

  # weibull two groups
  objFun_weib2 <- incubate:::objFunFactory(x = rweib_delayed(n=3, delay1 = 5, shape1 = 1.2), y = rweib_delayed(n=4, delay1=3, shape1 = 2.8, scale1 = .9),
                                           distribution = "weib", twoPhase = FALSE)
  extPars_weib2 <- rlang::env_get(rlang::fn_env(objFun_weib2), "extractPars")

  par_weib2 <- c(delay1.x = 1, shape1.x = 1.8, scale1.x = 34, delay1.y = 4.2, shape1.y = 3.4, scale1.y = 12)
  par_weib2.y <- setNames(par_weib2[-c(1:3)], nm = c("delay1", "shape1", "scale1"))
  par_weib2.y.tr <- setNames(c(par_weib2.y[1], log(par_weib2.y[-1])), nm = paste0(names(par_weib2.y), "_tr"))

  expect_identical(extPars_weib2(par_weib2, isOpt = FALSE, named = TRUE), par_weib2)
  expect_identical(extPars_weib2(par_weib2, isOpt = FALSE, named = TRUE, group = "y"), par_weib2.y)
  expect_identical(extPars_weib2(par_weib2, isOpt = FALSE, named = TRUE, group = "y", transform = TRUE), par_weib2.y.tr)

  expect_identical(extractParsTest(par_weib2, distribution = "weibull"), par_weib2)
  expect_identical(extractParsTest(par_weib2, distribution = "weibull", group = "y"), par_weib2.y)
  expect_identical(extractParsTest(par_weib2, distribution = "weibull", group = "y", transform = TRUE), par_weib2.y.tr)

  # extractParsTest specific: incomplete parameter vector
  par_weib2s <- par_weib2[-c(3, 6)] # drop scale
  expect_identical(extractParsTest(par_weib2s, distribution = "weibull", group = "y"),
                   setNames(par_weib2s[c(3,4)], nm = c("delay1", "shape1")))

  expect_identical(extractParsTest(par_weib2s, distribution = "weibull", group = "y", transform = TRUE),
                   setNames(c(par_weib2s[3], log(par_weib2s[4])), nm = c("delay1_tr", "shape1_tr")))

})



test_that("Fit delayed Exponentials", {

  # single group ------------------------------------------------------------

  # from a call to rexp_delayed(16, delay = 9, rate = 0.5)
  exp_d9 <- c(9.3758422, 9.436878, 9.440794, 9.63324, 9.6421,
          9.76594, 9.80527, 9.907327, 10.357, 10.371,
          10.596, 10.623, 11.1074, 11.575, 11.849, 16.38)

  fd_exp <- delay_model(exp_d9, distribution = "expon", profiled = FALSE)
  coef_exp <- coef(fd_exp)

  expect_type(fd_exp$data, type = 'double')
  expect_identical(length(fd_exp$data), expected = 16L)
  expect_identical(fd_exp$method, expected = "MPSE")
  expect_named(coef(fd_exp), expected = c('delay1', 'rate1'))
  expect_named(fd_exp$optimizer, expected = c('parOpt', "valOpt", 'profiled', "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
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
  # we need more iterations during optimization
  expect_gt(min(fd_exp_updW$optimizer$counts), expected = min(fd_exp$optimizer$counts))
  # objective function does not reach a better (=lower) value when starting with worse start value
  # (add a little safety margin as buffer)
  expect_gte(fd_exp_updW$criterion + 1e-05, expected = fd_exp$criterion)

  # MPSE fit with profiled=TRUE
  fd_exp_P <- delay_model(exp_d9, distribution = "expon", profiled = TRUE)
  expect_equal(coef(fd_exp_P), coef_exp, tolerance = .01) # similar coefficients
  expect_named(fd_exp_P$optimizer$parOpt, expected = "delay1_tr")

  # another data set
  exp_d5 <- c(5.05133760899019, 5.21641483083995, 5.4131497041096, 5.62188922194764,
              5.75496241915971, 5.98260715969298, 6.18675520061515, 7.41160508326157,
              9.36626494375473, 10.5935673083787, 10.8441290040006)
  fd_exp1_mpseNP <- delay_model(exp_d5, method = "MPSE", profiled = FALSE)
  fd_exp1_mpseP  <- delay_model(exp_d5, method = "MPSE", profiled = TRUE)
  expect_identical(fd_exp1_mpseNP$optimizer$convergence, expected = 0L)
  expect_identical(fd_exp1_mpseP$optimizer$convergence, expected = 0L)
  expect_named(fd_exp1_mpseNP$optimizer$parOpt, expected = c("delay1_tr", "rate1_tr"))
  expect_named(fd_exp1_mpseP$optimizer$parOpt, expected = "delay1_tr")
  expect_equal(coef(fd_exp1_mpseP), expected = coef(fd_exp1_mpseNP), tolerance = .01)


  # MLE fits -----------------------------------------------------------

  # MLEn = naive MLE
  fd_exp_MLEn_NP <- delay_model(exp_d9, distribution = 'expon', method = 'MLEn')
  fd_exp_MLEn_P <- delay_model(exp_d9, distribution = 'expon', method = 'MLEn', profiled = TRUE)

  expect_type(fd_exp_MLEn_NP$data, type = 'double')
  expect_identical(length(fd_exp_MLEn_NP$data), expected = length(exp_d9))
  expect_identical(fd_exp_MLEn_NP$method, expected = "MLEn")
  expect_identical(length(coef(fd_exp_MLEn_NP)), expected = 2L)
  expect_named(coef(fd_exp_MLEn_NP), expected = c("delay1", "rate1"))
  expect_named(fd_exp_MLEn_NP$optimizer$parOpt, expected = c("delay1_tr", "rate1_tr"))

  # MLEn uses analytical solution when in single group mode
  expect_named(fd_exp_MLEn_NP$optimizer, expected = c("parOpt", "valOpt", 'profiled', "methodOpt", "convergence", "message", "counts")) ##analytic with no "optim_args"))
  expect_identical(purrr::chuck(fd_exp_MLEn_NP, 'optimizer', 'methodOpt'), expected = "analytic")
  expect_identical(purrr::chuck(fd_exp_MLEn_NP, 'optimizer', 'convergence'), expected = 0L)
  local({
    fd_exp_MLEn_NPopt <- attr(fd_exp_MLEn_NP$objFun, which = 'opt', exact = TRUE)
    expect_type(fd_exp_MLEn_NPopt, type = 'list')
    expect_named(fd_exp_MLEn_NPopt, expected = c("par_orig", 'par', 'value', "methodOpt", 'convergence', 'message', 'counts'))
  })
  # check objective function
  # objective function is for transformed parameters:
  purrr::walk2(.x = runif(n=7, min=0.1, max=9),   #delay1
               .y = runif(n=7, min=0.001, max=2), #rate1
               .f = ~ expect_equal(- length(exp_d9) * (log(.y) - .y * (mean(exp_d9) - .x)),
                                   expected = fd_exp_MLEn_NP$objFun(pars = c(delay1=.x, rate1=.y), criterion = TRUE))
  )

  expect_type(fd_exp_MLEn_P$data, type = 'double')
  expect_identical(length(fd_exp_MLEn_P$data), expected = length(exp_d9))
  expect_identical(fd_exp_MLEn_P$method, expected = "MLEn")
  expect_true(fd_exp_MLEn_P$optimizer$profiled)
  expect_identical(length(coef(fd_exp_MLEn_P)), expected = 2L)
  expect_named(coef(fd_exp_MLEn_P), expected = c("delay1", "rate1"))
  expect_named(fd_exp_MLEn_P$optimizer$parOpt, expected = "delay1_tr")



  # delay1 estimate = min(x)
  expect_identical(coef(fd_exp_MLEn_NP)[['delay1']], expected = exp_d9[[1L]])
  # MLEn's later delay is compensated for by higher estimated rate
  expect_true(all(coef(fd_exp_MLEn_NP) >  coef_exp))
  expect_equal(coef(fd_exp_MLEn_NP)[['rate1']], expected = (mean(exp_d9) - min(exp_d9))**-1)


  # MLEc
  fd_exp_MLEc_NP <- delay_model(exp_d9, distribution = 'expon', method = 'MLEc')
  fd_exp_MLEc_P <- delay_model(exp_d9, distribution = 'expon', method = 'MLEc', profiled = TRUE)

  expect_type(fd_exp_MLEc_NP$data, type = 'double')
  expect_identical(length(fd_exp_MLEc_NP$data), expected = length(exp_d9))
  expect_identical(length(coef(fd_exp_MLEc_NP)), expected = 2L)
  expect_named(coef(fd_exp_MLEc_NP), expected = c("delay1", "rate1"))
  expect_named(fd_exp_MLEc_NP$optimizer$parOpt, expected = c("delay1_tr", "rate1_tr"))

  expect_named(fd_exp_MLEc_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", "convergence", "message", "counts", "optim_args"))
  expect_identical(purrr::chuck(fd_exp_MLEc_NP, 'optimizer', 'convergence'), expected = 0L)
  # no exact solution for MLEc
  expect_null(attr(fd_exp_MLEc_NP$objFun, which = 'opt', exact = TRUE))

  expect_type(fd_exp_MLEc_P$data, type = 'double')
  expect_identical(length(fd_exp_MLEc_P$data), expected = length(exp_d9))
  expect_identical(length(coef(fd_exp_MLEc_P)), expected = 2L)
  expect_named(coef(fd_exp_MLEc_P), expected = c("delay1", "rate1"))
  expect_named(fd_exp_MLEc_P$optimizer$parOpt, expected = "delay1_tr")
  # profiled variant is an easier optimization
  expect_true(all(fd_exp_MLEc_P$optimizer$counts < fd_exp_MLEc_NP$optimizer$counts))
  # quite similar coefficients
  expect_equal(coef(fd_exp_MLEc_P), expected = coef(fd_exp_MLEc_NP), tolerance = .01)

  # MLEw
  fd_exp_MLEw <- delay_model(exp_d9, distribution = "expon", method = "MLEw")
  expect_identical(fd_exp_MLEw$method, expected = "MLEw")
  expect_true(fd_exp_MLEw$optimizer$profiled)
  expect_type(fd_exp_MLEw$data, type = "double")
  expect_identical(length(fd_exp_MLEw$data), expected = length(exp_d9))
  expect_identical(fd_exp_MLEw$optimizer$convergence, expected = 0L)
  expect_named(coef(fd_exp_MLEw), expected = c("delay1", "rate1"))
  expect_named(fd_exp_MLEw$optimizer$parOpt, expected = "delay1_tr")
  expect_equal(coef(fd_exp_MLEw), coef(fd_exp_MLEc_P), tolerance = .1) # similar coefficients


  # 2nd group ---------------------------------------------------------------

  set.seed(20220429)
  exp_d10 <- rexp_delayed(27L, delay1 = 10, rate1 = .9)

  fd_exp2 <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon")
  coef_exp2 <- coef(fd_exp2)

  expect_type(fd_exp2$data, type = 'list')
  # data gets sorted
  expect_identical(fd_exp2$data$x, sort(exp_d9))
  expect_identical(fd_exp2$data$y, sort(exp_d10))
  expect_named(fd_exp2$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2, 'optimizer', 'convergence'), expected = 0L)
  expect_type(purrr::chuck(fd_exp2, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef_exp2[1:2]), expected = as.numeric(coef_exp), tolerance = .01)

  # no delay in 2nd group
  set.seed(20221124)
  exp_d0 <- rexp_delayed(13L, delay1=0, rate1=1.1)
  fd_exp2nd <- delay_model(x = exp_d9, y = exp_d0, distribution = "expon") #nd = no delay
  coef_exp2nd <- coef(fd_exp2nd)

  expect_identical(purrr::chuck(fd_exp2nd, "optimizer", "convergence"), expected = 0L)
  expect_equal(coef_exp2nd[[3]], 0) # no delay in 3rd group

  # MLE fits -----
  fd_exp2_MLEn_NP <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", method = "MLEn")
  coef_exp2_MLEn_NP <- coef(fd_exp2_MLEn_NP)
  expect_type(fd_exp2_MLEn_NP$data, type = "list")
  expect_identical(fd_exp2_MLEn_NP$data$y, sort(exp_d10))
  expect_named(fd_exp2_MLEn_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2_MLEn_NP, 'optimizer', 'convergence'), expected = 0L) # converged
  expect_type(purrr::chuck(fd_exp2_MLEn_NP, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef_exp2_MLEn_NP[1:2]), expected = as.numeric(coef(fd_exp_MLEn_NP)), tolerance = .01)
  # MLEn fits have larger (=later) delay and (to recover a higher rate)
  purrr::walk(.x = seq_along(coef(fd_exp2)), .f = ~ expect_gte(coef(fd_exp2_MLEn_NP)[.x], coef(fd_exp2)[.x]))

  fd_exp2_MLEn_P <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", method = "MLEn", profiled = TRUE)
  coef_exp2_MLEn_P <- coef(fd_exp2_MLEn_P)
  expect_type(fd_exp2_MLEn_P$data, type = "list")
  expect_identical(fd_exp2_MLEn_P$data$y, sort(exp_d10))
  expect_named(fd_exp2_MLEn_P$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2_MLEn_P, 'optimizer', 'convergence'), expected = 0L) # converged
  expect_type(purrr::chuck(fd_exp2_MLEn_P, 'optimizer', 'optim_args'), type = 'list')
  expect_equal(coef_exp2_MLEn_P, expected = coef_exp2_MLEn_NP, tolerance = .01) # parameters are quite similar

  fd_exp2_MLEc_NP <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", method = "MLEc")
  expect_type(fd_exp2_MLEc_NP$data, type = "list")
  expect_identical(fd_exp2_MLEc_NP$data$y, sort(exp_d10))
  expect_named(fd_exp2_MLEc_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2_MLEc_NP, 'optimizer', 'convergence'), expected = 0L) # converged
  expect_type(purrr::chuck(fd_exp2_MLEc_NP, 'optimizer', 'optim_args'), type = 'list')
  # coefficient do not change much when adding a 2nd independent group and no binding
  expect_equal(as.numeric(coef(fd_exp2_MLEc_NP)[1:2]), expected = as.numeric(coef(fd_exp_MLEc_NP)), tolerance = .01)

  fd_exp2_MLEc_P <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", method = "MLEc", profiled = TRUE)
  expect_type(fd_exp2_MLEc_P$data, type = "list")
  expect_identical(fd_exp2_MLEc_P$data$y, sort(exp_d10))
  expect_named(fd_exp2_MLEc_P$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2_MLEc_P, 'optimizer', 'convergence'), expected = 0L) # converged
  expect_type(purrr::chuck(fd_exp2_MLEc_P, 'optimizer', 'optim_args'), type = 'list')
  expect_equal(coef(fd_exp2_MLEc_P), expected = coef(fd_exp2_MLEc_NP), tolerance = .01)


  # bind delay
  fd_exp2b_NP <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = "delay1")
  coef_exp2b_NP <- coef(fd_exp2b_NP)

  expect_named(coef_exp2b_NP, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_NP, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b_NP[[1L]],
               expected = min(coef_exp2[grepl(pattern = "delay1", names(coef_exp2), fixed = TRUE)]),
               tolerance = .01)

  fd_exp2b_P <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = "delay1", profiled = TRUE)
  coef_exp2b_P <- coef(fd_exp2b_P)
  expect_named(coef_exp2b_P, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_identical(purrr::chuck(fd_exp2b_P, 'optimizer', 'convergence'), expected = 0L)
  expect_equal(coef_exp2b_P, expected = coef_exp2b_NP, tolerance = .05) # profiled=T/F: parameters are similar

  # bind delay with MLEn
  fd_exp2b_MLEn_NP <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = "delay1", method = "MLEn")
  coef_exp2b_MLEn_NP <- coef(fd_exp2b_MLEn_NP)

  expect_named(coef_exp2b_MLEn_NP, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_MLEn_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_MLEn_NP, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b_MLEn_NP[[1L]],
               expected = min(coef(fd_exp2_MLEn_NP)[grepl(pattern = "delay1", names(coef(fd_exp2_MLEn_NP)), fixed = TRUE)]),
               tolerance = .01)

  fd_exp2b_MLEn_P <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = "delay1", method = "MLEn", profiled = TRUE)
  coef_exp2b_MLEn_P <- coef(fd_exp2b_MLEn_P)
  expect_named(coef_exp2b_MLEn_P, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_MLEn_P$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_MLEn_P, 'optimizer', 'convergence'), expected = 0L)
  expect_equal(coef_exp2b_MLEn_P, expected = coef_exp2b_MLEn_NP, tolerance = 1e-5) # profiled=T/F: similar coefficients

  # bind delay with MLEc
  fd_exp2b_MLEc_NP <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = "delay1", method = "MLEc")
  coef_exp2b_MLEc_NP <- coef(fd_exp2b_MLEc_NP)

  expect_named(coef_exp2b_MLEc_NP, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_MLEc_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_MLEc_NP, 'optimizer', 'convergence'), expected = 0L)
  # the bound delay is near the minimum of the two delay estimates from the individual group fits
  expect_equal(coef_exp2b_MLEc_NP[[1L]],
               expected = min(coef(fd_exp2_MLEc_NP)[grepl(pattern = "delay1", names(coef(fd_exp2_MLEc_NP)), fixed = TRUE)]),
               tolerance = .01)

  fd_exp2b_MLEc_P <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = "delay1", method = "MLEc", profiled = TRUE)
  coef_exp2b_MLEc_P <- coef(fd_exp2b_MLEc_P)

  expect_named(coef_exp2b_MLEc_P, expected = c("delay1", "rate1.x", "rate1.y"))
  expect_named(fd_exp2b_MLEc_P$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_identical(purrr::chuck(fd_exp2b_MLEc_P, 'optimizer', 'convergence'), expected = 0L)
  expect_equal(coef_exp2b_MLEc_P, expected = coef_exp2b_MLEc_NP, tolerance = .005) # profiled=T/F: similar coefficients


  # bind delay + rate
  fd_exp2c_NP <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = c("delay1", "rate1"))
  coef_exp2c_NP <- coef(fd_exp2c_NP)
  expect_named(coef_exp2c_NP, expected = c("delay1", "rate1"))

  # bind: order of parameters does not matter
  fd_exp2c_NP2 <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = c("rate1", "delay1"))
  expect_identical(coef(fd_exp2c_NP2), expected = coef(fd_exp2c_NP))
  expect_identical(fd_exp2c_NP2$optimizer$convergence, expected = 0L)
  expect_identical(fd_exp2c_NP2$bind, fd_exp2c_NP$bind)

  # profiling and full bind throws a warning:
  #+if rate1 is to be profiled out but later inferred from the delay1-estimate and the observations per group it will not be the same
  fd_exp2c_P <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = c("delay1", "rate1"), profiled = TRUE)
  coef_exp2c_P <- coef(fd_exp2c_P)
  expect_identical(fd_exp2c_P$optimizer$convergence, expected = 0L)
  expect_equal(coef_exp2c_P, coef_exp2c_NP, tolerance = .03) # parameters are similar (but because of bind= & profiled rate1 is simply averaged)

  # data combined
  fd_exp2comb_NP <- delay_model(x = c(exp_d9, exp_d10), distribution = "expon")
  coef_exp2comb_NP <- coef(fd_exp2comb_NP)

  fd_exp2comb_P <- delay_model(x = c(exp_d9, exp_d10), distribution = "expon", profiled = TRUE)
  coef_exp2comb_P <- coef(fd_exp2comb_P)

  expect_named(coef_exp2comb_NP, expected = c("delay1", "rate1"))
  expect_identical(length(coef_exp2comb_NP), length(coef_exp2c_NP))
  # similar coefficients (but the criterion is not identical, weighted mean for both groups vs overall mean)
  expect_equal(coef_exp2c_NP[[1L]], coef_exp2comb_NP[[1L]], tolerance = .01)
  expect_equal(coef_exp2c_NP[[2L]], coef_exp2comb_NP[[2L]], tolerance = .04)

  expect_named(coef_exp2comb_P, expected = c("delay1", "rate1"))
  expect_equal(coef_exp2comb_P, coef_exp2comb_NP, tolerance = .01)

  # bind delay + rate for MLE
  fd_exp2c_MLEn <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = c("delay1", "rate1"), method = "MLEn")
  coef_exp2c_MLEn <- coef(fd_exp2c_MLEn)
  fd_exp2comb_MLEn <- delay_model(x = c(exp_d9, exp_d10), distribution = "expon", method = "MLEn")
  coef_exp2comb_MLEn <- coef(fd_exp2comb_MLEn)

  expect_named(coef_exp2c_MLEn, expected = c("delay1","rate1"))
  expect_identical(length(coef_exp2comb_MLEn), length(coef_exp2c_MLEn))
  purrr::walk(seq_along(coef_exp2c_MLEn), ~expect_equal(coef_exp2c_MLEn[[.x]], coef_exp2comb_MLEn[[.x]], tolerance = .05))

  # MLEc
  fd_exp2c_MLEc <- delay_model(x = exp_d9, y = exp_d10, distribution = "expon", bind = c("delay1", "rate1"), method = "MLEc")
  coef_exp2c_MLEc <- coef(fd_exp2c_MLEc)
  fd_exp2comb_MLEc <- delay_model(x = c(exp_d9, exp_d10), distribution = "expon", method = "MLEc")
  coef_exp2comb_MLEc <- coef(fd_exp2comb_MLEc)

  expect_named(coef_exp2c_MLEc, expected = c("delay1","rate1"))
  expect_identical(length(coef_exp2comb_MLEc), length(coef_exp2c_MLEc))
  purrr::walk(seq_along(coef_exp2c_MLEc), ~expect_equal(coef_exp2c_MLEc[[.x]], coef_exp2comb_MLEc[[.x]], tolerance = .05))

})



test_that("Fit delayed Weibull", {

  # single group ----

  # susquehanna is an example dataset within incubate
  fd_maxFl <- delay_model(susquehanna, distribution = "weib")
  coef_maxFl <- coef(fd_maxFl)

  expect_identical(purrr::chuck(fd_maxFl, 'optimizer', 'convergence'), expected = 0L)
  expect_lte(fd_maxFl$criterion, expected = 3.0987)
  expect_equal(coef_maxFl, expected = c(delay1=0.244, shape1=1.310, scale1=.202), tolerance = .005)

  # MLE-based fits to susquehanna --
  fd_maxFl_MLEn_NP <- delay_model(susquehanna, distribution = "weib", method = "MLEn")
  coef_maxFl_MLEn <- coef(fd_maxFl_MLEn_NP)
  expect_false(fd_maxFl_MLEn_NP$optimizer$profiled)
  expect_identical(purrr::chuck(fd_maxFl_MLEn_NP, "optimizer", "convergence"), expected = 0L)
  expect_named(fd_maxFl_MLEn_NP$optimizer, expected = c("parOpt", "valOpt", "profiled", "methodOpt", 'convergence', 'message', 'counts', 'optim_args'))
  expect_gte(coef_maxFl_MLEn[["delay1"]], coef_maxFl[["delay1"]])
  expect_lte(coef_maxFl_MLEn[["scale1"]], coef_maxFl[["scale1"]])

  # check the objective function
  purrr::pwalk(.l = list(delay1=runif(7, max=0.25),
                         shape1=runif(7, max=3),
                         scale1=runif(7, max=5)),
               .f = ~ {
                 xc <- susquehanna - ..1
                 expect_equal(fd_maxFl_MLEn_NP$objFun(c(delay1=..1, shape1=..2, scale1=..3), criterion = TRUE),
                              expected = -length(susquehanna) * (log(..2) - ..2 * log(..3) + (..2 - 1L) * mean(log(xc)) - mean(xc**..2)/..3**..2))
               })

  fd_maxFl_MLEn_P <- delay_model(susquehanna, distribution = "weibu", method = "MLEn", profiled = TRUE)
  expect_true(fd_maxFl_MLEn_P$optimizer$profiled)
  expect_identical(fd_maxFl_MLEn_P$optimizer$convergence, expected = 0L)
  expect_named(coef(fd_maxFl_MLEn_P), expected = c("delay1", "shape1", "scale1"))
  expect_equal(coef(fd_maxFl_MLEn_P), expected = coef_maxFl_MLEn, tolerance = .001)
  expect_equal(fd_maxFl_MLEn_P$objFun(pars = fd_maxFl_MLEn_P$par, criterion = TRUE), fd_maxFl_MLEn_NP$criterion, tolerance = .001)

  # MLEn profiling by hand, using the indirect criterion of min(f')
  llProfObjFun_ind <- function(theta){
    stopifnot( is.numeric(theta), length(theta) == 2L )
    a <- theta[[1L]]
    k <- theta[[2L]]
    (1/k + mean(log(susquehanna-a)) - sum(log(susquehanna - a) * (susquehanna - a)**k)/sum((susquehanna-a)**k))^2 +
      # first factor inverse of harmonic mean of (susquehanna - a)
      (mean(1/(susquehanna-a)) * sum((susquehanna-a)**k)/sum((susquehanna-a)**(k-1)) - k/(k-1))^2
  }

  # the indirect objective function variant (using min(f')) is indeed close to 0 at the fit for the coefficients
  expect_equal(llProfObjFun_ind(coef(fd_maxFl_MLEn_P)[1:2]), expected = 0, tolerance = .01)

  # manual optimization, using the indirect objective function
  opt_maxFl_MLEn_Pman <- optim(par = c(a=0.255, k=1.8), #c(a=0.165, k=exp(1.3847)), #c(a=0.25, k=1.5),
                               fn = llProfObjFun_ind, method = "L-BFGS-B",
                      lower = c(0, 1 + 1.49e-8), upper = c(min(susquehanna)-1.49e-8, +Inf))
  coef_maxFl_MLEnp_man <- purrr::set_names(c(opt_maxFl_MLEn_Pman$par, mean((susquehanna-opt_maxFl_MLEn_Pman$par["a"])^opt_maxFl_MLEn_Pman$par["k"])^(1/opt_maxFl_MLEn_Pman$par["k"])),
                                  nm = c("delay1", "shape1", "scale1"))
  # estimates are only **roughly** equal
  expect_equal(coef_maxFl_MLEnp_man, coef(fd_maxFl_MLEn_P), tolerance = .25)

  fd_maxFl_MLEw <- delay_model(x = susquehanna, distribution = "weib", method = "MLEw")
  expect_identical(fd_maxFl_MLEw$optimizer$convergence, expected = 0L)
  expect_named(coef(fd_maxFl_MLEw), expected = c("delay1", "shape1", "scale1"))
  expect_gte(fd_maxFl_MLEw$criterion, fd_maxFl_MLEn_NP$criterion) #MLEn directly optimizes the criterion

  #llProfObjFun(coef_maxFl_MLEnp_man[1:2])
  #llProfObjFun(coef(fd_maxFl_MLEnp)[1:2])

  # # visualize the optimization function landscape
  # objFunLS_maxFl_MLEnp <- tidyr::expand_grid(a = seq.int(0, .26, length.out = 17),
  #                                            k = seq.int(1.3, 3, length.out = 17)) %>%
  #   #head %>%
  #   rowwise() %>%
  #   dplyr::mutate(ll = llProfObjFun(theta = c(a,k))) %>%
  #   ungroup()
  #
  # ggplot(objFunLS_maxFl_MLEnp, aes(x = a, y = k, z = log(ll+.1), colour = after_stat(level))) +
  #   geom_contour() +
  #   scale_colour_distiller(palette = "YlGn", direction = -1)



  # pollution is an example dataset within incubate
  fd_poll <- delay_model(pollution, distribution = "weib")
  objFun_poll <- fd_poll$objFun
  coef_poll <- coef(fd_poll)

  expect_identical(purrr::chuck(fd_poll, 'optimizer', 'convergence'), expected = 0L)
  expect_identical(length(coef_poll), expected = 3L)
  expect_lt(objFun_poll(coef_poll, criterion = TRUE), expected = 3.75)
  expect_equal(coef_poll, expected = c(delay1=1085, shape1=0.95, scale1=6562), tolerance = .001)

  # MLE-based fits for pollution
  fd_poll_MLEn <- delay_model(pollution, distribution = "weibu", method = "MLEn", profile = FALSE)
  coef_poll_MLEn <- coef(fd_poll_MLEn)
  expect_identical(fd_poll_MLEn$optimizer$convergence, expected = 0L)
  expect_named(coef_poll_MLEn, expected = c("delay1", "shape1", "scale1"))
  # naive MLE gives later delay estimates than MPSE, close to minimal observation
  expect_gt(coef_poll_MLEn[["delay1"]], coef_poll[["delay1"]])
  expect_equal(coef_poll_MLEn[["delay1"]], expected = min(pollution), tolerance = 1e-4)
  # MLEn allows for shape estimates k below 1
  expect_lt(coef_poll_MLEn[["shape1"]], expected = 1)

  fd_poll_MLEnp <- delay_model(pollution, distribution = "weibu", method = "MLEn", profile = TRUE)
  coef_poll_MLEnp <- coef(fd_poll_MLEnp)
  expect_identical(fd_poll_MLEnp$optimizer$convergence, expected = 0L)
  expect_named(coef_poll_MLEnp, expected = c("delay1", "shape1", "scale1"))
  # MLEn and MLEnp give very similar parameter estimates
  expect_equal(coef_poll_MLEn, expected = coef_poll_MLEnp, tolerance = 1e-4)
  # profiled MLEn allows for shape estimates k below 1
  expect_lt(coef_poll_MLEnp[["shape1"]], expected = 1)

  fd_poll_MLEw <- delay_model(pollution, distribution = "weibu", method = "MLEw", profile = TRUE)
  expect_identical(fd_poll_MLEw$optimizer$convergence, expected = 0L)
  expect_true(fd_poll_MLEw$optimizer$profiled)
  expect_named(coef(fd_poll_MLEw), expected = c("delay1", "shape1", "scale1"))


  # Cousineau's numerical example, taken from Weibull with delay = 300, shape k = 2 and scale = 100
  x <- c(310, 342, 353, 365, 383, 393, 403, 412, 451, 456)
  densFun_wb <- incubate:::getDist(distribution = "weib", type = "density")

  fd_wbc_mlen <- delay_model(x = x, distribution = "weib", method = "MLEn")
  expect_equal(coef(fd_wbc_mlen), expected = c(delay1 = 274.8, shape1 = 2.80, scale1 = 126.0), tolerance = .001)
  fd_wbc_mlenp <- delay_model(x = x, distribution = "weib", method = "MLEn", profile = TRUE)
  expect_equal(coef(fd_wbc_mlenp), expected = c(delay1=280.9, shape1=2.62, scale1=119.0), tolerance = .05)
  # our implementation finds a smaller objective value (which we minimize)
  expect_lte(fd_wbc_mlenp$optimizer$valOpt, fd_wbc_mlenp$objFun(c(delay1=280.9, shape1=log(2.62))))
  # but log-likelihood is in fact quite similar (actually, Cousineau's solution is slightly worse)
  expect_equal(fd_wbc_mlenp$criterion, fd_wbc_mlenp$objFun(pars = c(delay1 = 280.9, shape1=2.62, scale1=119.0), criterion = TRUE), tolerance = .01)
  expect_equal(fd_wbc_mlenp$criterion, -sum(densFun_wb(x, delay1 = 280.9, shape1=2.62, scale1=119.0, log = TRUE)), tolerance = .01)
  fd_wbc_mlewp <- delay_model(x = x, distribution = "weib", method = "MLEw", profile = TRUE)
  expect_true(fd_wbc_mlewp$optimizer$profiled)
  expect_equal(coef(fd_wbc_mlewp), c(delay1=283.7, shape1=2.29, scale1=116.0), tolerance = .04)
  # our implementation finds a smaller objective value.
  expect_lte(fd_wbc_mlewp$optimizer$valOpt, fd_wbc_mlewp$objFun(c(delay1=283.7, shape1=log(2.29))))
  # but log-likelihood is in fact quite similar (actually, Cousineau's solution is slightly better)
  expect_equal(fd_wbc_mlewp$criterion, -sum(densFun_wb(x, delay1=283.7, shape1=2.29, scale1=116.0, log = TRUE)), tolerance = .01)


  # two groups -----

  fd_wb2 <- incubate::delay_model(x = susquehanna, y = pollution, distribution = "weib")
  coef_wb2 <- coef(fd_wb2)

  expect_identical(purrr::chuck(fd_wb2, 'optimizer', 'convergence'), expected = 0L)
  expect_identical(length(coef_wb2), expected = 2L*3L)
  expect_named(coef_wb2, expected = c("delay1.x", "shape1.x", "scale1.x", "delay1.y", "shape1.y", "scale1.y"))
  expect_equal(as.numeric(coef_wb2[1:3]), expected = as.numeric(coef_maxFl), tolerance = .001)
  # there might occur quite some differences in the coefficients
  expect_equal(as.numeric(coef_wb2[4:6]), expected = as.numeric(coef_poll), tolerance = .25)

  # MLE based fits
  fd_wb2_MLEn_NP <- incubate::delay_model(x = susquehanna, y = pollution, distribution = "weib", method = "MLEn")
  expect_identical(names(coef(fd_wb2_MLEn_NP)), expected = names(coef(fd_wb2)))
  expect_equal(coef(fd_wb2_MLEn_NP), expected = coef(fd_wb2), tolerance = .3)
  # MLEn has later delay estimates
  expect_gt(coef(fd_wb2_MLEn_NP)[["delay1.x"]], coef(fd_wb2)[["delay1.x"]])
  expect_gt(coef(fd_wb2_MLEn_NP)[["delay1.y"]], coef(fd_wb2)[["delay1.y"]])

  fd_wb2_MLEn_P <- delay_model(x = susquehanna, y = pollution, distribution = "weib", method = "MLEn", profile = TRUE)
  expect_identical(fd_wb2_MLEn_P$optimizer$convergence, expected = 0L)
  # roughly equal coefficients as compared to single fits
  expect_equal(coef(fd_wb2_MLEn_P, group = "x"), expected = coef(fd_maxFl_MLEn_P), tolerance = .2)
  expect_equal(coef(fd_wb2_MLEn_P, group = "y"), expected = coef(fd_poll_MLEnp), tolerance = .01)

  # MLEc
  fd_wb2_MLEc_NP <- delay_model(x = susquehanna, y = pollution, distribution = "weib", method = "MLEc")
  expect_named(coef(fd_wb2_MLEc_NP), expected = names(coef(fd_wb2)))
  expect_false(fd_wb2_MLEc_NP$optimizer$profiled)
  expect_equal(coef(fd_wb2_MLEc_NP), expected = coef(fd_wb2), tolerance = .3)
  expect_gt(coef(fd_wb2_MLEc_NP)[["delay1.x"]], coef(fd_wb2)[["delay1.x"]])
  expect_gt(coef(fd_wb2_MLEc_NP)[["delay1.y"]], coef(fd_wb2)[["delay1.y"]])

  fd_wb2_MLEc_P <- delay_model(x = susquehanna, y = pollution, distribution = "weib", method = "MLEc", profiled = TRUE)
  expect_true(fd_wb2_MLEc_P$optimizer$profiled)
  expect_identical(fd_wb2_MLEc_P$optimizer$convergence, expected = 0L)
  # MLEc profiling has very little impact on coefficient estimates
  expect_equal(coef(fd_wb2_MLEc_P), expected = coef(fd_wb2_MLEc_NP), tolerance = .01)

  # MLEw with two groups
  fd_wb2_MLEw_P <- incubate::delay_model(x = susquehanna, y = pollution, distribution = "weib", method = "MLEw", profile = TRUE)
  expect_identical(names(coef(fd_wb2_MLEw_P)), expected = names(coef(fd_wb2)))
  expect_equal(coef(fd_wb2_MLEw_P, group = "x"), expected = coef(fd_maxFl_MLEw), tolerance = .01)
  expect_equal(coef(fd_wb2_MLEw_P, group = "y"), expected = coef(fd_poll_MLEw), tolerance = .01)


  # two groups with binding
  set.seed(20210430)
  fd_wb2b <- delay_model(x = rweib_delayed(n=37, delay1 = 7, shape1 = 1.8, scale1 = 3),
                         y = rweib_delayed(n=51, delay1 = 5, shape1 = 1.2, scale1 = 1.5),
                         distribution = "weib", bind = "delay1")
  coef_wb2b <- coef(fd_wb2b)

  expect_identical(purrr::chuck(fd_wb2b, 'optimizer', 'convergence'), expected = 0L)
  #delay is bound
  expect_identical(length(coef_wb2b), expected = 2L*3L-1L)
  expect_named(coef_wb2b, c("delay1", "shape1.x", "scale1.x", "shape1.y", "scale1.y"))
  # expect a delay close to the minimum of the two true delay parameters
  expect_equal(coef_wb2b[[1L]], expected = 5, tolerance = .02)

  fd_wb2b_MLEn_NP <- delay_model(x = fd_wb2b$data$x, y = fd_wb2b$data$y, distribution = "weib", method = "MLEn", bind = "delay1")
  expect_identical(fd_wb2b_MLEn_NP$optimizer$convergence, expected = 0L)
  expect_named(coef(fd_wb2b_MLEn_NP), c("delay1", "shape1.x", "scale1.x", "shape1.y", "scale1.y"))
  expect_equal(coef(fd_wb2b_MLEn_NP), coef_wb2b, tolerance = .03)

  fd_wb2b_MLEn_P <- delay_model(x = fd_wb2b$data$x, y = fd_wb2b$data$y, distribution = "weib", method = "MLEn", profile = TRUE, bind = "delay1")
  expect_identical(fd_wb2b_MLEn_P$optimizer$convergence, expected = 0L)
  expect_equal(coef(fd_wb2b_MLEn_P), coef(fd_wb2b_MLEn_NP), tolerance = .01) # profiling has little effect on coefficients

  fd_wb2b_MLEc_NP <- delay_model(x = fd_wb2b$data$x, y = fd_wb2b$data$y, distribution = "weib", method = "MLEc", profile = FALSE, bind = "delay1")
  expect_identical(fd_wb2b_MLEc_NP$optimizer$convergence, expected = 0L)
  fd_wb2b_MLEc_P <- delay_model(x = fd_wb2b$data$x, y = fd_wb2b$data$y, distribution = "weib", method = "MLEc", profile = TRUE, bind = "delay1")
  expect_identical(fd_wb2b_MLEc_P$optimizer$convergence, expected = 0L)
  expect_named(coef(fd_wb2b_MLEc_P), c("delay1", "shape1.x", "scale1.x", "shape1.y", "scale1.y"))
  expect_equal(coef(fd_wb2b_MLEc_P), coef(fd_wb2b_MLEc_NP), tolerance = .001) # profiling has hardly an effect on the coefficients

  fd_wb2b_MLEw_P <- delay_model(x = fd_wb2b$data$x, y = fd_wb2b$data$y, distribution = "weib", method = "MLEw", profile = TRUE, bind = "delay1")
  expect_identical(fd_wb2b_MLEw_P$optimizer$convergence, expected = 0L)

  # bind shape
  fd_wb2bb <- delay_model(x = rweib_delayed(n=37, delay1 = 7, shape1 = 1.8, scale1 = 3),
                          y = rweib_delayed(n=51, delay1 = 5, shape1 = 1.5, scale1 = 1.5),
                          distribution = "weib", bind = "shape1")
  coef_wb2bb <- coef(fd_wb2bb)

  expect_identical(fd_wb2bb$optimizer$convergence, expected = 0L)
  expect_identical(length(coef_wb2bb), expected = 2L*3L-1L)
  expect_identical(names(coef_wb2bb), c("shape1", "delay1.x", "scale1.x", "delay1.y", "scale1.y"))

  fd_wb2bb_MLEn_NP <- delay_model(x = fd_wb2bb$data$x, y = fd_wb2bb$data$y, distribution = "weib", method = "MLEn", bind = "shape1")
  expect_identical(fd_wb2bb_MLEn_NP$optimizer$convergence, expected = 0L)
  expect_identical(names(coef(fd_wb2bb_MLEn_NP)), c("shape1", "delay1.x", "scale1.x", "delay1.y", "scale1.y"))
  expect_equal(coef(fd_wb2bb_MLEn_NP), coef_wb2bb, tolerance = .05)

  fd_wb2bb_MLEn_P <- delay_model(x = fd_wb2bb$data$x, y = fd_wb2bb$data$y, distribution = "weib", method = "MLEn", profile = TRUE, bind = "shape1")
  expect_identical(fd_wb2bb_MLEn_P$optimizer$convergence, expected = 0L)
  expect_true(fd_wb2bb_MLEn_P$optimizer$profiled)
  expect_identical(names(coef(fd_wb2bb_MLEn_P)), c("shape1", "delay1.x", "scale1.x", "delay1.y", "scale1.y"))
  expect_equal(coef(fd_wb2bb_MLEn_P), coef(fd_wb2bb_MLEn_NP), tolerance = .01) # profiling has hardly an effect

  fd_wb2bb_MLEc_NP <- delay_model(x = fd_wb2bb$data$x, y = fd_wb2bb$data$y, distribution = "weib", method = "MLEc", profile = FALSE, bind = "shape1")
  expect_identical(fd_wb2bb_MLEc_NP$optimizer$convergence, expected = 0L)
  expect_identical(names(coef(fd_wb2bb_MLEc_NP)), c("shape1", "delay1.x", "scale1.x", "delay1.y", "scale1.y"))
  fd_wb2bb_MLEc_P <- delay_model(x = fd_wb2bb$data$x, y = fd_wb2bb$data$y, distribution = "weib", method = "MLEc", profile = TRUE, bind = "shape1")
  expect_identical(fd_wb2bb_MLEc_P$optimizer$convergence, expected = 0L)
  expect_equal(coef(fd_wb2bb_MLEc_P), coef(fd_wb2bb_MLEc_NP), tolerance = .001) # profiling has hardly an effect on the coefficients

  fd_wb2bb_MLEw_P <- delay_model(x = fd_wb2bb$data$x, y = fd_wb2bb$data$y, distribution = "weib", method = "MLEw", profile = TRUE, bind = "shape1")
  expect_identical(fd_wb2bb_MLEw_P$optimizer$convergence, expected = 0L)
  expect_identical(names(coef(fd_wb2bb_MLEw_P)), c("shape1", "delay1.x", "scale1.x", "delay1.y", "scale1.y"))
  expect_equal(coef(fd_wb2bb_MLEw_P), coef(fd_wb2bb_MLEc_P), tolerance = .1) # MLEw and MLEc coefficients are only rougly similar
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
