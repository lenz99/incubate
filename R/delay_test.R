
#' Goodness-of-fit (GOF) test statistic.
#'
#' The GOF-test is performed for a fitted delay-model.
#' There are different GOF-tests implemented:
#' * __Moran GOF__ is based on spacings, like the MPSE-criterion itself.
#' * __Pearson GOF__ uses categories and compares observed to expected frequencies.
#'
#' @param delayFit delay_model fit
#' @param method character(1). which method to use for GOF. Default is 'moran'.
#' @return An `htest`-object containing the GOF-test result
#' @export
test_GOF <- function(delayFit, method = c('moran', 'pearson')){

  stopifnot( inherits(delayFit, what = 'incubate_fit') )
  method <- match.arg(method)

  if (delayFit$method != 'MPSE'){
    stop('Goodness-of-fit test only supported for models that are fit with maximum product of spacings estimation (MPSE)!')
  }

  twoGroup <- isTRUE(delayFit$twoGroup)
  data_name <- if (twoGroup) paste(names(delayFit$data), collapse = ' and ') else 'x'
  nObs <- if (twoGroup) lengths(delayFit$data) else length(delayFit$data)
  params_tr <- coef.incubate_fit(delayFit, transformed = TRUE)
  k <- length(params_tr)


  # required variables
  meth <- statist <- p_val <- NULL

  switch(method,
         moran = {
           # Moran's GOF test
           meth <- "Moran's Goodness-of-fit (GOF) test"

           EUL_MAS <- -digamma(1L)

           # Moran test statistic is negative sum of logged spacings, see Cheng & Stephens (1989)
           # @param mseCrit: the negative avg logged cumulative spacings, length 1 or 2
           # @param n nbr of observations, length 1 or 2
           # @param k nbr of parameters to be estimated
           # @return Moran's test statistic, length 1 or 2
           testStat_mo <- function(mseCrit, n, k){
             mo_m <- (n+1L) * (log(n+1L) + EUL_MAS) - .5 - (12L*(n + 1L))**-1L
             mo_v <- (n+1L) * (pi**2L / 6L - 1L) - .5 - (6L*(n + 1L))**-1L

             C1 <- mo_m - sqrt(.5 * n * mo_v)
             C2 <- sqrt(mo_v / (2L*n))

             # factor (n+1) to go from -avg to -sum
             (mseCrit * (n+1L) + .5 * k - C1) / C2
           }# fun


           statist <- if (twoGroup){ ##  && length(delayFit$bind) < length(oNames) # not needed!?
             # sum of two independent chi-sq. is chi-sq
             c(`X^2` = sum(testStat_mo(mseCrit = delayFit$objFun(pars = params_tr, aggregated = FALSE), #criterion per group
                                       n = nObs, k = k/2)) )
           } else {
             # single group
             c(`X^2` = testStat_mo(mseCrit = delayFit$val, n = nObs, k = k) )
           }

           p_val <- stats::pchisq(q = statist, df = sum(nObs), lower.tail = FALSE)
         },

         pearson = {
           # Pearson GOF-test
           meth <- "Pearson's Goodness-of-fit (GOF) test" # (per group) ## this is our standard

           # under H0, expect frequency counts according to uniform distribution
           #+nbr of classes as recommended by David S. Moore (Tests of Chi-squared Type, 1986)

           testStat_pe <- function(datr, nCl) {
             tab_transf <- tabulate(findInterval(datr,
                                   vec = seq.int(from=0L, to=1L, length.out = nCl+1L),
                                   rightmost.closed = TRUE, all.inside = TRUE), nbins = nCl)
             sum((tab_transf - mean(tab_transf))**2L) / mean(tab_transf)
           }


           nCl <- NA_real_ #dummy value to start with

           statist <- if (twoGroup) {
             nCl <- pmax.int(k/2 + 2L, ceiling(2L * nObs**.4))
             # sum of two chi-square test statistics
             c(`X^2` = sum(purrr::map2_dbl(.x = transform(delayFit), .y = nCl, .f = testStat_pe)))
           } else {
             nCl <- max(k + 2L, ceiling(2L * nObs**.4))
             c(`X^2` = testStat_pe(datr = transform(delayFit), nCl = nCl))
             # inference based on Chi-square distribution.
             # use adjusted degrees of freedom (loose one df for each parameter estimated)
           }
           p_val <- stats::pchisq(q = statist, df = sum(nCl) - k - (1L+twoGroup), lower.tail = FALSE)
         },

         AD =, ad =, anderson = {
           # EDF-based GOF-test
           # Anderson-Darling (AD) test statistic
           # cf Stephens, Tests based on EDF Statistics p.101, (4.2)

           meth <- 'Anderson-Darling Goodness-of-fit (GOF) test (per group)'

           testStat_ad <- function(datr, n){
             i <- seq_along(datr)
             # in fact, A2 utilizes rev-order in its 2nd summation term
             -n - mean((2L * i - 1L) * log(datr) + (2L*n - (2L*i-1L))*log(1-datr))
           }

           A2 <- if (twoGroup) {
             purrr::map2_dbl(.x = transform(delayFit), .y = nObs, .f = testStat_ad)
           } else {
             testStat_ad(datr = transform(delayFit), n = nObs)
           }

           p_val <- if ( delayFit$distribution == 'exponential' ){
             # modification for Exponential (cf Stephens, Table 4.14, p.138)
             # the correction factor approaches 1 from above.
             # We keep the number of N as all observations, independent of the number of parameters estimated in the null-model.
             # QQQ Should we increase N by the number of parameters p estimated less 2 ( p -2 because 2 parameters are estimated in standard delayed exponential)
             A2_mod <- A2 * pmax.int(1L, 1L + 5.4 / nObs - 11 / nObs**2L)

             .ad_pval[['exponential']](A2_mod)

           } else {
             stopifnot( delayFit$distribution == 'weibull' )

             # P-value for Weibull based on Lockhart, 1994 (Table 1)
             # interpolation model on logits using critical value and inverse of shape parameter
             params_ntr <- coef.incubate_fit(delayFit, transformed = FALSE)
             .ad_pval[['weibull']](A2, params_ntr[grepl('shape', names(params_ntr), fixed = TRUE)])
           }

           if (twoGroup){
             A2 <- paste(signif(A2, 4), collapse = ' and ')
             # use Liptak to bring both P-values of AD-tests per group together
             p_val <- stats::pnorm(sum(sqrt(nObs) * stats::qnorm(p_val)) / sqrt(sum(nObs)))
           }

           statist <- c(`A^2` = A2)

         },
         # catch all
         stop('This GOF-test method is not supported!')
  )



  # return test object

  # stats:::print.htest recognizes:
  #+parameter, alternative, null.value, conf.int, estimate
  structure(
    list(method = meth, data.name = data_name,
         statistic = statist, p.value = as.numeric(p_val)),
    class = 'htest')
}


#' Test the difference for delay model parameter(s) between two uncorrelated groups, based on maximum product of spacings estimation (MPSE).
#'
#' It is in fact a model comparison between a null model where the parameters are enforced to be equal and an unconstrained full model.
#' As test statistic we use twice the difference in best (=lowest) objective function value, i.e. 2 * (`val_0` - `val_1`).
#' This is reminiscent of a likelihood ratio test statistic albeit the objective function is not a negative log-likelihood
#' but the negative of the maximum product spacing metric.
#'
#' High values of this difference speak against the null-model (i.e. high `val_0` indicates bad fit under 0-model and low values of `val_1` indicate a good fit under the more general model1.
#' The test is implemented as a parametric bootstrap test, i.e. we
#'
#' 1. take given null-model fit as ground truth
#' 2. regenerate data according to this model.
#' 3. recalculate the test statistic
#' 4. appraise the observed test statistic in light of the generated distribution under H0
#'
#'
#' @param x data from reference/control group.
#' @param y data from the treatment group.
#' @param distribution character(1). Name of the parametric delay distribution to use.
#' @param twoPhase logical(1). Do we model two phases per group? Default is `FALSE`, i.e. a single delay phase per group.
#' @param param character. Names of parameters to test difference for. Default value is `'delay1'`.
#' @param R numeric(1). Number of bootstrap samples to evaluate the distribution of the test statistic.
#' @param ties character. How to handle ties in data vector of a group?
#' @param type character. Which type of tests to perform?
#' @param chiSqApprox logical flag. In bootstrap, should we calculate the approximate degrees of freedom for the distribution of the test statistic under H0?
#' @param verbose numeric. How many details are requested? Higher value means more details. 0=off, no details.
#' @return list with the results of the test. Element P contains the different P-values, for instance from parametric bootstrap
#' @export
test_diff <- function(x, y=stop('Provide data for group y!'), distribution = c("exponential", "weibull"), twoPhase = FALSE, param = "delay1", R = 400,
                      ties = c('density', 'equidist', 'random', 'error'), type = c('all', 'bootstrap', 'gof', 'moran', 'pearson', 'lr', 'lr_pp'),
                      chiSqApprox = FALSE, verbose = 0) {
  STRICT <- TRUE #accept models only if they converged flawlessly, i.e., convergence=0
  TOL_CRIT <- 1e-7
  distribution <- match.arg(arg = distribution)
  ties <- match.arg(arg = ties)
  type <- tolower(type)
  type <- match.arg(arg = type)
  onames <- getDist(distribution = distribution, type = "param", twoPhase = FALSE, twoGroup = FALSE, transformed = FALSE)
  stopifnot( is.numeric(x), length(x) > length(onames), is.numeric(y), length(y) > length(onames) )
  stopifnot( is.numeric(R), length(R) == 1L, R >= 1L )
  stopifnot( is.character(param) )
  if (is.logical(verbose)) verbose <- as.numeric(verbose)
  if ( is.null(verbose) || ! is.numeric(verbose) || ! is.finite(verbose) ) verbose <- 0
  verbose <- verbose[[1L]]

  # parameters to test differences
  if (any(grepl(pattern = "_tr", param, fixed = TRUE))){
    stop("Parameter names in param= refer to the distribution parameters and not to the transformed parameters of the objective function.", call. = FALSE)
  }

  # translate convenience names (for single phase) to canonical names
  unNmbrdIdx <- !grepl(pattern = "[12]", param, fixed = FALSE)
  if (any(unNmbrdIdx)){
    param[unNmbrdIdx] <- paste0(param[unNmbrdIdx], "1") #interpret un-numbered parameters as referring to phase 1
    if (verbose > 0L) cat("The unnumbered parameter names in param= are taken to refer to initial phase and are translated to canonical parameter names.\n")
  }

  # retain only valid names in canonical order
  param <- intersect(onames, param)

  if (!length(param)){
    stop("Provide valid parameter names from the distribution to test for differences in two groups.", call. = FALSE)
  }

  # bitmask for test types
  testMask <- purrr::set_names(logical(5L), nm = c('bootstrap', 'pearson', 'moran', 'lr', 'lr_pp'))

  switch(EXPR = type,
         all = { testMask <- testMask | TRUE },
         # bootstrap + LR-tests (for stankovic results) #use better flags? like doBootstrap=, doGOF=, doLR=?!
         bootstrap = {testMask[c('bootstrap', 'lr', 'lr_pp')] <- TRUE},
         gof = {testMask[c('pearson', 'moran')] <- TRUE},
         moran = {testMask['moran'] <- TRUE},
         pearson = {testMask['pearson'] <- TRUE},
         lr = {testMask['lr'] <- TRUE},
         lr_pp = {testMask['lr_pp'] <- TRUE},
         stop('This type of tests is not supported!')
  )

  # Test statistic calculated from the given data and the model specification.
  #
  # The test statistic takes non-negative values.
  # High values of the test statistic speak in favour of H1:
  # @return list containing value of test statistic and null model fit. Or `NULL` in case of trouble.
  testStat <- function(x, y) {
    fit0 <- delay_model(x = x, y = y, distribution = distribution, twoPhase = twoPhase, method = 'MPSE', ties = ties, bind = param)
    fit1 <- delay_model(x = x, y = y, distribution = distribution, twoPhase = twoPhase, method = 'MPSE', ties = ties)

    if (is.null(fit0) || is.null(fit1)) return(invisible(NULL))

    # if the more restricted model (fit0) yields better fit (=lower criterion) than the more general model (fit1)
    #+we are in trouble, possibly due to non-convergence, e.g. convergence code 52
    #+we refit the general fit1 again using parameter-values from fit0
    if ( fit0[['val']] + TOL_CRIT < fit1[['val']] &&
         !is.null(fit1oa <- purrr::pluck(fit1, 'optimizer', 'optim_args')) ) {
      if (verbose > 0) warning('Restricted model with better fit than unrestricted model.', call. = FALSE)
      # re-run fit1 with start values based on fitted parameters of reduced model fit0
      stopifnot( is.list(fit1oa), 'par' %in% names(fit1oa) )

      coef0 <- coef.incubate_fit(fit0, transformed = TRUE)
      pn1 <- names(fit1[["optimizer"]][["parOpt"]])
      # take over optimization coefficients for start values of fit1
      # QQQ Would match() or pmatch() help avoid the for-loop?
      for (na0 in names(fit0[["optimizer"]][["parOpt"]])) fit1oa[["par"]][startsWith(pn1, prefix = na0)] <- coef0[[na0]]

      fit1oa[['control']][['parscale']] <- scalePars(parV = fit1oa[["par"]])
      fit1 <- update.incubate_fit(fit1, optim_args = fit1oa)

      if (is.null(fit1) || fit0[['val']] + TOL_CRIT < fit1[['val']]) {
        warning('Restricted model with better fit than unrestricted model even after refit of the unrestricted model!', call. = FALSE)
        return(invisible(NULL))
      }
    }# fi bad fit1

    # convergence of re-fits?
    # keep also fits with error-code 52: in my tests all those fits *looked* actually OK..
    # XXX try harder for non-convergence when testStat is called for the first time (to have basis model)
    if ( STRICT && (purrr::chuck(fit0, 'optimizer', 'convergence') != 0 || purrr::chuck(fit1, 'optimizer', 'convergence') != 0) ) return(invisible(NULL))

    # higher values of T speak in favour of H1:
    #   1. fit0 has high value (=bad fit)
    #   2. fit1 has low value (=good fit)
    list(val = 2L * max(0L, fit0[["val"]] - fit1[["val"]]),
         fit0 = fit0, fit1 = fit1)
  } #fn testStat

  # observed test statistic
  ts_obs <- testStat(x, y)
  if ( is.null(ts_obs) || ! is.list(ts_obs) || ! is.numeric(ts_obs[['val']]) || ts_obs[['val']] < -TOL_CRIT ){
    stop("Delay model failed for restricted null-model or free full model", call. = FALSE)
  }
  fit0 <- ts_obs[["fit0"]]
  fit1 <- ts_obs[["fit1"]]


  # P-values based GOF-tests
  #+ H0: simpler/restricted model 0 is sufficient
  #+ the GOF-test solely builds on fit0
  #+ take fitted parameters for both groups under null-model
  #+ and transform the observed data for both groups via cumulative distribution functions

  # spacings-based GOF-test
  GOF_mo0 <- if (testMask[['moran']]) test_GOF(delayFit = fit0, method = 'moran')
  GOF_mo1 <- if (testMask[['moran']]) test_GOF(delayFit = fit1, method = 'moran')
  #if (verbose > 0L) cat("Moran test stat: ", GOF_mo$statistic, "\n")

  # Pearson GOF-test based on Chi-square distribution.
  # under H0, expect counts according to uniform distribution
  GOF_pears0 <- if (testMask[['pearson']]) test_GOF(delayFit = fit0, method = 'pearson')
  GOF_pears1 <- if (testMask[['pearson']]) test_GOF(delayFit = fit1, method = 'pearson')



  t0_dist <- P_boot <- chisq_df_hat <- NULL

  if (testMask[['bootstrap']]){
    # parametric bootstrap:
    # generate R samples (x, y) by random sampling from the fitted H0-model (e.g. common delay),
    #+where all nuisance parameters are at there best fit
    # calculate the test statistic on the simulated data
    # estimate P as proportion of simulated test statistics that exceed the observed test statistic t_obs

    ranFun <- getDist(distribution, type = "r")
    # arguments to the random function generation
    ranFunArgsX <- c(list(n=length(x)), coef(fit0, group = "x", transformed = FALSE))
    ranFunArgsY <- c(list(n=length(y)), coef(fit0, group = "y", transformed = FALSE))

    retL <- 1L+(verbose>0L)
    t0_dist <- future.apply::future_vapply(X = seq.int(R), FUN.VALUE = double(retL),
                                           FUN = function(dummy){

                                             # generate new data according to given fitted null-model
                                             # sort is not needed here, as it goes through the whole pipeline (factory method)
                                             ts_boot <- testStat(x = rlang::exec(ranFun, !!! ranFunArgsX),
                                                                 y = rlang::exec(ranFun, !!! ranFunArgsY))
                                             if (is.null(ts_boot)) rep.int(NA_real_, times = retL) else
                                                 c(ts_boot[['val']],
                                                   purrr::chuck(ts_boot, 'fit0', 'optimizer', 'convergence'))[seq_len(retL)]

                                           }, future.seed = TRUE)

    if (verbose > 0L){
      stopifnot( NROW(t0_dist) == 2L )
      fit0_conv <- t0_dist[2L,]
      cat(glue(
        'Proportion of model failures: {as_percent(length(which(is.na(fit0_conv)))/length(fit0_conv))}',
        'Proportion of conv =  0: {as_percent(length(which(fit0_conv == 0))/ length(fit0_conv))}',
        'Proportion of conv = 52: {as_percent(length(which(fit0_conv == 52))/length(fit0_conv))}',
        .sep = '\n'), '\n')
      t0_dist <- t0_dist[1L, , drop=TRUE] #retain only ts_boot[['val']]
    }
    t0_dist <- t0_dist[is.finite(t0_dist)]

    if (chiSqApprox && length(t0_dist) > 101L){
      try(expr = {chisq_df_hat <- coef(MASS::fitdistr(x = t0_dist, densfun = "chi-squared",
                                                      start = list(df = length(param)),
                                                      method = "Brent", lower = .001, upper = 401))},
          silent = TRUE)
    }

    P_boot <- (1L + sum(t0_dist >= ts_obs[["val"]])) / (length(t0_dist)+1L)
  } # bootstrap

  ## P-value from Log-rank tests
  dat_2gr <- tibble::tibble(evtime = c(x,y),
                            group = rep.int(c("x", "y"), times = c(length(x), length(y))))
  P_lr <- if (testMask[['lr']]) stats::pchisq(q = survival::survdiff(survival::Surv(evtime) ~ group, rho = 0, data = dat_2gr)$chisq,
                        df = 1L, lower.tail = FALSE)
  # Peto & Peto modified Gehan-Wilcoxon test
  P_lr_pp <- if (testMask[['lr_pp']]) stats::pchisq(q = survival::survdiff(survival::Surv(evtime) ~ group, rho = 1, data = dat_2gr)$chisq,
                           df = 1L, lower.tail = FALSE)


  structure(
    purrr::compact(list(
      #fit0 = fit0, fit1 = fit1, ##for debugging only
      t_obs = ts_obs[["val"]],
      testDist = t0_dist,
      R = if (testMask[['bootstrap']]) length(t0_dist),
      chisq_df_hat = chisq_df_hat,
      param = param,
      P = purrr::compact( # compact = cleanse NULL entries
        list(
          bootstrap = P_boot,
          moran = as.vector(GOF_mo0$p.value),
          moran1 = as.vector(GOF_mo1$p.value),
          pearson = as.vector(GOF_pears0$p.value),
          pearson1 = as.vector(GOF_pears1$p.value),
          lr = P_lr,
          lr_pp = P_lr_pp
        )
      )
    )), class = "incubate_test")
}

#' @export
print.incubate_test <- function(x, ...){
  params <- paste(x$param, collapse = ' & ')
  P_boot_str <- if (is.numeric(x$P$bootstrap)) format.pval(x$P$bootstrap) else '-'
  cat(glue("Test for difference in {if (length(x$param) > 1) 'parameters' else 'parameter'} {params} between two groups.",
           "Alternative hypothesis: {params} {c('is', 'are')[[1L + (length(x$param) > 1)]]} different between the two groups.",
           "Parametric Bootstrap P-value: {P_boot_str}", .sep = "\n"), '\n')
}


#' @export
plot.incubate_test <- function(x, y, title, subtitle, ...){
  stopifnot(inherits(x, "incubate_test"))

  rlang::check_installed(pkg = 'ggplot2', reason = 'to get plots', version = '3.3')

  testDist <- x[["testDist"]]

  if (is.null(testDist)) {
    message('No bootstrap test to visualize.')
    return(invisible(NULL))
  }

  R <- purrr::pluck(x, "R", .default = "-")

  if (missing(title)) title <- glue("Distribution of test statistic under H0 for parameter {paste(x$param, collapse = ' & ')}")
  if (missing(subtitle) && ! is.null(x$P$bootstrap)){
    subtitle <- glue('Sampling distribution, based on {R} parametric bootstrap draws. ',
                     'Bootstrap P-value = {format.pval(x$P$bootstrap, eps = 1e-3)}')
    #"Approximated by a chi-square distribution with df={signif(x[['chisq_df_hat']], 2)}.")
  }


  p <- ggplot2::ggplot(tibble::tibble(testDist = testDist),
                       mapping = ggplot2::aes_(x = ~testDist, y = ~ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(bins = 11L + ceiling(sqrt(R)))

  # extract maximum density value
  ymax <- max(ggplot2::layer_data(p)[['y']])
  ymax <- ceiling(max(ymax + .1, ymax * 1.01))

  if (! is.null(x[['chisq_df_hat']]) && is.numeric(x[['chisq_df_hat']])) p <- p + ggplot2::geom_function(inherit.aes = FALSE,
                                                                                                         fun = stats::dchisq, args = list(df = x[['chisq_df_hat']]),
                                                                                                         col = "red", linetype = "dotted")

  p +
    ggplot2::geom_vline(xintercept = x[["t_obs"]], linetype = "dashed", colour = "grey") +
    ggplot2::coord_cartesian(ylim = c(0L, ymax)) +
    ggplot2::labs(x = "Test statistic",
                  title = title, subtitle = subtitle)

}


#' Power simulation function for a two-group comparison of the delay parameter.
#'
#' There are two ways of operation:
#' 1. `power=NULL` Given sample size `n` it simulates the power.
#' 2. `n=NULL` Given a power an iterative search is started to find a suitable `n` within a specified range.
#'
#' In any case, the distribution, the parameters that are tested for, the type of test and the effect size (`eff=`) need to be specified.
#' The more power simulation rounds (parameter `nPowerSim=`) the more densely the space of data according to the specified model is sampled.
#'
#' Note that this second modus (when `n` is estimated) is computationally quite heavy.
#' The iterative search for `n` uses some heuristics and the estimated sample size might actually give a different power-level.
#' It is important to check the stated power in the output. The search algorithm comes to results closer to the power aimed at
#' when the admissible range for sample size (`nRange=`) is chosen sensibly.
#' In case the estimated sample size and the achieved power is too high it might pay off to rerun the function with an adapted admissible range.
#'
#' @param distribution character. Which assumed distribution is used for the power calculation.
#' @param twoPhase logical(1). Do we model two phases per group? Default is `FALSE`, i.e. a single delay phase per group.
#' @param eff list. The two list elements contain the model parameters (as understood by the delay-distribution functions provided by this package) for the two groups.
#' @param param character. Parameter name(s) which are to be tested for difference and for which to simulate the power. Default value is `'delay1'`.
#' @param test character. Which test to use for this power estimation?
#' @param n integer. Number of observations per group for the power simulation or `NULL` when n is to be estimated for a given power.
#' @param power numeric. `NULL` when power is to be estimated for a given sample size or a desired power is specified (and `n` is estimated).
#' @param r numeric. Ratio of both groups sizes, ny / nx. Default value is 1, i.e., balanced group sizes. Must be positive.
#' @param sig.level numeric. Significance level. Default is 0.05.
#' @param nPowerSim integer. Number of simulation rounds. Default value 1600 yields a standard error of 0.01 for power if the true power is 80%.
#' @param R integer. Number of bootstrap samples for test of difference in parameter within each power simulation. It affects the resolution of the P-value for each simulation round. A value of around `R=200` gives a resolution of 0.5% which might be enough for power analysis.
#' @param nRange integer. Admissible range for sample size when power is pre-specified and sample size is requested.
#' @return List of results of power simulation. Or `NULL` in case of errors.
#' @export
power_diff <- function(distribution = c("exponential", "weibull"), twoPhase = FALSE, param = "delay1",
                       test = c('bootstrap', 'pearson', 'moran', 'lr', 'lr_pp'),
                       eff = stop("Provide parameters for both groups that reflect the effect!"),
                       n = NULL, r = 1, sig.level = 0.05, power = NULL, nPowerSim = 1600, R = 201,
                       nRange = c(5, 50)){

  tol_pow <- .001
  distribution <- match.arg(distribution)
  if (! missing(test)) test <- tolower(test)
  test <- match.arg(arg = test)
  ranFun <- getDist(distribution, type = "r")
  onames <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, transformed = FALSE)

  # parameters to test differences and for which power is requested
  if (any(grepl(pattern = "_tr", param, fixed = TRUE))){
    stop("Parameter names in param= refer to the distribution parameters and not to the transformed parameters of the objective function.", call. = FALSE)
  }

  # translate convenience names (for single phase) to canonical names
  unNmbrdIdx <- !grepl(pattern = "[12]", param, fixed = FALSE)
  if (any(unNmbrdIdx)){
    param[unNmbrdIdx] <- paste0(param[unNmbrdIdx], "1") #interpret un-numbered parameters as referring to phase 1
    if (verbose > 0L) cat("Unnumbered parameter names in param= are taken to refer to initial phase and are translated to canonical parameter names.\n")
  }

  # only valid names in canonical order
  param <- intersect(onames, param)

  if (!length(param)){
    stop("Provide valid parameter names from the distribution to test for differences in two groups.", call. = FALSE)
  }
  param <- match.arg(param, choices = onames)

  stopifnot( is.null(n) || (is.numeric(n) && length(n) == 1L && is.finite(n) ))
  stopifnot( is.null(power) || (is.numeric(power) && length(power) == 1L && power > 0L && power < 1L ))
  if (is.null(n) + is.null(power) != 1L) stop('Either set `n=NULL` or `power=NULL`!')

  stopifnot( length(sig.level) == 1L, is.numeric(sig.level), is.finite(sig.level), sig.level > 0L, sig.level < 1L )
  stopifnot(is.numeric(r), length(r) == 1L, r > 0L)
  stopifnot( length(nPowerSim) == 1L, is.numeric(nPowerSim), nPowerSim >= 3L )
  stopifnot( length(R) == 1L, is.numeric(R), R >= 3L )
  stopifnot( length(nRange) == 2L, is.numeric(nRange), nRange[[1L]] > 1L, nRange[[2L]] > nRange[[1L]] )
  nPowerSim <- ceiling(nPowerSim)
  R <- ceiling(R)

  stopifnot( is.list(eff), length(eff) == 2L )
  parx <- eff[[1L]]
  pary <- eff[[2L]]

  stopifnot( is.numeric(parx), is.numeric(pary) )
  stopifnot( length(parx) == length(onames), length(pary) == length(onames))
  parx <- purrr::set_names(parx, onames)
  pary <- purrr::set_names(pary, onames)


  simulatePower <- function(nx, ny, B = nPowerSim, R){
    nx <- ceiling(nx); ny <- ceiling(ny)

    # repeatedly test for difference in parameter on bootstrapped data
    P_dist <- future.apply::future_vapply(X = seq_len(B), FUN.VALUE = double(1L),
                                          FUN = function(dummy) {
                                            # generate data according to chosen model
                                            #+and with the specified effect
                                            datx <- purrr::exec(ranFun, !!! c(n=nx, parx))
                                            daty <- purrr::exec(ranFun, !!! c(n=ny, pary))

                                            P_val <- NA_real_
                                            try(expr = {
                                              P_val <- purrr::pluck(test_diff(x = datx, y = daty,
                                                                              distribution = distribution, twoPhase = twoPhase,
                                                                              param = param, type = test, R = R),
                                                                    "P", test, .default = NA_real_)
                                            }, silent = TRUE)
                                            P_val
                                          }, future.seed = TRUE)


    P_dist <- P_dist[is.finite(P_dist)]

    if ( !length(P_dist) ){
      warning("No valid power simulation results.")
      return(invisible(NULL))
    }

    if ( length(P_dist) < 100L )
      warning("Low resultion for power estimate.")

    if (length(P_dist))
      sum(P_dist < sig.level) / length(P_dist) else
        NA_real_
  } #fun


  nx <- ny <- -1
  powerGrid <- NULL

  if (is.null(power)){

    # easy case: estimate power once!
    nx <- ceiling(n)
    ny <- ceiling(r * n)
    if ( nx < length(onames) || ny < length(onames) ){
      warning("Too few observations to fit parameters.")
      return(invisible(NULL))
    }

    power <- simulatePower(nx = nx, ny = ny, B = nPowerSim, R = R)

  } else {

    # estimate n for specified power
    stopifnot( is.null(n) )

    Bmax1 <- 200L
    Rmax1 <- 100L
    B1 <- min(Bmax1, nPowerSim)
    R1 <- min(Rmax1, R)
    i2 <- -1L

    # 1st iteration
    nx_cand1 <- unique(ceiling(seq.int(from = nRange[[1L]], to = nRange[[2L]], length.out = 5L)))
    NBR_CAND1 <- length(nx_cand1)
    pow_cand1 <- rep_len(-1, length.out = NBR_CAND1)

    # if single n remains, return the power for it (no search for n necessary)
    if (NBR_CAND1 == 1L) return(power_diff(distribution, twoPhase = twoPhase, param, test, eff,
                                           n = nx_cand1[[1L]], power = NULL,
                                           r = r, sig.level = sig.level,nPowerSim = nPowerSim, R = R))


    for (i1 in seq_along(nx_cand1)) {
      nxc <- nx_cand1[[i1]]
      pow_cand1[[i1]] <- simulatePower(nx = nxc, ny = nxc * r, B = B1, R = R1)

      if (pow_cand1[[i1]] >= power - tol_pow) break
    } #rof

    # store preliminary power estimates
    powerGrid <- tibble::tibble(
      nx = nx_cand1[pow_cand1 > 0],
      ny = ceiling(nx * r),
      power = pow_cand1[pow_cand1 > 0],
      iter = 1L,
      B = B1,
      R = R1
    )

    if (NROW(powerGrid) <= 1L) {
      stop('Failed to find power estimates within specified range!', call. = FALSE)
    }

    REFINE <- TRUE #NROW(powerGrid) >= 2L

    # check first iteration
    if (i1 == 1L){
      warning('Smallest allowed n already exceeds requested power!', call. = FALSE)
      REFINE <- FALSE
    }

    # check last iteration
    if (i1 == NBR_CAND1 && pow_cand1[[NBR_CAND1]] > -1 && pow_cand1[[NBR_CAND1]] < power - tol_pow){
      warning(glue('Failed to reach requested power with maximally allowed n: {nx_cand1[[NBR_CAND1]]} ',
                   'yields a power of {as_percent(pow_cand1[[NBR_CAND1]])}.'))
      REFINE <- FALSE
    }


    if (!REFINE){
      nx <- ceiling(nx_cand1[[i1]])
      ny <- ceiling(nx_cand1[[i1]] * r)
      power <- if (B1 < nPowerSim || R1 < R) stats::weighted.mean(x = c(pow_cand1[[i1]], simulatePower(nx, ny, B = nPowerSim, R = R)),
                                                                  w = c(B1, nPowerSim)) else
                                                                    pow_cand1[[i1]]
    } else {
      powerMod <- if (NROW(powerGrid) == 2L) stats::lm(power ~ nx, data = powerGrid) else stats::lm(power ~ poly(nx, 2), data = powerGrid)
      powerPred <- tibble::tibble(nx = seq.int(from = nRange[[1L]], to = nRange[[2L]], by = 1L),
                                  predpower = stats::predict.lm(powerMod, newdata = data.frame(nx = nx)),
                                  diffpower = .data$predpower - power)
      # examine close neighbourhood of predicted best n
      powerPredInd <- intersect(seq_len(NROW(powerPred)), c(-1L, 0L, 1L) + which.max(powerPred$diffpower >= 0L))


      nx_cand2 <- powerPred$nx[powerPredInd]
      NBR_CAND2 <- length(nx_cand2)
      pow_cand2 <- rep_len(-1, length.out = NBR_CAND2)

      for (i2 in seq_along(nx_cand2)){
        nxc <- nx_cand2[[i2]]
        pow_cand2[[i2]] <- simulatePower(nx = nxc, ny = nxc * r, B = nPowerSim, R = R)

      } #rof

      powerGrid2 <- tibble::tibble(
        nx = nx_cand2[pow_cand2 > 0],
        ny = ceiling(nx * r),
        power = pow_cand2[pow_cand2 > 0],
        iter = 2L,
        B = nPowerSim,
        R = R)

      stopifnot( any(powerGrid2$power >= power - tol_pow) )
      nx <- powerGrid2$nx[which.max(powerGrid2$power >= power - tol_pow)]
      ny <- ceiling(nx * r)
      power <- powerGrid2$power[which(powerGrid2$nx == nx)]

      # store 2nd round (refinement) power estimates
      powerGrid <- rbind(powerGrid, powerGrid2)

    }
    stopifnot( nx > 0L, ny > 0L, power > 0L)

  }

  purrr::compact(
    list(name = "Difference in delayed model for time-to-event data in two groups",
         distribution = distribution,
         twoPhase = twoPhase,
         param = param,
         test = test,
         eff = eff,
         sig.level = sig.level,
         nx = nx, ny = ny, N = nx + ny,
         #P_dist = P_dist, ##debug
         powerGrid = powerGrid,
         power = power
    )
  )
}

