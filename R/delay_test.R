
#' Goodness-of-fit (GOF) test statistic.
#' The GOF-test is performed for a fitted delay-model.
#' There are different GOF-tests implemented:
#' * __Moran GOF__ is based on spacings, like the MSE-criterion itself.
#' * __Pearson GOF__ uses categories and compares observed to expected frequencies.
#' * __Anderson-Darling GOF__ makes inference based on the cumulative distribution function.
#'
#' @param delayFit delay_model fit
#' @param method character(1). which method to use for GOF. Default is 'moran'.
#' @return An `htest`-object containing the GOF-test result
#' @export
test_GOF <- function(delayFit, method = c('moran', 'pearson', 'AD', 'ad', 'anderson')){

  stopifnot( inherits(delayFit, what = 'incubate_fit') )
  method <- match.arg(method)

  if (delayFit$method != 'MSE'){
    stop('Goodness-of-fit test only supported for models that are fit with maximum spacing estimation (MSE)!')
  }

  twoGr <- isTRUE(delayFit$twoGroup)
  data_name <- if (twoGr) paste(names(delayFit$data), collapse = ' and ') else 'x'
  nObs <- if (twoGr) lengths(delayFit$data) else length(delayFit$data)
  params <- coef(delayFit)
  k <- length(params)


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


           statist <- if (twoGr){ ##  && length(delayFit$bind) < length(oNames) # not needed!?
             # sum of two independent chi-sq. is chi-sq
             c(`X^2` = sum(testStat_mo(mseCrit = delayFit$objFun(pars = params, aggregated = FALSE), #criterion per group
                                       n = nObs, k = k/2)) )
           } else {
             # single group
             c(`X^2` = testStat_mo(mseCrit = delayFit$val,
                                   n = nObs, k = k) )
           }

           p_val <- stats::pchisq(q = statist, df = sum(nObs), lower.tail = FALSE)
         },

         pearson = {
           # Pearson GOF-test
           meth <- "Pearson's Goodness-of-fit (GOF) test" # (per group) ## this is our standard

           # under H0, expect frequency counts according to uniform distribution
           #+nbr of classes as recommended by David S. Moore (Tests of Chi-squared Type, 1986)

           testStat_pe <- function(datr, nCl){
             tab_transf <- tabulate(findInterval(datr,
                                   vec = seq.int(from=0L, to=1L, length.out = nCl+1L),
                                   rightmost.closed = TRUE, all.inside = TRUE), nbins = nCl)
             c(`X^2` = sum((tab_transf - mean(tab_transf))**2L) / mean(tab_transf))
           }

           if (twoGr) {
             nCl <- pmax.int(k/2 + 2L, ceiling(2L * nObs**.4))
             # sum of two chi-square test statistics
             statist <- c(`X^2` = sum(purrr::map2_dbl(.x = delayFit$data_transf, .y = nCl, .f = testStat_pe)))
             p_val <- stats::pchisq(q = statist, df = sum(nCl) - k - 2L, lower.tail = FALSE)
           } else {
             nCl <- max(k + 2L, ceiling(2L * nObs**.4))
             statist <- testStat_pe(datr = delayFit$data_transf, nCl = nCl)
             # inference based on Chi-square distribution.
             # use adjusted degrees of freedom (loose one df for each parameter estimated)
             p_val <- stats::pchisq(q = statist, df = nCl - k - 1L, lower.tail = FALSE)
           }

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

           A2 <- if (twoGr) {
             purrr::map2_dbl(.x = delayFit$data_transf, .y = nObs, .f = testStat_ad)
           } else {
             testStat_ad(datr = delayFit$data_transf, n = nObs)
           }

           p_val <- if ( identical(delayFit$distribution, 'exponential') ){
             # modification for Exponential (cf Stephens, Table 4.14, p.138)
             # the correction factor approaches 1 from above.
             # We keep the number of N as all observations, independent of the number of parameters estimated in the null-model.
             # QQQ Should we increase N by the number of parameters p estimated less 2 ( p -2 because 2 parameters are estimated in standard delayed exponential)
             A2_mod <- A2 * pmax.int(1L, 1L + 5.4 / nObs - 11 / nObs**2L)

             .ad_pval[['exponential']](A2_mod)

           } else {
             stopifnot( identical(delayFit$distribution, 'weibull') )
             # P-value for Weibull based on Lockhart, 1994 (Table 1)
             # interpolation model on logits using critical value and inverse of shape parameter
             .ad_pval[['weibull']](A2, params[grepl('shape', names(params), fixed = TRUE)])
           }

           if (twoGr){
             A2 <- paste(signif(A2, 4), collapse = ' and ')
             # use Liptak to bring both P-values of AD-tests per group together
             p_val <- pnorm(sum(sqrt(nObs) * qnorm(p_val)) / sqrt(sum(nObs)))
           }

           statist <- c(`A^2` = A2)

         },
         # catch all
         stop('This method is not supported!')
  )



  # return test object

  # stats:::print.htest recognizes:
  #+parameter, alternative, null.value, conf.int, estimate
  structure(
    list(method = meth, data.name = data_name,
         statistic = statist, p.value = as.numeric(p_val)),
    class = 'htest')
}


#' Test the difference for delay model parameter(s) between two uncorrelated groups, based on maximum spacing estimation (MSE).
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
#' @param distribution character[1]. Name of the parametric delay distribution to use.
#' @param param character. Names of parameters to test difference for. Default value is `'delay'`.
#' @param R numeric[1]. Number of bootstrap samples to evaluate the distribution of the test statistic.
#' @param ties character. How to handle ties in data vector of a group?
#' @param verbose numeric. How many details are requested?
#' @return list with the results of the test. Element P contains the different P-values, for instance from parametric bootstrap
#' @export
test_diff <- function(x, y=stop('Provide data for group y!'), distribution = c("exponential", "weibull"), param = "delay", R = 400,
                      ties = c('equidist', 'density', 'random'), verbose = 0) {
  STRICT <- TRUE #keep only conv=0 model fits?
  distribution <- match.arg(distribution)
  ties <- match.arg(ties)
  par_names <- getDist(distribution = distribution, type = "param")
  stopifnot( is.numeric(x), length(x) > length(par_names), is.numeric(y), length(y) > length(par_names) )
  stopifnot( is.numeric(R), length(R) == 1L, R >= 1L )
  stopifnot( is.character(param), length(param) >= 1L, nzchar(param) )

  # Test statistic calculated from the given data and the model specification.
  #
  # The test statistic takes non-negative values.
  # High values of the test statistic speak in favour of H1:
  # @return list containing value of test statistic and null model fit. Or `NULL` in case of trouble.
  testStat <- function(x, y) {
    #fit0 <- fit1 <- NULL
    fit0 <- delay_model(x = x, y = y, distribution = distribution, method = 'MSE', ties = ties, bind = param)
    fit1 <- delay_model(x = x, y = y, distribution = distribution, method = 'MSE', ties = ties)

    if (is.null(fit0) || is.null(fit1)) return(invisible(NULL))

    # if more restricted model fit0 yields better fit (=lower criterion) than more general fit1
    #+we are in trouble, possibly due to non-convergence, e.g. convergence code 52
    #+we refit the general fit1 again using parameter-values from fit0
    if ( fit0[['val']] < fit1[['val']] ){
      warning('Restricted model with better fit than unrestricted model.', call. = FALSE)
      # re-run fit1 with start values based on fitted parameters of reduced model fit0
      fit1oa <- attr(fit1[['objFun']], 'optim_args', exact = TRUE)
      stopifnot( is.list(fit1oa), 'par' %in% names(fit1oa) )

      # QQQ does match() or pmatch help here to avoid for-loop?
      pn1 <- names(fit1[['par']])
      for (na0 in names(fit0[['par']])) fit1oa[['par']][startsWith(pn1, prefix = na0)] <- coef(fit0)[[na0]]

      ## XXX abs should not be necessary because we allow only for non-negative parameters
      # Better apply lower limits for parameters?!!
      newparsc <- abs(fit1oa[['par']])
      newparsc[which(newparsc < 1e-7)] <- 1e-7
      newparsc[which(newparsc > 1e8)] <- 1e8
      fit1oa[['control']][['parscale']] <- newparsc

      fit1 <- update(fit1, optim_args = fit1oa)

      if (is.null(fit1) || fit0[['val']] < fit1[['val']]) return(invisible(NULL))
    }# fi bad fit1

    # convergence of re-fits?
    # keep also fits with error-code 52: in my tests all those fits *looked* actually OK..
    # XXX try harder for non-convergence when testStat is called first?
    if ( STRICT && (fit0[['convergence']] != 0 || fit1[['convergence']] != 0) ) return(invisible(NULL))

    # higher values of T speak in favour of H1:
    #   1. fit0 has high value (=bad fit)
    #   2. fit1 has low value (=good fit)
    list(val = 2L * (fit0[["val"]] - fit1[["val"]]),
         fit0 = fit0)
  }

  # observed test statistic
  ts_obs <- testStat(x, y)
  if (is.null(ts_obs)){
    stop("Delay model failed for restricted null-model or free full model", call. = FALSE)
  }
  fit0 <- ts_obs[["fit0"]]


  # P-values based GOF-tests
  #+ H0: simpler/restricted model 0 is sufficient
  #+ the GOF-test solely builds on fit0
  #+ take fitted parameters for both groups under null-model
  #+ and transform the observed data for both groups via cumulative distribution functions

  # spacings-based GOF-test
  GOF_mo <- test_GOF(delayFit = fit0, method = 'moran')
  #if (verbose > 0L) cat("Moran test stat: ", GOF_mo$statistic, "\n")

  # Pearson GOF-test based on Chi-square distribution.
  # under H0, expect counts according to uniform distribution
  GOF_pears <- test_GOF(delayFit = fit0, method = 'pearson')

  # EDF-based GOF-test
  # Anderson-Darling (AD) test statistic, cf Stephens, Tests based on EDF Statistics p.101, (4.2)
  GOF_ad <- test_GOF(delayFit = fit0, method = 'AD')


  # parametric bootstrap:
  # generate R samples (x, y) by random sampling from the fitted H0-model (e.g. common delay),
  #+where all nuisance parameters are at there best fit
  # calculate the test statistic on the simulated data
  # estimate P as proportion of simulated test statistics that exceed the observed test statistic t_obs

  ranFun <- getDist(distribution, type = "r")
  # arguments to the random function generation
  ranFunArgsX <- c(list(n=length(x)), coef(fit0, group = "x"))
  ranFunArgsY <- c(list(n=length(y)), coef(fit0, group = "y"))

  t0_dist <- future.apply::future_vapply(X = seq.int(R), FUN.VALUE = double(1L+(verbose>0L)),
                                         FUN = function(dummy){

                                           # generate new data according to given fitted null-model
                                           # sort is not needed here, as it goes through the whole pipeline (factory method)
                                           ts_boot <- testStat(x = rlang::exec(ranFun, !!! ranFunArgsX),
                                                               y = rlang::exec(ranFun, !!! ranFunArgsY))
                                           if (is.null(ts_boot)) return(rep(NA_real_, 1L+(verbose>0L)))

                                           if (verbose > 0L)
                                             c(ts_boot[["val"]], ts_boot[['fit0']][['convergence']]) else
                                               ts_boot[['val']]


                                         }, future.seed = TRUE)

  if (verbose > 0L){
    fit0_conv <- t0_dist[2L,]
    cat('Proportion of model failures:', sprintf('%6.1f%%', 100L*length(which(is.na(fit0_conv)))/length(fit0_conv)), '\n')
    cat('Proportion of convergence= 0:', sprintf('%6.1f%%', 100L*length(which(fit0_conv == 0))/length(fit0_conv)), '\n')
    cat('Proportion of convergence=52:', sprintf('%6.1f%%', 100L*length(which(fit0_conv == 52))/length(fit0_conv)), '\n')
    t0_dist <- t0_dist[1L,]
  }
  t0_dist <- t0_dist[is.finite(t0_dist)]
  chisq_df_hat <- NULL
  try(expr = {chisq_df_hat <- coef(MASS::fitdistr(x = t0_dist, densfun = "chi-squared",
                                                  start = list(df = length(param)),
                                                  method = "Brent", lower = .001, upper = 401))},
      silent = TRUE)

  P_boot <- (1L + sum(t0_dist >= ts_obs[["val"]])) / (length(t0_dist)+1L)


  ## P-value from Log-rank tests
  dat_2gr <- tibble::tibble(evtime = c(x,y),
                            group = rep(c("x", "y"), times = c(length(x), length(y))))
  P_lr <- stats::pchisq(q = survival::survdiff(survival::Surv(evtime) ~ group, rho = 0, data = dat_2gr)$chisq,
                        df = 1L, lower.tail = FALSE)
  # Peto & Peto modified Gehan-Wilcoxon test
  P_lr_pp <- stats::pchisq(q = survival::survdiff(survival::Surv(evtime) ~ group, rho = 1, data = dat_2gr)$chisq,
                           df = 1L, lower.tail = FALSE)


  structure(
    list(#fit0 = fit0, fit1 = fit1, ##debug
      t_obs = ts_obs[["val"]],
      testDist = t0_dist,
      R = length(t0_dist),
      chisq_df_hat = chisq_df_hat,
      param = param,
      P = list(boot = P_boot,
               gof_mo = as.vector(GOF_mo$p.value),
               gof_pearson = as.vector(GOF_pears$p.value),
               gof_ad = as.vector(GOF_ad$p.value),
               lr = P_lr,
               lr_pp = P_lr_pp)
    ), class = "incubate_test")
}

#' @export
print.incubate_test <- function(x, ...){
  params <- paste(x$param, collapse = ' & ')
  cat('Test for difference in parameter ', params, 'between two groups.\n')
  cat('Alternative hypothesis: ', params, 'are different between the two groups.\n')
  cat('Parametric Bootstrap P-value: ', format.pval(x$P$boot))
  cat('\n')
}

#' @export
plot.incubate_test <- function(x, y, title, subtitle, ...){
  stopifnot(inherits(x, "incubate_test"))

  teststat <- x[["testDist"]]

  if (missing(title)) title <- glue::glue('Distribution of test statistic under H0 for parameter {dQuote(x$param)}')
  if (missing(subtitle)) subtitle <- glue::glue('Sampling distribution, based on {length(teststat)} parametric bootstrap draws. ',
                                                'Bootstrap P-value = {format.pval(x$P$boot, eps = 1e-3)}')
  #"Approximated by a chi-square distribution with df={signif(x[['chisq_df_hat']], 2)}.")


  p <- dplyr::tibble(teststat = teststat) %>%
    ggplot2::ggplot(ggplot2::aes(x = teststat, y = ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(bins = 11L + ceiling(sqrt(length(teststat))))

  # extract maximum density value
  ymax <- ggplot2::layer_data(p) %>%
    dplyr::pull(y) %>%
    {ceiling(max(. + .1, . * 1.01))}

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
#' Given an effect size and a sample size `n` it simulates the power.
#' The more power simulation rounds (parameter `nPowerSim=`) the more densely the space of data according to the specified model is sampled.
#'
#' @param eff list of length 2. Model parameters (as understood by the delay-distribution functions provided by this package) for each of the two groups.
#' @param param character. Parameter name for which to simulate the power.
#' @param n integer. Number of observations per group for the power simulation. Can be two different numbers, control group and then treatment group.
#' @param nPowerSim integer. Number of simulation rounds. Default value 1600 yields a standard error of 0.01 for power if the true power is 80%.
#' @param R integer. Number of bootstrap samples to assess difference in parameter within each power simulation. It affects the resolution of the P-value for each simulation round. A value of around `R=200` gives a resolution of 0.5% which is enough for power analysis.
#' @return List of results of power simulation. Or `NULL` in case of errors.
#' @export
power_diff <- function(distribution = c("exponential", "weibull"), param = "delay",
                       eff = stop("Provide parameters for each group (reference group first) that reflect the effect!"),
                       n, sig.level = 0.05, nPowerSim = 16e2, R = 2e2){

  distribution <- match.arg(distribution)
  ranFun <- getDist(distribution, type = "r")
  par_names <- getDist(distribution, type = "param")
  param <- match.arg(param, choices = par_names)

  stopifnot( is.list(eff), length(eff) == 2L )
  stopifnot( length(sig.level) == 1L, is.numeric(sig.level), is.finite(sig.level) )
  stopifnot( is.numeric(n), length(n) > 0L, length(n) <= 2L, all(is.finite(n) & n >= 1L) )
  stopifnot( length(nPowerSim) == 1L, is.numeric(nPowerSim), nPowerSim >= 3L )
  stopifnot( length(R) == 1L, is.numeric(R), R >= 5L )
  nPowerSim <- ceiling(nPowerSim)
  R <- ceiling(R)

  if (length(n) == 1L) n <- rep(n, 2L)

  n_ctrl <- n[[1L]]
  n_trtm <- n[[2L]]


  if (any(n < length(par_names))){
    warning("Too few observations to fit parameters.")
    return(invisible(NULL))
  }

  par_ctrl <- eff[[1L]]
  par_trtm <- eff[[2L]]

  stopifnot( is.numeric(par_ctrl), is.numeric(par_trtm) )
  stopifnot( length(par_ctrl) == length(par_names), length(par_trtm) == length(par_names))
  par_ctrl <- purrr::set_names(par_ctrl, par_names)
  par_trtm <- purrr::set_names(par_trtm, par_names)


  # repeatedly test for difference in parameter on bootstrapped data
  P_dist <- future.apply::future_vapply(X = seq(nPowerSim), FUN.VALUE = double(1L),
                                        FUN = function(dummy) {
                                          # generate data according to chosen model
                                          #+and with the specified effect
                                          dat_ctrl <- purrr::exec(ranFun, !!! c(n=n_ctrl, par_ctrl))
                                          dat_trtm <- purrr::exec(ranFun, !!! c(n=n_trtm, par_trtm))

                                          # if (is.null(fit_ctrl) || is.null(fit_trtm) ||
                                          #     fit_ctrl$convergence > 0L || fit_trtm$convergence > 0L ) NA else
                                          P_val <- NA_real_
                                          try(expr = {
                                            P_val <- test_diff(x = dat_ctrl, y =dat_trtm,
                                                               distribution = distribution, param = param, R = R) %>%
                                              purrr::pluck("P", "boot", .default = NA_real_) }, silent = TRUE)
                                          P_val
                                        }, future.seed = TRUE)


  P_dist <- P_dist[is.finite(P_dist)]

  if ( !length(P_dist) ){
    warning("No valid power simulation results.")
    return(invisible(NULL))
  }

  if ( length(P_dist) < 100L )
    warning("Low resultion for power estimate.")

  # structure(
  list(id = "delayed:2groups", name = "Difference in delayed model for time-to-event data in two groups",
       distribution = distribution, param = param,
       eff = eff, sig.level = sig.level, n = n, N = sum(n),
       P_dist = P_dist, ##debug
       power = if (length(P_dist)) sum(P_dist < sig.level) / length(P_dist) else NA_real_
  )
  #   , class = "sscn"
  # )
}

