# mkuhn, 2021-04-06
# MPS fitting
# Further ideas:
# 1/ how to deal with censoring? In particular, censoring before the first observed event
# 2/ currently, a single early event has dramatic consequences on the estimate of delay.
#+The delay parameter is assumed to be fixed for the whole population of subjects.
#+Could we think of an inter-individual variability of delay? Or: associate a variability of delay linked with the rate parameters?
#+i.e. we want to introduce a probabilistic element to the delay.. think for Weibull (Exp follows)


# distribution functions ----

#' Density of delayed exponential distribution.
#' @export
dexp_delayed <- function(x, delay, rate = 1, ...) dexp(x = x - delay, rate, ...)
#' @export
pexp_delayed <- function(q, delay, rate = 1, ...) pexp(q = q - delay, rate, ...)
#' @export
rexp_delayed <- function(n, delay, rate = 1) rexp(n = n, rate = rate) + delay

#' Density of delayed Weibull distribution
#' @export
dweib_delayed <- function(x, delay, shape, scale = 1, ...) dweibull(x = x - delay, shape, scale, ...)
#' @export
pweib_delayed <- function(q, delay, shape, scale = 1, ...) pweibull(q = q - delay, shape, scale, ...)
#' @export
rweib_delayed <- function(n, delay, shape, scale = 1) rweibull(n = n, shape = shape, scale = scale) + delay


#' Get delay distribution function
#' @param distribution character[1]. delay distribution.
#' @param type character[1]. type of function, cdf: cumulative distribution function, density or random function
#' @param twoGroup logical[1]. Do we have two groups?
#' @param bind character. Names of parameters that are bind between the two groups.
#' @return selected distribution function or parameter names
getDist <- function(distribution = c("exponential", "weibull"), type = c("cdf", "prob", "density", "random", "param"),
                    twoGroup = FALSE, bind = NULL) {
  distribution <- match.arg(distribution)
  type <- match.arg(type)

  switch(type,
         # cumulative distribution function
         prob = ,
         cdf     = c(pexp_delayed, pweib_delayed),
         # density function
         density = c(dexp_delayed, dweib_delayed),
         # random draw function
         random  = c(rexp_delayed, rweib_delayed),
         param   = {
           par_list <- list(exponential = c("delay", "rate"),
                            weibull = c("delay", "shape", "scale"))
           if (twoGroup) {
             # bind only parameters from chosen distribution
             myPars <- par_list[[1L + (distribution == 'weibull')]]
             bind <- intersect(myPars, bind) #intersect: enforces original order
             par_gr <- purrr::map(par_list, ~ setdiff(., bind))
             # bind parameters first
             purrr::map(par_gr, ~ c(bind,
                                    paste(rep(., 2L), rep(c("x", "y"), each = length(.)),
                                          sep = ".")))
           } else par_list

         },
         stop(glue::glue("Unknown attribute of distribution {distribution}."))
  )[[1L + (distribution == 'weibull')]]
}


# delay estimation ----

#' Factory method for negative maximum spacing estimation objective function.
#'
#' Negative or infinite values are discarded.
#' @param x numeric. observations
#' @param y numeric. observations in second group.
#' @param distribution character(1). delayed distribution family
#' @param bind character(1). parameter names that are bind together (i.e. equated) between both groups
#' @return objective function
geomSpaceFactory <- function(x, y=NULL, distribution = c("exponential", "weibull"), bind=NULL) {
  twoGr <- ! is.null(y)
  stopifnot( is.numeric(x), length(x) > 0L, ! twoGr || is.numeric(y) && length(y) > 0L )
  stopifnot( is.null(bind) || is.character(bind) && length(bind) >= 1L )
  # enforce that bind respects the canonical order of dist-parameters
  oNames <- getDist(distribution, type = "param", twoGroup = FALSE, bind = NULL)
  bind <- intersect(oNames, bind)

  ind_neg_x <- which(x < 0L)
  if (length(ind_neg_x)){
    warning("Negative values in data of x! These are dropped.")
    x <- x[-ind_neg_x]
  }

  # drop NA and ±Inf & sort
  x <- sort(x[is.finite(x)])

  if (!length(x)) {
    warning("No valid data in x! Only non-negative and finite real values are valid.")
    return(invisible(NULL))
  }

  if ( twoGr ){
    ind_neg_y <- which(y < 0L)
    if (length(ind_neg_y)){
      warning("Negative values in data of y! These are dropped.")
      y <- y[-ind_neg_y]
    }

    # drop NA and ±Inf & sort
    y <- sort(y[is.finite(y)])

    if (!length(y)) {
      warning("No valid data in y! Only non-negative and finite real values are valid.")
      return(invisible(NULL))
    }
  } #fi twoGr

  distribution <- match.arg(distribution)

  par_start <- NULL
  par_names <- getDist(distribution = distribution, type = "param",
                       twoGroup = twoGr, bind = bind)

  optim_args <- NULL

  #XXX par_start with twoGr & bind: only qucik*dirty solution
  #+ see function getPar_weib below or come up with something like:
  # getStartVals <- function(twoGr, bind){
  #   if (distribution == 'exponential')
  # }
  if (distribution == 'exponential') {
    # upper bound for *D*elay
    upD.x <- max(1e-6, min(x)-1L/length(x), min(x)*.9999)
    if (twoGr){

      par_start <- if (is.null(bind)) {
        c(max(0L, min(x)-2L/length(x)), 1L/mean(x - min(x) + 2L/length(x)),
          max(0L, min(y)-2L/length(y)), 1L/mean(y - min(y) + 2L/length(y)))
      } else if (identical(bind, 'delay')) {
        c(max(0L, min(min(x)-2L/length(x), min(y)-2L/length(y))),
          1L/mean(x - min(x) + 2L/length(x)),
          1L/mean(y - min(y) + 2L/length(y)))
      } else if (identical(bind, 'rate')) {
        c((1L/mean(x - min(x) + 2L/length(x)) + 1L/mean(y - min(y) + 2L/length(y)))/2L,
          max(0L, min(x)-2L/length(x)), max(0L, min(y)-2L/length(y)))
      } else {
        stopifnot( identical(bind, c("delay", "rate")) )
        c( max(0L, min(c(x, y))-2L/(length(x)+length(y))),
           1L/mean(c(x, y) - min(c(x, y)) + 2L/(length(x)+length(y))) )
      }

      myLower <- if (is.null(bind)) c(0L, 1e-9, 0L, 1e-9) else if (identical(bind, 'delay')) c(0L, 1e-9, 1e-9) else
        if (identical(bind, 'rate')) c(1e-9, 0L, 0L) else c(0L, 1e-9)

      upD.y <- max(min(y)-1L/length(y), min(y)*.9999)
      myUpper <- if (is.null(bind)) c(upD.x, +Inf, upD.y, +Inf) else
        if (identical(bind, 'delay')) c(min(upD.x, upD.y), +Inf, +Inf) else
          if (identical(bind, 'rate')) c(+Inf, upD.x, upD.y) else c(min(upD.x, upD.y), +Inf)

      optim_args <- list(par = par_start,
                         method = "L-BFGS-B",
                         # QQQ something like exp(trunc(log(par))) where par is the start parameters
                         control = list(parscale = pmin(1e7, pmax(1e-7, par_start))),
                         lower = purrr::set_names(myLower, par_names),
                         upper = purrr::set_names(myUpper, par_names))

    } else { # single group
      par_start <- c(max(0L, min(x)-2L/length(x)), 1L/mean(x - min(x) + 2L/length(x)))
      optim_args <- list(par = par_start,
                         method = "L-BFGS-B",
                         # QQQ something like exp(trunc(log(par))) where par is the start parameters
                         control = list(parscale = pmin(1e7, pmax(1e-7, par_start))),
                         lower = purrr::set_names(c(0L, 1e-9), par_names),
                         upper = purrr::set_names(c(upD.x, +Inf), par_names)
      )
    }

  } else {
    stopifnot( distribution == 'weibull' )

    # start values from 'Weibull plot'
    #+using the empirical distribution function
    ## in MASS::fitdistr they simplify:
    # lx <- log(x)
    # m <- mean(lx)
    # v <- var(lx)
    # shape <- 1.2/sqrt(v)
    # scale <- exp(m + 0.572/shape)
    getPar_weib <- function(obs){
      # contract: obs is sorted!
      # take out extreme values for robustness (against possible outliers)
      #+ when at least 10 observations
      xx <- obs[floor(length(obs)*.09):ceiling(length(obs)*.91)] #assume sorted data
      start_y <- log(-log(1-(seq_along(xx)-.3)/(length(xx)+.4)))
      # cf. lm.fit(x = cbind(1, log(obs)), y = start_y))$coefficients
      start_shape <- cor(log(xx), start_y) * sd(start_y) / sd(log(xx))
      start_scale <- exp(mean(log(xx) - mean(start_y) / start_shape))

      list(
        start = c(max(0L, min(obs) - 2L/length(obs)), start_shape, start_scale),
        lower = c(0L, 1e-9, 1e-9),
        upper = c(max(min(obs)-1L/length(obs), min(obs)*.9999), +Inf, +Inf)
        )
    }

    par0_x <- getPar_weib(x)
    if (twoGr){
      par0_y <- getPar_weib(y)
      start0_x <- par0_x[["start"]]
      start0_y <- par0_y[["start"]]

      par_start <- if (is.null(bind)) c(start0_x, start0_y) else
        if (identical(bind, 'delay')) c(min(start0_x[1L], start0_y[1L]), start0_x[-1L], start0_y[-1L]) else
          stop("Currently, for Weibull only bind=NULL and bind='delay' are supported!")

      lower_x <- par0_x[["lower"]]
      lower_y <- par0_y[["lower"]]
      myLower <- if (is.null(bind)) c(lower_x, lower_y) else
        if (identical(bind, 'delay')) c(min(lower_x[1L], lower_y[1L]), lower_x[-1L], lower_y[-1L]) else
            stop("Currently, for Weibull only bind=NULL and bind='delay' are supported!")

      upper_x <- par0_x[["upper"]]
      upper_y <- par0_y[["upper"]]
      myUpper <- if (is.null(bind)) c(upper_x, upper_y) else
        if (identical(bind, 'delay')) c(max(upper_x[1L], upper_y[1L]), upper_x[-1L], upper_y[-1L]) else
          stop("Currently, for Weibull only bind=NULL and bind='delay' are supported!")

      optim_args <- list(par = par_start,
                         method = "L-BFGS-B",
                         # QQQ something like exp(trunc(log(par))) where par is the start parameters
                         control = list(parscale = pmin(1e7, pmax(1e-7, par_start))),
                         lower = purrr::set_names(myLower, par_names),
                         upper = purrr::set_names(myUpper, par_names)
      )

    } else {  # single group

      par_start <- par0_x[["start"]]
      optim_args <- list(par = par_start,
                         method = "L-BFGS-B",
                         # QQQ something like exp(trunc(log(par))) where par is the start parameters
                         control = list(parscale = pmin(1e7, pmax(1e-7, par_start))),
                         lower = purrr::set_names(par0_x[["lower"]], par_names),
                         upper = purrr::set_names(par0_x[["upper"]], par_names) )
    }
  } # weibull

  # extract the parameters for the specified group
  # cf coef.mps_fit (below)
  # directly uses twoGr, oNames, bind from the factory
  getPars <- function(parL, group = "x") {
    if (! twoGr) return(parL)

    # extract all group parameters and restore original name (e.g. remove ".x")
    parL.gr <- purrr::set_names(parL[grepl(pattern = paste0(".", group), x = names(parL), fixed = TRUE)],
                                nm = setdiff(oNames, bind))

    # contract: bind is intersected and has right order
    # contract: bind comes first in par
    c(parL.gr, parL[bind])[oNames]
  }

  #' negative maximum spacing estimation objective function.
  #' Estimate parameters by minimizing this function.
  #' @param par parameter vector.
  #' XXX Here, all parameters are being optimized.
  #' +Make this more flexible: allow to set parameters to given fixed values
  negMSE <- function(par){
    stopifnot( length(par_names) == length(par) )
    parL <- purrr::set_names(as.list(par), par_names)

    # contract: x is sorted!
    cumDiffs.x <- diff(c(0L,
                         purrr::exec(getDist(distribution, type = "cdf"),
                                     !!! c(list(q=x), getPars(parL, group = "x"))),
                         1L))

    #XXX Does «mini» uniform distribution within rounding margin work as well for ties?
    #+ this would be easier for automatic differentiation (AD) and we would not need densFun at all!?

    # use densFun for ties
    ind_zx <- which(cumDiffs.x == 0L)
    if (length(ind_zx))
      cumDiffs.x[ind_zx] <- purrr::exec(getDist(distribution, type = "dens"),
                                        !!! c(list(x = x[pmax(ind_zx-1L, 1L)]), getPars(parL, group = "x")))

    # respect the machine's numerical lower limit
    cumDiffs.x[cumDiffs.x < .Machine$double.xmin] <- .Machine$double.xmin

    if (twoGr){
      # contract: y is sorted!
      cumDiffs.y <- diff(c(0L,
                           purrr::exec(getDist(distribution, type = "cdf"),
                                       !!! c(list(q=y), getPars(parL, group = "y"))),
                           1L))

      #XXX Does «mini» uniform distribution within rounding margin work as well for ties?
      #+ this would be easier for automatic differentiation (AD) and we would not need densFun at all!?

      # use densFun for ties
      ind_zy <- which(cumDiffs.y == 0L)
      if (length(ind_zy))
        cumDiffs.y[ind_zy] <- purrr::exec(getDist(distribution, type = "dens"), !!! c(list(x = y[pmax(ind_zy-1L, 1L)]),
                                                                                      getPars(parL, group = "y")))

      # respect the machine's numerical lower limit
      cumDiffs.y[cumDiffs.y < .Machine$double.xmin] <- .Machine$double.xmin
    }

    #XXX for twoGr:
    #+does it make a difference to first merge x and y and then do the cumDiffs, log and mean
    - if (twoGr) weighted.mean(c(mean(log(cumDiffs.x)), mean(log(cumDiffs.y))),
                               w = c(length(x), length(y))) else
                                 mean(log(cumDiffs.x))
  }

  # add "optim_args" & distribution as attributes to the objective function
  attr(negMSE, which = "optim_args") <- c(list(fn = negMSE), optim_args) #optim_args
  attr(negMSE, which = "distribution") <- distribution
  attr(negMSE, which = "twoGroup") <- twoGr
  attr(negMSE, which = "bind") <- bind

  negMSE
}


#' Parameter fitting through numerical optimization.
#'
#' The objective function carries the given data in its environment.
#' `stats::optim` does the optimization.
#' @param objFun objective function
#' @param optim_args list of own arguments for optimization. If `NULL` it uses the default optim arguments associated to the objective function.
#' @param verbose integer that indicates the level of verboseness. Default 0 is quiet.
#' @return optimization object or `NULL` in case of errors during optimization
delay_fit <- function(objFun, optim_args = NULL, verbose = 0) {

  if (is.null(objFun)) return(invisible(NULL))
  stopifnot( is.function(objFun) )
  distribution <- attr(objFun, which = "distribution", exact = TRUE)
  if (is.null(optim_args)) optim_args <- attr(objFun, which = "optim_args", exact = TRUE)

  stopifnot( is.numeric(verbose), length(verbose) == 1L )


  optObj <- NULL

  try(expr = {optObj <- purrr::exec(stats::optim, !!! optim_args)},
      silent = TRUE)

  if (is.null(optObj)){
    if (verbose >= 1) warning("MSE-optimization failed during delay fit!")
  } else if ( isTRUE(optObj$convergence > 0L) ){
    # do a 2nd attempt of optim in case it did not converge in the 1st place
    if (verbose >= 2) message("No proper convergence during 1st optimization in delay fit. Re-try with different parameter scaling.")

    stopifnot( "control"  %in% names(optim_args),
               "parscale" %in% names(optim_args$control) )

    # Use parameter values of non-converged fit as new start values (and adapt parscale accordingly)?
    #+We assume that for objFun smaller is better (minimization problem)
    if ( isTRUE(is.numeric(optObj$par) && all(is.finite(optObj$par)) && optObj$value < objFun(optim_args$par)) ){
      newpar <- optObj$par
      newparscale <- pmin(pmax(newpar, 1e-8), 1e+8)

      optim_args <- optim_args %>%
        purrr::assign_in(where = "par", value = newpar) %>%
        purrr::assign_in(where = c("control", "parscale"), value = newparscale)


      try(expr = {optObj <- purrr::exec(stats::optim, !!! optim_args)},
          silent = TRUE)

      #XXX should we update the optim_args: start values when we used the default optim_args?
      if ( isTRUE(optObj$convergence > 0L && verbose >= 1) ) warning("No proper convergence after re-try.")
    }## fi rescaling for 2nd attempt
  }## fi 2nd attempt necessary?

  optObj
}

#' Fit a delayed Exponential or Weibull model to one or two given sample(s).
#'
#' Maximum product spacing is used to fit the parameters.
#' Numerical optimization is done by `stats::optim`.
#' @param x numeric. observations of 1st group. Can also be a list of data from two groups.
#' @param y numeric. observations from 2nd group
#' @param bind character. parameter names that are bind together in 2-group situation.
#' @param verbose integer. level of verboseness. Default 0 is quiet.
#' @return `mps_fit` object that contains the information of the delayed model fit. Or `NULL` if optimization failed (e.g. too few observations).
#' @export
delay_model <- function(x, y = NULL, distribution = c("exponential", "weibull"), bind=NULL, verbose = 0) {
  if (is.list(x) && length(x) == 2L){
    y <- x[[2L]]
    x <- x[[1L]]
  }

  twoGr <- ! is.null(y)
  distribution <- match.arg(distribution)
  objFun <- geomSpaceFactory(x = x, y = y, distribution = distribution, bind = bind)

  optObj <- delay_fit(objFun, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  parV <- purrr::set_names(optObj$par, nm = getDist(distribution, type = "param",
                                                    twoGroup = twoGr, bind = bind))

  structure(
    list(
      data = if (twoGr) list(x = x, y = y) else x,
      distribution = distribution,
      twoGroup = twoGr,
      bind = bind,
      objFun = objFun,
      par = parV,
      val = objFun(parV),
      convergence = optObj$convergence
    ), class = "mps_fit")
}

#' @export
print.mps_fit <- function(x){
  cat( glue::glue_data(x, .sep = "\n",
                 "Fit a delayed {distribution} via Maximum Product Spacing for {if (twoGroup) 'two independent groups' else 'a single group'}.",
                 "Data: {if (twoGroup) paste(lengths(data), collapse = ' and ') else length(data)} observations, ranging from {paste(signif(range(data), 4), collapse = ' to ')}",
                 "Fitted coefficients: {coef(x) %>% signif(5) %>% paste(paste('\n  ', names(.)), ., sep = ': ', collapse = ' ')}\n\n")
  )
}

#' Coefficients of a delay-model fit.
#' @param group character string to request the canonical parameter for one group
#' @export
coef.mps_fit <- function(object, group = NULL){

  #stopifnot( inherits(object, "mps_fit") )
  par <- object[["par"]]

  if (! object[["twoGroup"]] || is.null(group)) par else {

    # original parameter names of distribution
    oNames <- getDist(object[["distribution"]], type = "param", twoGroup = FALSE, bind = NULL)
    # contract: bind was intersected with parameter names and, hence, has right order
    bind <- object[["bind"]]
    # extract all group parameters and restore original name (e.g. remove ".x")
    par.gr <- purrr::set_names(par[grepl(pattern = paste0(".", group), x = names(par), fixed = TRUE)],
                               nm = setdiff(oNames, bind))

    # restore original order
    c(par.gr, par[bind])[oNames]
  }
}

#' @export
summary.mps_fit <- function(object, ...){
  print(object)
}


#' @export
plot.mps_fit <- function(x, y, title, subtitle, ...){
  stopifnot( inherits(x, "mps_fit") )
  # parameter y comes from the plot-generic. y is not used here.

  cumFun <- getDist(x[["distribution"]], type = "cdf")

  p <- grNames <- NULL

  # catch the one-group case!
  if ( isTRUE(x[["twoGroup"]]) ){
    stopifnot( is.list(x[["data"]]) )

    grNames <- names(x[["data"]])

    p <- ggplot2::ggplot(data = x[["data"]] %>% unlist %>%
                           tibble::enframe(name = "group") %>%
                           dplyr::mutate(group = substr(x = group, 1L, 1L)),
                         mapping = ggplot2::aes(x = value, col = group)) +
      # add estimated delay model(s)
      purrr::map(grNames,
                 .f = ~ ggplot2::geom_function(mapping = ggplot2::aes(col = .x), inherit.aes = FALSE,
                                               fun = cumFun,
                                               args = coef(x, group = .x), linetype = "dashed"))
  } else {
    grNames <- "x"
    p <- ggplot2::ggplot(data = tibble::tibble(value=x[["data"]]),
                         mapping = ggplot2::aes(x = value)) +
      # add estimated delay model
      ggplot2::geom_function(inherit.aes = FALSE,
                             fun = cumFun,
                             args = coef(x, group = grNames), linetype = "dashed")
  }


  if (missing(title)) title <- glue::glue_data(x, "Fitted {distribution} delay {c('model', 'models')[1+twoGroup]}")
  if (missing(subtitle)) subtitle <- paste(purrr::map_chr(grNames,
                                                          ~ paste(names(coef(x, group = .)), signif(coef(x, group = .), 4),
                                                                  sep = ": ", collapse = ", ")),
                                           collapse = " - ")


  p +
    ggplot2::stat_ecdf(pad=TRUE) +
    ggplot2::xlim(0L, NA) +
    ggplot2::labs(x = "Time", y = "Cumulative prop. of events") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::scale_y_reverse()
}



# delay test ----

#' Test the difference for delay model parameter(s) between two uncorrelated groups.
#'
#' It is in fact a model comparison between a null model where the parameters are enforced to be equal and an unconstrained full model.
#' As test statistic we use the difference in best (=lowest) objective function value (`val_0` - `val_1`).
#' This is reminiscent of a likelihood ratio test statistic albeit the objective function is not a negative log-likelihood
#' but the negative of the maximum product spacing metric.
#'
#' High values of this difference speak against the null-model (i.e. high `val_0` indicates bad fit under 0-model and low values of `val_1` indicate a good fit under the more general model1.
#' The test is implemented as a parametric bootstrap test, i.e.
#'
#' 1. we take given null-model fit as ground truth
#' 2. we regenerate data according to this model.
#' 3. we recalculate the test statistic
#' 4. we appraise the observed test statistic in light of the generated distribution under H0
#'
#'
#' @param x data from reference/control group.
#' @param y data from the treatment group.
#' @param distribution character[1]. Name of the parametric delay distribution to use.
#' @param param character[1]. Parameter to test difference for. Default value is `'delay'`.
#' @param R numeric[1]. Number of bootstrap samples to evaluate the distribution of the test statistic.
#' @return test result: bootstrap P-value
#' @export
test_delay_diff <- function(x, y, distribution = c("exponential", "weibull"), param = "delay", R = 400) {
  distribution <- match.arg(distribution)
  par_names <- getDist(distribution = distribution, type = "param")
  stopifnot( is.numeric(x), length(x) > length(par_names), is.numeric(y), length(y) > length(par_names) )
  stopifnot( is.numeric(R), length(R) == 1L, R >= 1L )
  stopifnot( is.character(param), length(param) >= 1L )

  # Test statistic calculated from the given data and the model specification.
  #
  # The test statistic takes non-negative values.
  # Higher values of the test statistic speak in favour of H1:
  # @return list containing value of test statistic and null model fir
  testStat <- function(x, y) {
    #fit0 <- fit1 <- NULL
    fit0 <- delay_model(x = x, y = y, distribution = distribution, bind = param) #}, silent = TRUE)
    fit1 <- delay_model(x = x, y = y, distribution = distribution) #}, silent = TRUE)

    # XXX do I need to check convergence of re-fits?
    # keep also fits with error-code 52: in my tests all those fits *looked* actually OK..

    if (is.null(fit0) || is.null(fit1)) invisible(NULL) else
        # higher values of T speak in favour of H1:
        #   1. fit0 has high value (=bad fit)
        #   2. fit1 has low value (=good fit)
        list(val = 2L * (fit0[["val"]] - fit1[["val"]]),
             fit0 = fit0)
  }

  # observed test statistic
  ts_obs <- testStat(x, y)
  if (is.null(ts_obs)){
    stop("Delay model failed for restricted null-model or free full model")
  }
  fit0 <- ts_obs[["fit0"]]

  # P-values based GOF-tests
  #+ H0: simpler model is sufficient and the GOF-test solely builds on fit0:
  #+ take fitted parameters for both groups under null-model
  #+ and apply cumulative distribution functions on the observed data for both groups
  cFun <- getDist(distribution = distribution, type = "cdf")
  # sorted transformed observations. sorting is needed for Anderson-Darling (AD)-test
  transf_obs <- sort(c(rlang::exec(cFun, !!! c(list(q = x), coef(fit0, group = "x"))),
                       rlang::exec(cFun, !!! c(list(q = y), coef(fit0, group = "y"))) ))

  # number of observations
  N <- length(transf_obs)
  nx <- length(x); ny <- length(y)
  stopifnot( N == nx + ny )

  # Pearson GOF-test based on Chi-square distribution.
  # under H0, expect counts according to uniform distribution
  # nbr of classes as recommended by David S. Moore (Tests of Chi-squared Type, 1986)
  gof_nClasses <- max(length(coef(fit0)) + 2L, ceiling(2L * N**.4))
  transf_tab <- tabulate(findInterval(transf_obs, vec = seq(0L, 1L, length.out = gof_nClasses+1L),
                                      rightmost.closed = TRUE, all.inside = TRUE), nbins = gof_nClasses)
  # use adjusted degrees of freedom (loose one df for each parameter estimated)
  P_gof_pearson <- pchisq(q = sum((transf_tab - mean(transf_tab))**2L) / mean(transf_tab),
                          df = gof_nClasses - length(coef(fit0)) - 1L, lower.tail = FALSE)

  # EDF-based GOF-test
  # Anderson-Darling (AD) test statistic, cf Stephens, Tests based on EDF Statistics p.101, (4.2)
  i <- seq_along(transf_obs)
  A2 <- -N - mean((2L * i - 1L) * log(transf_obs) + (2L*N + 1L - 2L*i)*log(1-transf_obs))


  P_gof_ad <- if (identical(distribution, 'exponential')){
    # modification for Exponential (cf Stephens, Table 4.14, p.138)
    # the correction factor approaches 1 from above.
    # We keep the number of N as all observations, independent of the number of parameters estimated in the null-model.
    #+XXX should we increase N by the number of parameters p estimated less 2 ( p -2 because 2 parameters are estimated in standard delayed exponential)
    A2_mod <- A2 * max(1L, 1L + 5.4 / N - 11 / N**2L)

    # P-value for A2_mod based on upper tail percentage points (table 4.14, p.138)
    # use linear extrapolation on log-scaled sig_level to get P-values
    # A2_mod_crit <- tribble(~sig_level, ~crit_val,
    #   .25, .736,
    #   .15, .916,
    #   .10, 1.062,
    #   .05, 1.321,
    #   .025, 1.591,
    #   .01, 1.959)
    # ggplot(A2_mod_crit, mapping = aes(x = crit_val, y=sig_level)) +
    #   scale_y_log10() +
    #   geom_point() +
    #   geom_line(size = .25)
    # model with log-siglevel as response and critical value as predictor
    # coef(lm(log(sig_level) ~ crit_val, data = A2_mod_crit))
    # min(1L, exp(.51014 - 2.62843 * A2_mod))

    # linear model with logits as response and sqrt(crit) as predictor
    # coef(lm(log(sig_level / (1-sig_level)) ~ sqrt(crit_val), data = A2_mod_crit))
    # stats::plogis(4.42397965 - 6.42703744 * sqrt(A2_mod))

    # P-value comes from back-transformed linear model with logits as response and crit + sqrt(crit) as predictors
    # coef(lm(log(sig_level / (1-sig_level)) ~ crit_val + sqrt(crit_val), data = A2_mod_crit))
    stats::plogis(3.8829139 - 0.4357149 * A2_mod -5.4427475 * sqrt(A2_mod))

  } else {
    stopifnot( identical(distribution, "weibull") )
    # P-value for Weibull based on Lockhart, 1994 (Table 1)
    # interpolation model on logits using critical value and inverse of shape parameter
    # fm_c3 <- lm(log(sig_level/(1-sig_level)) ~ (crit_val + sqrt(crit_val)) * c, data = crit_weib)
    # coef(fm_c3)
    # the P-value for the AD-test depends on the shape parameter estimate
    #+cf the MC-study of Lockhart for finite samples but there it assumes a single shape parameter estimate,
    #+XXX we might have two! for the time being we use the harmonic mean (the smallest mean)
    #+which might overstate shape_inverse and is hence conservative for the P-value
    shape_pars <- coef(fit0)[grepl("shape", names(coef(fit0)))]
    stopifnot( length(shape_pars) <= 2L )
    # harmonic mean of shape parameters
    shape_hm <- if (length(shape_pars) == 1L) shape_pars[[1L]] else 2L * prod(shape_pars) / sum(shape_pars)
    shape_inv <- min(0.5, shape_inv) # maximal value of 0.5, see Lockhart, 1994
    stats::plogis(5.358389287 -2.589581570 * A2 -8.640666568 * sqrt(A2) + 0.635487306 * shape_inv +
                    # interaction effects
                    3.626348262  * A2 * shape_inv -1.422402077 * sqrt(A2) * shape_inv)
  }


  ## QQQ add GOF Moran's test: Cheng & Stephens (1989)


  # parametric bootstrap:
  # generate R samples (x, y) by random sampling from the fitted H0-model (e.g. common delay),
  #+where all nuisance parameters are at there best fit
  # calculate the test statistic on the simulated data
  # estimate P as proportion of simulated test statistics that exceed the observed test statistic t_obs

  ranFun <- getDist(distribution, type = "r")
  # arguments to the random function generation
  ranFunArgsX <- c(list(n=length(x)), coef(fit0, group = "x"))
  ranFunArgsY <- c(list(n=length(y)), coef(fit0, group = "y"))

  t0_dist <- future.apply::future_vapply(X = seq(R), FUN.VALUE = double(1L),
                                         FUN = function(dummy){

    # generate new data according to given fitted null-model
    # sort is not needed here, as it goes through the whole pipeline (factory method)
    #
    # ts_boot <- NA_real_
    # try(expr = {ts_boot <- testStat(x = rlang::exec(ranFun, !!! ranFunArgsX),
    #                                 y = rlang::exec(ranFun, !!! ranFunArgsY))[["val"]]},
    #     silent = TRUE)
    ts_boot <- testStat(x = rlang::exec(ranFun, !!! ranFunArgsX),
                        y = rlang::exec(ranFun, !!! ranFunArgsY))[["val"]]
    if (is.null(ts_boot)) ts_boot <- NA_real_

    ts_boot

  }, future.seed = TRUE)

  t0_dist <- t0_dist[is.finite(t0_dist)]
  chisq_df_hat <- coef(MASS::fitdistr(x = t0_dist, densfun = "chi-squared",
                                start = list(df = length(param)),
                                method = "Brent", lower = .001, upper = 401))

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
                  gof_pearson = P_gof_pearson,
                  gof_ad = P_gof_ad,
                  lr = P_lr,
                  lr_pp = P_lr_pp)
    ), class = "test_delay")
}


#' @export
plot.test_delay <- function(x, y, title, subtitle, ...){
  stopifnot(inherits(x, "test_delay"))

  teststat <- x[["testDist"]]

  if (missing(title)) title <- glue::glue("Distribution of test statistic under H0 for parameter {dQuote(x$param)}")
  if (missing(subtitle)) subtitle <- glue::glue("Sampling distribution, based on {length(teststat)} parametric bootstrap drawings. ",
                                                "Approximated by a chi-square distribution with df={signif(x[['chisq_df_hat']], 2)}.")


  p <- dplyr::tibble(teststat = teststat) %>%
    ggplot2::ggplot(ggplot2::aes(x = teststat, y = ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(bins = 11L + ceiling(sqrt(length(teststat))))

  # extract maximum density value
  ymax <- ggplot2::layer_data(p) %>%
    dplyr::pull(y) %>%
    {ceiling(max(. + .1, . * 1.01))}

  p + ggplot2::geom_vline(xintercept = x[["t_obs"]], linetype = "dashed", colour = "grey") +
    ggplot2::geom_function(inherit.aes = FALSE,
                           fun = stats::dchisq, args = list(df = x[['chisq_df_hat']]),
                           col = "red", linetype = "dotted") +
    ggplot2::coord_cartesian(ylim = c(0L, ymax)) +
    ggplot2::labs(x = "Test statistic",
                  title = title, subtitle = subtitle)

}


#' Power simulation function for a two-group comparison of the delay parameter.
#'
#' Given an effect size and a sample size `n` it simulates the power.
#' The higher number of power simulation rounds the more densely the space of data according to the specified model is sampled.
#'
#' @param eff list. With model parameters for each of the two groups. Control group is given first.
#' @param param character. Parameter name for which to simulate the power.
#' @param n integer. Number of observations per group for the power simulation. Can be two different numbers, control group and then treatment group.
#' @param nPowerSim integer. Number of simulation rounds. Default value 1600 yields a standard error of 0.01 for power if the true power is 80 percent.
#' @param R integer. Number of bootstrap samples to assess difference in parameter within each power simulation.
#' @return `sscn`-object with power estimate from simulations. Or `NULL` in case of errors.
#' @export
ssc_delay_sim_power <- function(distribution = c("exponential", "weibull"), param = "delay",
                             eff = stop("Provide parameters for each group that reflect the effect!"),
                             n, sig.level = 0.05, nPowerSim = 16e2, R = 4e2){

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
    message("Too few observations to fit parameters.")
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
    try(P_val <- test_delay_diff(x = dat_ctrl, y =dat_trtm,
                    distribution = distribution, param = param, R = R)[["P"]],
        silent = TRUE)
    P_val
  }, future.seed = TRUE)


  P_dist <- P_dist[is.finite(P_dist)]

  if ( !length(P_dist) )
    warning("No valid power simulation results.") else if ( length(P_dist) < 100L )
      warning("Low resultion for power estimate.")

  structure(
    list(id = "delay:2groups", name = "Delay: Difference in delay for time-to-event data in two groups",
         eff = eff, sig.level = sig.level, n = n, N = sum(n),
         P_dist = P_dist, ##debug
         ##dat_mean_c = m_ctrl, ##debug
         power = if (length(P_dist)) sum(P_dist < sig.level) / length(P_dist) else NA_real_
    ),
    class = "sscn"
  )
}

