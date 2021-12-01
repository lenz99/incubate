
#' Extract the parameters for the specified group.
#' This is an internal helper function
#' used in `coef.mps_fit` and in the factory method `geomSpaceFactory` below
#' @param par named parameters (as simple vector or as list both work)
#' @return parameter vector for the relevant group
getPars <- function(par, group = "x", twoGr, oNames, bind) {
  if ( ! twoGr || is.null(group) ) return(par)

  stopifnot( is.character(group), nzchar(group) )

  # extract all group parameters and restore original name (e.g. remove ".x")
  par.gr <- purrr::set_names(par[grepl(pattern = paste0(".", group), x = names(par), fixed = TRUE)],
                              nm = setdiff(oNames, bind))

  # restore original order
  # contract: bind is intersected and has right order
  # contract: bind comes first in par
  c(par.gr, par[bind])[oNames]
}

#' Factory method for negative maximum spacing estimation objective function.
#'
#' Negative or infinite values are discarded.
#' @param x numeric. observations
#' @param y numeric. observations in second group.
#' @param distribution character(1). delayed distribution family
#' @param bind character. parameter names that are bind together (i.e. equated) between both groups
#' @param ties character. How to handle ties within data of a group.
#' @return objective function
geomSpaceFactory <- function(x, y=NULL, distribution = c("exponential", "weibull"), bind=NULL,
                             ties=c('equidist', 'density', 'random'), verbose = 0L) {

  # setup ----
  twoGr <- ! is.null(y)
  stopifnot( is.numeric(x), length(x) > 0L, ! twoGr || is.numeric(y) && length(y) > 0L )
  stopifnot( is.null(bind) || is.character(bind) && length(bind) >= 1L )
  ties <- match.arg(ties)
  # tie-breaking via density is difficult for two-group situation
  stopifnot( ! (ties == 'density' && twoGr) )

  #standard ('original') names of distribution
  oNames <- getDist(distribution, type = "param", twoGroup = FALSE, bind = NULL)
  #bind: intersect with oNames enforces the canonical order of dist-parameters!
  bind <- intersect(oNames, bind)

  distribution <- match.arg(distribution)


  # data preparation ----


  # Break ties in case of ties='break'
  # @param obs: data vector
  # @return sorted, cleaned up data vector or NULL in case of trouble
  preprocess <- function(obs) {


    ind_neg <- which(obs < 0L)
    if (length(ind_neg)){
      warning("Negative values in data", deparse(substitute(obs)), "! These are dropped.", call. = FALSE)
      obs <- obs[-ind_neg]
    }# fi

    # drop NA and ±Inf & sort
    obs <- sort(obs[is.finite(obs)])

    if (!length(obs)) {
      warning("No valid data! Only non-negative and finite real values are valid.", call. = FALSE)
      return(invisible(NULL))
    }# fi

    if (is.null(obs)) return(NULL)

    # tie break
    if (ties == 'density' ) return(obs) ##|| ties == 'groupedML') # groupedML not implemented yet


    diffobs <- diff(obs)
    stopifnot( all(diffobs >= 0L) ) # i.e. sorted obs

    tiesDiffInd <- which(diffobs == 0L) # < .Machine$double.xmin

    if (length(tiesDiffInd)){
      #rl <- rle(diff(tiesDiffInd))
      if (verbose > 0L){
        cat(length(tiesDiffInd) + sum(diff(tiesDiffInd)>1L) + 1L, 'tied observations in',
            sum(diff(tiesDiffInd)>1L) + 1L, #length(which(rl$values > 1L))+1L,
            'group(s) within data vector.\n')
      }

      roundOffPrecision <- estimRoundingError(obs, maxObs = 1000L)
      if (verbose > 0L){
        cat("Round-off error has magnitude", roundOffPrecision, "\n")
      }

      # rounding radius can't be wider than smallest observed diff.
      # plogis to mitigate the effect of sample size: the larger the sample the more we can «trust» the observed minimal diff
      # obs[1L] = min(obs) = diff of minimal obs with 0
      rr <- .5 * min(stats::plogis(q = length(obs), scale = 11) * diffobs[which(diffobs > 0L)],
                     # rounding precision here
                     roundOffPrecision, obs[1L], na.rm = TRUE)

      ## modify tied observations per group of ties
      startInd <- endInd <- 1L
      repeat {
        #proceed to end of tie-group
        while (endInd < length(tiesDiffInd) && tiesDiffInd[endInd+1L] == tiesDiffInd[endInd] + 1L) {endInd <- endInd+1L}
        #include adjacent index to complete tie-group
        obsInd <- c(tiesDiffInd[startInd:endInd], tiesDiffInd[endInd]+1L)
        stopifnot( sd(obs[obsInd]) == 0L ) #check: tie-group
        obs[obsInd] <- obs[obsInd] + if (ties == 'random') {
          # sort ensures that data after tie-break is still sorted from small to large
          sort(runif(n = length(obsInd), min = -rr, max = +rr)) } else {
            # ties == 'equidist'
            # use evenly spaced observations to break tie as proposed by Cheng (1989) on Moran test statistic
            #+They first use the ties = 'density' approach for initial estimation of parameters for Moran's statistic
            seq.int(from = -rr, to = +rr, length.out = length(obsInd))
          }
        startInd <- endInd <- endInd+1L
        if ( startInd > length(tiesDiffInd) ) break
      } #repeat

      if (verbose > 1L){ cat("New data: ", paste(obs, collapse = ", "), "\n")}
    } #fi tiesdiff

    # we have broken all ties
    stopifnot( !any(diff(obs)==0L) )

    obs
  } #fn preprocess

  if (is.null({x <- preprocess(obs = x)})) return(invisible(NULL))
  if (isTRUE(twoGr) && is.null({y <- preprocess(obs = y)})) return(invisible(NULL))


  # optimization arguments -----
  par_names <- getDist(distribution = distribution, type = "param",
                       twoGroup = twoGr, bind = bind)



  # get start values for a group of observations
  # @return list with par and upper limit for delay
  getParSetting.gr <- function(obs){
    # contract: obs is sorted!

    DELAY_MIN <- 1e-9

    parV <- switch (EXPR = distribution,
                    # min(obs) = obs[1L]
                    exponential = c( max(DELAY_MIN, obs[1L] - 2/length(obs)),
                                     mean(obs - obs[1L] + 2/length(obs))**-1L ),
                    weibull = {
                      # start values from 'Weibull plot'
                      #+using the empirical distribution function
                      ## in MASS::fitdistr they simplify:
                      # lx <- log(x)
                      # m <- mean(lx)
                      # v <- var(lx)
                      # shape <- 1.2/sqrt(v)
                      # scale <- exp(m + 0.572/shape)
                      # take out extreme values for robustness (against possible outliers)
                      #+ when at least 20 observations
                      obs_f <- obs[max(1L, floor(length(obs)*.02)):ceiling(length(obs)*.95)] #assume sorted data
                      start_y <- log(-log(1-(seq_along(obs_f)-.3)/(length(obs_f)+.4)))
                      # cf. lm.fit(x = cbind(1, log(obs)), y = start_y))$coefficients
                      start_shape <- cor(log(obs_f), start_y) * sd(start_y) / sd(log(obs_f))
                      start_scale <- exp(mean(log(obs_f) - mean(start_y) / start_shape))


                      c( max(DELAY_MIN, obs[1L] - 2/length(obs)),
                         start_shape,
                         start_scale )
                    },
                    stop("Unknown distribution provided!", call. = FALSE)
    )

    list(
      par = parV,
      delay_upper = max(DELAY_MIN, min(obs) - .1/length(obs), min(obs)*.99999)
    )
  }# fn getParSetting.gr


  # parameter bounds: lower & upper
  lowerVec <- upperVec <- purrr::set_names(rep(NA_real_, length(par_names)),
                                           nm = par_names)

  PAR_LOW <- 1e-13
  PAR_BOUNDS <- list(delay = c(lower = 0, upper = NA_real_),
                     rate  = c(lower = PAR_LOW, upper = +Inf),
                     shape = c(lower = PAR_LOW, upper = +Inf),
                     scale = c(lower = PAR_LOW, upper = +Inf))


  # alas, purrr::iwalk did not work for me here
  for (na in names(PAR_BOUNDS)) {
    idx <- startsWith(par_names, prefix = na)
    if (any(idx)) {
      lowerVec[idx] <- purrr::chuck(PAR_BOUNDS, na, 'lower')
      upperVec[idx] <- purrr::chuck(PAR_BOUNDS, na, 'upper')
    } #fi
  } #rof

  par0_x <- getParSetting.gr(x)
  parV <-
    if (! twoGr) {
      upperVec['delay'] <- par0_x[['delay_upper']]
      par0_x[['par']]

    } else {
      stopifnot( twoGr )

      # all parameters are bound
      if ( length(bind) == length(oNames) ) {
        # treat x and y as a single group for start value heuristic
        par0_xy <- getParSetting.gr(c(x,y))

        upperVec['delay'] <- par0_xy[['delay_upper']]
        par0_xy[['par']]


      } else { #twoGr, not all params bound!

        par0_y <- getParSetting.gr(y)

        start_x <- par0_x[['par']]
        start_y <- par0_y[['par']]

        # set upper bound for delay parameter(s)!
        if ('delay' %in% bind) {
          upperVec['delay'] <- min(par0_x[['delay_upper']], par0_y[['delay_upper']])
        } else {
          upperVec['delay.x'] <- par0_x[['delay_upper']]
          upperVec['delay.y'] <- par0_y[['delay_upper']]
        } # fi


        if (distribution == 'exponential') {

          if (is.null(bind)) c(start_x, start_y) else
            if (identical(bind, 'delay')) c(min(start_x[[1L]], start_y[[1L]]), start_x[-1L], start_y[-1L]) else
              # harmonic mean of rates (go towards higher variability setting)
              if (identical(bind, 'rate')) c(2L * start_x[[2L]] * start_y[[2L]] / (start_x[[2L]] + start_y[[2L]]), start_x[-2L], start_y[-2L])


        } else {
          stopifnot( distribution == 'weibull' )

          switch( EXPR = paste(bind, collapse = '+'),
                  delay = {
                    c(min(start_x[[1L]], start_y[[1L]]), start_x[-1L], start_y[-1L])
                  },
                  shape = {
                    # geometric mean of shapes
                    c(sqrt(start_x[[2L]] * start_y[[2L]]), start_x[-2L], start_y[-2L])
                  },
                  scale = {
                    # arithmetic mean of scales (corresponds to harmonic mean of rates)
                    c((start_x[[3L]] + start_y[[3L]])/2L, start_x[-3L], start_y[-3L])
                  },
                  `delay+shape` = {
                    c(min(start_x[[1L]], start_y[[1L]]), sqrt(start_x[[2L]] * start_y[[2L]]), start_x[[3]], start_y[[3]])
                  },
                  `delay+scale` = {
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[3L]] + start_y[[3L]])/2L, start_x[[2]], start_y[[2]])
                  },
                  `shape+scale` = {
                    c(sqrt(start_x[[2L]] * start_y[[2L]]), (start_x[[3L]] + start_y[[3L]])/2L, start_x[[1L]], start_y[[1L]])
                  },
                  # default: bind=NULL
                  {
                    stopifnot( is.null(bind) )
                    c(start_x, start_y)
                  })

        } # weibull
      } #twoGr, not all params bound!
    } # twoGrp

  stopifnot( ! any(is.na(upperVec), is.na(lowerVec)) )

  optim_args <- list(
    par = parV,
    lower = lowerVec,
    upper = upperVec,
    method = "L-BFGS-B",
    # QQQ something like exp(trunc(log(par))) where par is the start parameters
    control = list(parscale = pmin.int(1e11, pmax.int(1e-11, parV)))
  )




  # objective function ----

  # log spacings:
  # calculate the differences in EDF (for given parameters in group) of adjacent observations on log scale
  # @return n+1 cumulative diffs on log-scale
  getCumDiffs <- function(pars, group){
    obs <- get(group)
    pars.gr <- getPars(pars, group = group, twoGr = twoGr, oNames = oNames, bind = bind)

    # calculate spacings
    # contract: data is sorted!
    cumDiffs <- diff(c(0L,
                       purrr::exec(getDist(distribution, type = "cdf"), !!! c(list(q=obs), pars.gr)),
                       1L))


    # use densFun for ties
    # we check difference of obs directly (not cumDiffs)
    #+because cumDiffs can be 0 even if obs are different, in particular for non-suitable parameters!
    ind_t <- which(diff(obs) == 0L)
    if ( length(ind_t) ){
      stopifnot( ties == 'density' ) # other tie-strategies have already dealt with ties in preprocess
      # increase index by 1 to get from diff(obs)-indices to cumDiffs-indices
      cumDiffs[1L+ind_t] <- purrr::exec(getDist(distribution, type = "dens"), !!! c(list(x = obs[ind_t]), pars.gr))
    } #fi

    # respect the machine's numerical lower limit
    cumDiffs[which(cumDiffs < .Machine$double.xmin)] <- .Machine$double.xmin

    log(cumDiffs)

  }# fn getCumDiffs


  #' negative maximum spacing estimation objective function.
  #' Estimate parameters by minimizing this function.
  #' @param pars parameter vector.
  #' @param aggregated logical. For two group case, if `FALSE` return individual mean log cum-diffs per group
  negMSE <- function(pars, aggregated = TRUE){
    stopifnot( length(par_names) == length(pars) )
    pars <- purrr::set_names(pars, par_names)

    - if (! twoGr) {
      mean(getCumDiffs(pars, group = "x"))
    } else {
      #twoGr:
      #the approach to first merge x and y and then do the cumDiffs, log and mean does *not* work out
      #because the parameters should be optimized within group.
      #merged data lead to frequent non-convergence or visually bad fits
      res <- c(mean(getCumDiffs(pars, group = "x")), mean(getCumDiffs(pars, group = "y")))
      if (aggregated)
        weighted.mean(res, w = c(length(x), length(y)))
      else
        res
    }

  }

  # add preprocessed data as attribute
  attr(negMSE, which = 'x') <- x
  if (twoGr) { attr(negMSE, which = 'y') <- y}
  # add "optim_args" & distribution as attributes to the objective function
  attr(negMSE, which = "optim_args") <- c(list(fn = negMSE), optim_args) #optim_args
  attr(negMSE, which = "distribution") <- distribution
  attr(negMSE, which = "twoGroup") <- twoGr
  # when bind is NULL there will be no attributed named 'bind'
  attr(negMSE, which = "bind") <- bind
  # method how to handle ties
  attr(negMSE, which = "ties") <- ties

  negMSE
}


#' Parameter fitting according to MSE through numerical optimization.
#'
#' The objective function carries the given data in its environment and it is to be minimized.
#' R's standard routine `stats::optim` does the optimization, using numerical derivatives.
#' @param objFun objective function to be minimized
#' @param optim_args list of own arguments for optimization. If `NULL` it uses the default optim arguments associated to the objective function.
#' @param verbose integer that indicates the level of verboseness. Default 0 is quiet.
#' @return optimization object including a named parameter vector or `NULL` in case of errors during optimization
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
    if (verbose > 0L) warning("MSE-optimization failed during delay fit!", call. = FALSE)
  } else if ( isTRUE(optObj$convergence > 0L) ){
    # do a 2nd attempt of optim in case it did not converge in the first place
    if (verbose > 1L) message("No proper convergence during 1st optimization in delay fit. Re-try with different parameter scaling.")

    stopifnot( "control"  %in% names(optim_args),
               "parscale" %in% names(optim_args$control) )

    # Use parameter values of non-converged fit as new start values (and adapt parscale accordingly)?
    #+The objFun is to be minimized,  smaller is better!
    if ( isTRUE(is.numeric(optObj$par) && all(is.finite(optObj$par)) && optObj$value < objFun(optim_args$par)) ){
      optim_args[['par']] <- optObj$par  # purrr::assign_in(where = "par", value = optObj$par)
      newparsc <- abs(optim_args[['par']])
      newparsc[which(newparsc < 1e-7)] <- 1e-7
      newparsc[which(newparsc > 1e8)] <- 1e8
      optim_args[['control']][['parscale']] <- newparsc

      optObj <- NULL
      try(expr = {optObj <- purrr::exec(stats::optim, !!! optim_args)},
          silent = TRUE)

      #XXX should we update the optim_args: start values when we used the default optim_args?
      if ( is.null(optObj) || isTRUE(optObj$convergence > 0L && verbose > 0L) ) warning("No proper convergence after re-try.", call. = FALSE)
    }## fi rescaling for 2nd attempt
  }## fi 2nd attempt necessary?

  # set names to parameter vector
  if (! is.null(optObj)){
    optObj$par <- purrr::set_names(optObj$par,
                                   nm = getDist(distribution, type = "param",
                                                twoGroup = attr(objFun, which = "twoGroup", exact = TRUE),
                                                bind = attr(objFun, which = "bind", exact = TRUE)))
  }

  optObj
}

#' Fit a delayed Exponential or Weibull model to one or two given sample(s).
#'
#' Maximum product spacing is used to fit the parameters.
#' Numerical optimization is done by `stats::optim`.
#' @param x numeric. observations of 1st group. Can also be a list of data from two groups.
#' @param y numeric. observations from 2nd group
#' @param bind character. parameter names that are bind together in 2-group situation.
#' @param ties character. How to handle ties.
#' @param optim_args list. optimization arguments to use. Use `NULL` to use the data-dependent default values.
#' @param verbose integer. level of verboseness. Default 0 is quiet.
#' @return `mps_fit` object that contains the information of the delayed model fit. Or `NULL` if optimization failed (e.g. too few observations).
#' @export
delay_model <- function(x, y = NULL, distribution = c("exponential", "weibull"), bind=NULL,
                        ties=c('equidist', 'density', 'random'), optim_args=NULL, verbose = 0) {

  # unpack x if it is a list of two vectors
  if (is.list(x)){
    stopifnot( length(x) == 2L )
    y <- x[[2L]]
    x <- x[[1L]]
  }

  twoGr <- ! is.null(y)
  distribution <- match.arg(distribution)
  ties <- match.arg(ties)
  objFun <- geomSpaceFactory(x = x, y = y, distribution = distribution, bind = bind, ties = ties, verbose = verbose)
  # set preprocessed data
  x <- attr(objFun, 'x', exact = TRUE)
  if (twoGr) { y <- attr(objFun, 'y', exact = TRUE)}

  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  # get CDF-transformed data
  oNames <- getDist(distribution, type = 'param')
  data_tr <- purrr::exec(getDist(distribution, type = "cdf"),
                         !!! c(list(q=x), getPars(par = optObj$par, group = 'x', twoGr = twoGr, oNames = oNames, bind = bind)))
  if (twoGr) data_tr <-
    list(x = data_tr,
         y = purrr::exec(getDist(distribution, type = "cdf"),
                         !!! c(list(q=y), getPars(par = optObj$par, group = 'y', twoGr = twoGr, oNames = oNames, bind = bind))))
  structure(
    list(
      data = if (twoGr) list(x = x, y = y) else x,
      data_transf = data_tr, # store CDF-transformed data
      distribution = distribution,
      twoGroup = twoGr,
      bind = bind,
      objFun = objFun,
      par = optObj$par,
      val = optObj$value, ##objFun(optObj$par),
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
#' @return named coefficient vector
#' @export
coef.mps_fit <- function(object, group = NULL){

  #stopifnot( inherits(object, "mps_fit") )

  # original parameter names of distribution
  oNames <- getDist(object[["distribution"]], type = "param", twoGroup = FALSE, bind = NULL)
  getPars(object[["par"]], group = group, twoGr = object[["twoGroup"]], oNames = oNames,
          # contract: bind was intersected with parameter names and, hence, has right order
          bind = object[["bind"]])
}

#' @export
summary.mps_fit <- function(object, ...){
  print(object)
}

#' Refit a MPS-fit with specified optimization arguments.
#' If more things need to be changed use `delay_model`.
#' @export
update.mps_fit <- function(object, optim_args, verbose = 0, ...){

  objFun <- object[['objFun']]

  ## fit model with given optim_args
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  structure(
    list(
      data = object[['data']],
      distribution = object[['distribution']],
      twoGroup = object[['twoGroup']],
      bind = object[['bind']],
      objFun = objFun,
      par = optObj$par,
      val = optObj$value, ##objFun(optObj$par),
      convergence = optObj$convergence
    ), class = "mps_fit")
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

#' Simulate Data from Fitted Model
#' @param object MPS-fit object
#' @param nsim number of simulations
#' @param seed currently unused! XXX
#' @return list of simulated data
#' @export
simulate.mps_fit <- function(object, nsim = 1, seed = NULL, ...){
  stopifnot( inherits(object, 'mps_fit'))

  ranFun <- getDist(object$distribution, type = "r")
  nObs <- if (isTRUE(object$twoGroup)) lengths(object$data) else length(object$data)

  # arguments to the random function generation
  ranFunArgsX <- as.list(c(n=nObs[[1L]], coef(object, group = "x")))
  ranFunArgsY <- if (isTRUE(object$twoGroup)) as.list(c(n=nObs[[2L]], coef(object, group = "y")))

  simExpr <- if (isTRUE(object$twoGroup))
    expression(list(x=rlang::exec(ranFun, !!! ranFunArgsX),
                    y=rlang::exec(ranFun, !!! ranFunArgsY))) else
                      expression(rlang::exec(ranFun, !!! ranFunArgsX))

  if (nsim > 1000L){
    future.apply::future_replicate(n = nsim, expr = eval(simExpr), simplify = FALSE, future.seed = TRUE)
  } else {
    replicate(n = nsim, expr = eval(simExpr), simplify = FALSE)
  }
}


#' Confidence intervals for parameters of MPS-model fits.
#'
#' Basic bootstrap confidence limits are generated. At least R=1000 bootstrap replications are recommended.
#' @param R bootstrap replications
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.mps_fit <- function(object, parm, level = 0.95, R = 99L, ...){
  stopifnot( inherits(object, 'mps_fit'))

  stopifnot( is.numeric(level), length(level) == 1L, level < 1L, level > 0L)
  cf <- coef(object)
  pnames <- names(cf)
  stopifnot( is.numeric(cf), is.character(pnames) )

  if (missing(parm)) parm <- pnames else
    if (is.numeric(parm)) parm <- pnames[parm]
  parm <- intersect(pnames, parm)

  if (is.null(parm) || ! length(parm) || any(! nzchar(parm)) ){
    warning('Invalid parameter name given in argument parm=')
    return(invisible(NULL))
  }

  stopifnot( is.character(parm), length(parm) >= 1L )

  a <- (1L - level) / 2L
  a <- c(a, 1L - a)


  # basic bootstrap confidence limits

  # get bootstrap distribution of coefficients in data from fitted model

  # for performance reasons, we 'inline' the simulate code into future_vapply, cf test_delay_diff
  ranFun <- getDist(object$distribution, type = "r")
  nObs <- if (isTRUE(object$twoGroup)) lengths(object$data) else length(object$data)

  # arguments to the random function generation
  ranFunArgsX <- as.list(c(n=nObs[[1L]], coef(object, group = "x")))
  ranFunArgsY <- if (isTRUE(object$twoGroup)) as.list(c(n=nObs[[2L]], coef(object, group = "y")))

  coefMat <- future.apply::future_vapply(X = seq.int(R), FUN.VALUE = double(length(cf)),
                                         FUN = function(dummy){
                                           x <- rlang::exec(ranFun, !!! ranFunArgsX)
                                           y <- if (isTRUE(object$twoGroup)) rlang::exec(ranFun, !!! ranFunArgsY)

                                           coef(delay_model(x=x, y=y, distribution = object$distribution, bind = object$bind))
                                         }, future.seed = TRUE)

  # more clear and shorter but less efficient!
  # coefMat <- object %>%
  #   simulate(nsim = R) %>%
  #   future.apply::future_vapply(FUN.VALUE = double(length(cf)),
  #                               FUN = function(da){
  #                                 coef(delay_model(x=da, distribution = object$distribution, bind = object$bind))
  #                               })

  # cf Davison, p28
  ci <- 2L * cf - t(apply(coefMat, 1L, quantile, probs = rev(a), na.rm = TRUE))
  colnames(ci) <- rev(colnames(ci))

  # enforce parameter bounds also for CI
  # all parameters are non-negative!
  ci[which(ci<0L)] <- 0L


  ci[parm, , drop = FALSE]
}
