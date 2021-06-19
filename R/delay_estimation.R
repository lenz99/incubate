
#' Extract the parameters for the specified group.
#' This is an internal helper function
#' used in `coef.mps_fit` and in the factory method `geomSpaceFactory` below
#' @param par named parameters (as simple vector or as list both works)
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
#' @param bind character(1). parameter names that are bind together (i.e. equated) between both groups
#' @return objective function
geomSpaceFactory <- function(x, y=NULL, distribution = c("exponential", "weibull"), bind=NULL) {
  twoGr <- ! is.null(y)
  stopifnot( is.numeric(x), length(x) > 0L, ! twoGr || is.numeric(y) && length(y) > 0L )
  stopifnot( is.null(bind) || is.character(bind) && length(bind) >= 1L )

  # intersect enforces the canonical order of dist-parameters in bind!
  oNames <- getDist(distribution, type = "param", twoGroup = FALSE, bind = NULL) #standard ('original') names
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

    if ( !length(y) ) {
      warning("No valid data in y! Only non-negative and finite real values are valid.")
      return(invisible(NULL))
    }
  } #fi twoGr

  distribution <- match.arg(distribution)

  par_names <- getDist(distribution = distribution, type = "param",
                       twoGroup = twoGr, bind = bind)

  optim_args <- NULL



  # get parameter setting for a group of observations
  # @return list with par, lower and upper information per parameter
  getParSetting.gr <- function(obs){
    # contract: obs is sorted!

    LOW_MIN <- 1e-9

    switch (EXPR = distribution,
            exponential = {

              list(
                par = c(max(LOW_MIN, min(obs)-2L/length(obs)), 1L/mean(obs - min(obs) + 2L/length(obs))),
                lower = purrr::set_names(c(0L, LOW_MIN), oNames),
                upper = purrr::set_names(c(max(LOW_MIN, min(obs)-1L/length(obs), min(obs)*.9999), +Inf), oNames)
              )

            },
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

              list(
                par = c(max(LOW_MIN, min(obs) - 2L/length(obs)), start_shape, start_scale),
                lower = purrr::set_names(c(0L, LOW_MIN, LOW_MIN), oNames),
                upper = purrr::set_names(c(max(LOW_MIN, min(obs)-1L/length(obs), min(obs)*.9999), +Inf, +Inf), oNames)
              )},
            stop("Unknown distribution provided!")
    )
  }


  par0_x <- getParSetting.gr(x)

  optim_args <-
    if (! twoGr) {
      append(par0_x, list(method = "L-BFGS-B",
                          # QQQ something like exp(trunc(log(par))) where par is the start parameters
                          control = list(parscale = pmin.int(1e7, pmax.int(1e-7, par0_x[["par"]]))) )
      )

    } else {
      stopifnot( twoGr )

      # all parameters are bound: treat x and y as a single group
      if ( length(bind) == length(oNames) ) {
        par0_xy <- getParSetting.gr(c(x,y))
        append(par0_xy,
               list(method = "L-BFGS-B",
                    control = list(parscale = pmin.int(1e7, pmax.int(1e-7, par0_xy[["par"]]))) )
        )

      } else { #twoGr, not all params bound!

        par0_y <- getParSetting.gr(y)

        start2 <- lower2 <- upper2 <- NULL

        start_x <- par0_x[["par"]]
        start_y <- par0_y[["par"]]
        lower_x <- par0_x[["lower"]]
        lower_y <- par0_y[["lower"]]
        upper_x <- par0_x[["upper"]]
        upper_y <- par0_y[["upper"]]


        #XXX check here if we need min or max for lower/upper?! Think.
        #XXX actually, only delay needs special attention, rest is simple! Can we simplify the code?!
        if (distribution == 'exponential') {

          start2 <- if (is.null(bind)) c(start_x, start_y) else
            if (identical(bind, 'delay')) c(min(start_x[1L], start_y[1L]), start_x[-1L], start_y[-1L]) else
              if (identical(bind, 'rate')) c((start_x[2L] + start_y[2L])/2L, start_x[-2L], start_y[-2L])

          lower2 <- if (is.null(bind)) c(lower_x, lower_y) else
            if (identical(bind, 'delay')) c(max(lower_x[1L], lower_y[1L]), lower_x[-1L], lower_y[-1L]) else
              if (identical(bind, 'rate')) c(min(lower_x[2L], lower_y[2L]), lower_x[-2L], lower_y[-2L])

          upper2 <- if (is.null(bind)) c(upper_x, upper_y) else
            if (identical(bind, 'delay')) c(min(upper_x[1L], upper_y[1L]), upper_x[-1L], upper_y[-1L]) else
              if (identical(bind, 'rate')) c(min(upper_x[2L], upper_y[2L]), upper_x[-2L], upper_y[-2L])


        } else {
          stopifnot( distribution == 'weibull' )

          start2 <- if (is.null(bind)) c(start_x, start_y) else
            if (identical(bind, 'delay')) c(min(start_x[1L], start_y[1L]), start_x[-1L], start_y[-1L]) else
              if (identical(bind, 'shape')) c('XXX') else
                stop("Currently, for Weibull only bind=NULL and bind='delay' are supported!")

          lower2 <- if (is.null(bind)) c(lower_x, lower_y) else
            if (identical(bind, 'delay')) c(max(lower_x[1L], lower_y[1L]), lower_x[-1L], lower_y[-1L]) else
              stop("Currently, for Weibull only bind=NULL and bind='delay' are supported!")

          upper2 <- if (is.null(bind)) c(upper_x, upper_y) else
            if (identical(bind, 'delay')) c(min(upper_x[1L], upper_y[1L]), upper_x[-1L], upper_y[-1L]) else
              stop("Currently, for Weibull only bind=NULL and bind='delay' are supported!")

        } # weibull

        list(par = start2,
             method = "L-BFGS-B",
             # QQQ something like exp(trunc(log(par))) where par is the start parameters
             control = list(parscale = pmin.int(1e7, pmax.int(1e-7, start2))),
             lower = purrr::set_names(lower2, par_names),
             upper = purrr::set_names(upper2, par_names)
        )
      }
    } # twoGrp



  # calculate the differences in EDF (for given parameters in group) of adjacent observations on log scale
  getCumDiffs <- function(data, pars, group){
    pars.gr <- getPars(pars, group = group, twoGr = twoGr, oNames = oNames, bind = bind)

    # contract: x is sorted!
    cumDiffs <- diff(c(0L,
                       purrr::exec(getDist(distribution, type = "cdf"),
                                   !!! c(list(q=data), pars.gr)),
                       1L))

    # use densFun for ties
    ind_z <- which(cumDiffs == 0L)
    if ( length(ind_z) ){
      #ind_z[which(ind_z == 1L)] <- 2L #at least 2 to avoid idx 0 later when using x[ind_zx - 1]
      ind_z[which(ind_z > length(data))] <- length(data) # cap at length of data, then we use ind_z to directly address data
      cumDiffs[ind_z] <- purrr::exec(getDist(distribution, type = "dens"),
                                        !!! c(list(x = x[ind_z]), pars.gr))
    } #fi

    # respect the machine's numerical lower limit
    cumDiffs[which(cumDiffs < .Machine$double.xmin)] <- .Machine$double.xmin

    log(cumDiffs)
  }

  #' negative maximum spacing estimation objective function.
  #' Estimate parameters by minimizing this function.
  #' @param pars parameter vector.
  negMSE <- function(pars){
    stopifnot( length(par_names) == length(pars) )
    pars <- purrr::set_names(pars, par_names)

    #QQQ for twoGr:
    #+does it make a difference to first merge x and y and then do the cumDiffs, log and mean
    - if (! twoGr) mean(getCumDiffs(x, pars, group = "x")) else
      weighted.mean(c(mean(getCumDiffs(x, pars, group = "x")), mean(getCumDiffs(y, pars, group = "y"))),
                    w = c(length(x), length(y)))

  }

  # add "optim_args" & distribution as attributes to the objective function
  attr(negMSE, which = "optim_args") <- c(list(fn = negMSE), optim_args) #optim_args
  attr(negMSE, which = "distribution") <- distribution
  attr(negMSE, which = "twoGroup") <- twoGr
  attr(negMSE, which = "bind") <- bind

  negMSE
}


#' Parameter fitting according to MSE through numerical optimization.
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
#' @return named coefficient vector
#' @export
coef.mps_fit <- function(object, group = NULL){

  #stopifnot( inherits(object, "mps_fit") )

  # original parameter names of distribution
  oNames <- getDist(object[["distribution"]], type = "param", twoGroup = FALSE, bind = NULL)
  # contract: bind was intersected with parameter names and, hence, has right order
  bind <- object[["bind"]]

  getPars(object[["par"]], group = group, twoGr = object[["twoGroup"]], oNames = oNames, bind = bind)
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

