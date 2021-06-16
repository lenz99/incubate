
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

  #XXX par_start with twoGr & bind: only quick*dirty solution
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
      #+ when at least 20 observations
      xx <- obs[max(1L, floor(length(obs)*.02)):ceiling(length(obs)*.95)] #assume sorted data
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

    if (twoGr) {
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
    } # single group
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
  } #getPars

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

