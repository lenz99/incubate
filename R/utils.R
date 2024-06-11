# mkuhn, 2021-10-11
# Utility functions used within this package


TOL_NUM <- sqrt(.Machine$double.eps)
DELAY_MIN <- .Machine$double.eps #.Machine$double.xmin even smaller ##1e-9

#' Format a number as percentage.
#'
#' Internal helper function that is not exported.
#' @param x numeric vector to be formatted as percentage
#' @param digits requested number of decimal digits of the percentage
#' @return number formatted as percentage character
as_percent <- function(x, digits = 1) {
  stopifnot( is.numeric(digits) )
  sprintf(fmt = paste0('%.', as.integer(digits),'f%%', x*100))
}

#' Checks if arguments are numerically close.
#'
#' The function is vectorized and R's recycling rules apply.
#' @param x numeric first vector
#' @param y numeric second vector
#' @return logical vector if arguments from x and y are close
#' @seealso [dplyr::near()]
near <- function(x, y) {
  abs(x-y) < TOL_NUM
}

#' Minimize an objective function with alternative optimizer
#'
#' The primary optimization routine is BFGS from `stats::optim`.
#' If this fails for some reason we try an alternative which is implemented here.
#' It can use the derivative-free minimizaiton through `bobyqa` or the PORT-routine `nlminb`.
#'
#' This is only a thin wrapper to the chosen alternative optimizer.
#' @param objFun function to minimize
#' @param start vector of start values for parameters
#' @param lower numeric. lower bound for parameters (boxed constraint)
#' @param upper numeric. upper bound for parameters (boxed constraint)
#' @param verbose numeric. Verbosity level
#' @param method Specifies which optimizer to use
#' @return optimization object with some common entries like `parOpt`, `valOpt` `convergence`, `methodOpt` and `counts`. Or `NULL` in case of failure.
minObjFunAlt <- function(objFun, start, lower = -Inf, upper = +Inf, verbose = 0,
                          method = c("bobyqa", "nlminb")) {
  optObj <- NULL

  try({
    switch(method,
           nlminb = {
             optObj <- stats::nlminb(start = start, objective = objFun,
                                     lower = lower, upper = upper,
                                     control = list(trace = verbose))
             optObj$counts <- optObj$evaluations
             optObj$methodOpt <- "PORT (nlminb)"
           },
           bobyqa = {
             rhob <- min(.95, min(abs(upper-lower))/2 * .999, .2 * max(abs(start)), na.rm = TRUE)
             optObj <- minqa::bobyqa(par = start, fn = objFun,
                                     lower = lower, upper = upper,
                                     control = list(iprint = verbose,
                                                    # trust region setting: see ?bobyqa
                                                    rhobeg = rhob,
                                                    rhoend = rhob / 1e6))
             optObj$value <- optObj$fval
             optObj$message <- optObj$msg
             optObj$counts <- optObj$feval
             optObj$methodOpt <- "minqa::bobyqa"
             optObj$convergence <- optObj$ierr
           },
           stop("This optimizer-method is not supported!", call. = FALSE)
    )
  }, silent = TRUE)
  optObj
}


#' Calculate parameter scaling for optimization routine.
#'
#' The scale per parameter corresponds to the step width within the optimization path.
#' @param parV named numeric parameter vector for optimization
#' @param lowerB numeric. lower bound for parameter scales
#' @param upperB numeric. upper bound for parameter scales
#' @return numeric. vector of parameter scaling
scalePars <- function(parV, lowerB = 1e-5, upperB = 1e5){
  if (is.null(lowerB)) lowerB <- -Inf
  if (is.null(upperB)) upperB <- +Inf

  stopifnot(is.numeric(parV), is.numeric(lowerB), is.numeric(upperB))
  stopifnot(length(lowerB) == 1L || length(lowerB) == length(parV))
  stopifnot(length(upperB) == 1L || length(upperB) == length(parV))

  # scale vector: default value is 1
  scVect <- rlang::rep_along(along = parV, x = 1)

  # # non-log parameters get scaling depending on their initial value
  # idx.nonLog <- which(startsWith(names(parV), "delay1") & parV > 0)
  # scVect[idx.nonLog] <- parV[idx.nonLog]^.1 #10th root pushes towards 1
  # delay1 parameter is now also on log-scale
  idx.nonLog <- which(startsWith(names(parV), "delay1"))
  scVect[idx.nonLog] <- 1+abs(parV[idx.nonLog])/2

  # enforce upper and lower bounds
  pmax.int(lowerB, pmin.int(upperB, scVect))
}


#' Estimate rounding error based on given sample of metric values
#' The idea is to check at which level of rounding the sample values do not change.
#' @param obs numeric. Metric values from a sample to estimate the corresponding rounding error
#' @param roundDigits integer. Which level of rounding to test? Negative numbers round to corresponding powers of 10
#' @param maxObs integer. How many observations to consider at most? If the provided sample has more observations a sub-sample is used.
#' @return estimated rounding error
estimRoundingError <- function(obs, roundDigits = seq.int(from = -4L, to = 6L), n_obs = 100L) {
  stopifnot( is.numeric(obs) )

  # in case of Surv: take into account only times of event or censoring
  if (inherits(obs, "Surv")) {
    obs <- obs[, 1L]
  }

  # drop NA and Inf
  obs <- obs[is.finite(obs)]

  stopifnot( is.numeric(n_obs), length(n_obs) == 1L )
  n_obs <- trunc(n_obs)

  if (n_obs > 1L && length(obs) > n_obs) {
    obs <- obs[round(seq.int(from = 1L, to = length(obs), length.out = n_obs))]
  }

  # digits to round to
  roundDigits <- unique(trunc(roundDigits))

  #XXX think here: fails for small obs close to 0!
  if (all(abs(obs) < .1)) obs <- 1L + obs

  rDigInd <- purrr::map_lgl(.x = roundDigits,
                            .f = function(.x) all(abs(obs - round(obs, digits = .x)) < 2L * 10L**min(-2L, -.x-1L)) )

  10L**-if (!any(rDigInd)) max(roundDigits)+1L else if (all(rDigInd)) min(roundDigits)-1L else roundDigits[which.max(rDigInd)]
}


#' Check and prepare the survival response(s)
#'
#' Allowed censoring types are right-, left-, and interval-censoring.
#' If `y0` is not `NULL` this function will return either both numeric, non-Surv or both Surv-objects of the same type.
#' @param x0 response as numeric or [survival::Surv] using left, right or interval-coding
#' @param y0 response as numeric or [survival::Surv] using left, right or interval-coding
#' @param simplify logical. Should the result be as simple as possible? If `FALSE`, result will be in any case [survival::Surv] objects.
#' @return a list of the two responses, either both as [survival::Surv] or plain numeric (if no censorings and `simplify=TRUE`)
prepResponseVar <- function(x0, y0=NULL, simplify = TRUE) {

  if (missing(x0)) stop("Input for x0 is mandatory!", call. = FALSE)
  stopifnot(is.numeric(x0),  is.null(y0) || is.numeric(y0))

  isSurv.x <- is.Surv(x0)
  isSurv.y <- is.Surv(y0)

  # check if both are numeric, non-Surv (also y0=NULL)
  if (!isSurv.x && !isSurv.y) {
    # if simplify=FALSE: right-censored, all observed
    return(if (simplify) list(x=x0, y=y0) else list(x=Surv(time = x0), y=if (is.null(y0)) NULL else Surv(time = y0)))
  }

  # from here on: there are some Surv-objects

  # check if all are observed
  allobsvd.x <- !isSurv.x || all(x0[, "status"] == 1)
  allobsvd.y <- !isSurv.y || all(y0[, "status"] == 1)

  if (simplify && allobsvd.x && allobsvd.y) {
    return(list(x=if (isSurv.x) x0[,1L] else x0, y=if (isSurv.y) y0[,1L] else y0))
  }

  # from here on: return all Surv-objects (because there are censorings somewhere or simplify=FALSE)

  survType.x <- attr(x0, which = "type", exact = TRUE)
  survType.y <- attr(y0, which = "type", exact = TRUE)

  SURV_TYPES_ALLOWED <- c("right", "left", "interval")

  # when x is Surv..
  if (isSurv.x) {
    stopifnot( is.character(survType.x), nzchar(survType.x) )
    if (! survType.x %in% SURV_TYPES_ALLOWED ){
      stop("Survival-objects must be of type {",
           paste(SURV_TYPES_ALLOWED, collapse = ", "), "}!", call. = FALSE)
    }
    if (is.null(y0)) {
      # y0=NULL remains unchanged
      return(list(x=x0, y=NULL))
    } else if (isSurv.y) {
      # two Surv-objects
      if (! identical(survType.x, survType.y)) stop("Provide Surv-objects of same type!", call. = FALSE)
      return(list(x=x0, y=y0))
    } else {
      # y is numeric, non-Surv: coerce to Surv of appropriate type!
      return(list(x = x0, y = switch(survType.x,
                                     right = Surv(y0),
                                     left = Surv(y0, event = rep_len(1L, length.out = length(y0)), type = "left"),
                                     interval = Surv(y0, time2 = NA, event = rep_len(1L, length.out = length(y0)), type = "interval"),
                                     stop("This type of censoring is not supported!", call. = FALSE))))
    }
  } else {
    # x is numeric, non-Surv. y is Surv
    stopifnot(is.character(survType.y), nzchar(survType.y))
    if (!survType.y %in% SURV_TYPES_ALLOWED) {
      stop("Survival-objects must be of type {",
           paste(SURV_TYPES_ALLOWED, collapse = ", "), "}!",
           call. = FALSE)
    }

    return(list(x = switch(survType.y,
                           right = Surv(x0),
                           left = Surv(x0, event = rep_len(1L, length.out = length(x0)), type = "left"),
                           interval = Surv(x0, time2 = NA, event = rep_len(1L, length.out = length(x0)), type = "interval"),
                           stop("This type of censoring is not supported!", call. = FALSE)),
                y = y0))

  } #esle
}
