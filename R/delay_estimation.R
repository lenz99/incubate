


#' Extract certain parameters from a given named parameter vector.
#'
#' This routine handles parameter vectors for one or two groups, it knows about bound parameters (cf. bind=) and
#' it can also transform the given parameters from the standard form (as used in the delayed distribution functions)
#' to the internal transformed parametrization as used by the optimization function, and vice versa.
#'
#' When parameters for a single group are requested the parameters are always returned in the canonical order of the distribution.
#' Otherwise, the order of parameters is unchanged.
#' It parses the parameter names to find out if it is twoPhase, twoGroup and if `transform=TRUE` which direction to transform.
#'
#' This is an internal helper function
#' used in [delay_model()], in the factory method [objFunFactory()] and in [coef.incubate_fit()].
#' @param parV numeric parameter vector (for single group or two groups)
#' @param distribution character name. Exponential or Weibull
#' @param group character. Extract only parameters for specified group.
#' @param transform logical. Do we want to transform the parameter vector?
#' @return requested parameters as named numeric vector
extractPars <- function(parV, distribution = c('exponential', 'weibull'), group = NULL, transform = FALSE) {

  stopifnot( is.numeric(parV), all(nzchar(names(parV))) )
  distribution <- match.arg(distribution)

  pNames <- names(parV)

  # isTransformed when all parameter names contain '_tr'
  trNames <- grepl(pattern = "_tr", pNames, fixed = TRUE)
  stopifnot( sum(trNames) %in% c(0L, length(parV)) ) #all or nothing
  isTransformed <- all(trNames)
  twoPhase <- any(grepl(pattern = "delay2", pNames, fixed = TRUE))
  # parameter names of single group where we start from!
  oNames <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, transformed = isTransformed)

  idx.nongrp <- which(!grepl(pattern = paste0("[.][xy]$"), pNames))
  #Is there any .x or .y parameter name?
  #+when two group but all parameters bound then isTwoGroup = FALSE
  isTwoGroup <- length(idx.nongrp) < length(parV)

  # transform parameters if requested
  # if no transformation required, simply extract the relevant parameters
  if ( transform ){

    # transform parameter vector for a single group
    # @param parV1: parameter vector for a single group
    # @return transformed parameter vector
    transformPars1 <- function(parV1){

      # parameter transformation matrices
      PARAM_TRANSF_M <- list(exponential = matrix(c(1, 0, 0, 0,
                                                    0, 1, 0, 0,
                                                    -1, 0, 1, 0,
                                                    0, 0, 0, 1), nrow = 4L, byrow = TRUE,
                                                  dimnames = list(c("delay1_tr", "rate1_tr", "delay2_tr", "rate2_tr"))),
                             weibull = {
                               m <- diag(nrow = 6L) #XXXX parameter transformation not factual
                               rownames(m) <- paste0(c("delay1", "shape1", "scale1", "delay2", "shape2", "scale2"), "_tr")
                               m})[[distribution]]
      PARAM_TRANSF_MINV <- list(exponential = matrix(c(1, 0, 0, 0,
                                                       0, 1, 0, 0,
                                                       1, 0, 1, 0,
                                                       0, 0, 0, 1), nrow = 4L, byrow = TRUE,
                                                     dimnames = list(c("delay1", "rate1", "delay2", "rate2"))),
                                weibull = {
                                  m <- diag(nrow = 6L)
                                  rownames(m) <- c("delay1", "shape1", "scale1", "delay2", "shape2", "scale2")
                                  m
                                })[[distribution]]


      PARAM_TRANSF_F <- list(exponential = c(identity, log, log, log),
                             weibull = c(identity, identity, identity, identity, identity, identity))[[distribution]]
      PARAM_TRANSF_FINV <- list(exponential = c(identity, exp, exp, exp),
                                weibull = c(identity, identity, identity, identity, identity, identity))[[distribution]]

      stopifnot( length(parV1) <= NCOL(PARAM_TRANSF_M), NCOL(PARAM_TRANSF_M) == length(PARAM_TRANSF_F) )


      purrr::set_names(
        if (isTransformed) {
          as.numeric(PARAM_TRANSF_MINV[seq_along(parV1), seq_along(parV1)] %*%
                       as.numeric(.mapply(FUN = function(f, x) f(x),
                                          dots = list(PARAM_TRANSF_FINV[seq_along(parV1)], parV1),
                                          MoreArgs = NULL)))
        } else {
          as.numeric(.mapply(FUN = function(f, x) f(x),
                             dots = list(PARAM_TRANSF_F[seq_along(parV1)],
                                         as.numeric(PARAM_TRANSF_M[seq_along(parV1), seq_along(parV1)] %*% parV1)),
                             MoreArgs = NULL))},

        nm = rownames(if (isTransformed) PARAM_TRANSF_MINV else PARAM_TRANSF_M)[seq_along(parV1)])
    }

    # update parameter vector according to transformation
    parV <- if (isTwoGroup) {
      # target parameter names
      pNames_trgt <- if (isTransformed) {
        sub(pattern = "_tr", replacement = "", pNames, fixed = TRUE) # remove '_tr' string
      } else {
        purrr::map_chr(strsplit(pNames, split = ".", fixed = TRUE),
                       .f = function(s)
                         if (length(s) == 1L)
                           paste0(s, "_tr") else
                             paste0(s[[1L]], "_tr.", s[[2L]]))
      }

      h <- purrr::map(.x = c("x", "y"),
                      .f = ~ {
                        idx.grp <- which(grepl(pattern = paste0("[.]", .x, "$"), pNames))
                        # drop group naming: last two letters
                        pNms_grp <- substr(pNames[idx.grp], start = 1L, stop = nchar(pNames[idx.grp], type = "chars") - 2L)
                        # apply transform on 1-group parameter vector in canonical order!
                        transformPars1(parV1 = c(parV[idx.nongrp], purrr::set_names(parV[idx.grp], nm = pNms_grp))[oNames]) })

      purrr::set_names(c(h[[1L]][pNames_trgt[idx.nongrp]], unlist(purrr::map(h, .f = ~ .x[! names(.x) %in% pNames_trgt[idx.nongrp]]))),
                       nm = pNames_trgt)
    } else transformPars1(parV1 = parV)

    # reflect transformation in meta-data
    isTransformed <- ! isTransformed
    pNames <- names(parV)
    oNames <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, transformed = isTransformed)
  }# transform


  # return parameter vector
  # single group extraction required?
  if ( ! isTwoGroup || is.null(group) || length(group) != 1L || ! group %in% c('x', 'y') ) parV else {

    stopifnot( is.character(group), nzchar(group) )

    # extract all group parameters and restore original name (e.g. remove ".x")
    parV.gr <- purrr::set_names(parV[grepl(pattern = paste0(".", group), x = pNames, fixed = TRUE)],
                                nm = setdiff(oNames, pNames[idx.nongrp]))

    # restore original order
    # contract: bind is intersected (only valid names, canonical order)
    c(parV.gr, parV[pNames[idx.nongrp]])[oNames]
  }

}


#' Calculate parameter scaling for optimization routine.
#'
#' The scale per parameter corresponds to the step width within the optimization path.
#' @param parV named numeric parameter vector
#' @param lowerB numeric. lower bound for parameter scales
#' @param upperB numeric. upper bound for parameter scales
#' @return vector of parameter scaling
scalePars <- function(parV, lowerB = 1e-5, upperB = 1e5){
  if (is.null(lowerB)) lowerB <- -Inf
  if (is.null(upperB)) upperB <- +Inf

  stopifnot( is.numeric(parV), is.numeric(lowerB), is.numeric(upperB))
  stopifnot( length(lowerB) == 1L || length(lowerB) == length(parV) )
  stopifnot( length(upperB) == 1L || length(upperB) == length(parV) )

  # scale vector: default value is 1
  scVect <- rep.int(1, times = length(parV))

  idx.nonLog <- which(startsWith(names(parV), "delay1") & parV > 0)
  scVect[idx.nonLog] <- (parV[idx.nonLog] + .01)**.3 # ca. 3th root

  # enforce upper and lower bounds
  pmax.int(lowerB, pmin.int(upperB, scVect))
}


#' Factory method for objective function, either according to maximum product of spacings estimation ('MPSE')
#' or according to standard maximum likelihood estimation ('MLE0').
#'
#' Given the observed data this factory method produces an objective function
#' which is either the negative of the MPSE-criterion H or the negative log-likelihood for MLE.
#'
#' The objective function takes a vector of model parameters as argument.
#'
#' @details
#' From the observations, negative or infinite values are discarded during pre-processing.
#' In any case, the objective function is to be minimized.
#'
#' @param x numeric. observations
#' @param y numeric. observations in second group.
#' @param method character(1). Specifies the method for which to build the objective function. Default value is `MPSE`. `MLE0` is the standard MLE-method, calculating the likelihood function as the product of density values
#' @param distribution character(1). delayed distribution family
#' @param bind character. parameter names that are bind together (i.e. equated) between both groups
#' @param twoPhase logical flag. Do we allow for two delay phases where event rate may change? Default is `FALSE`, i.e., a single delay phase.
#' @param ties character. How to handle ties within data of a group.
#' @param verbose integer flag. How much verbosity in output? The higher the more output. Default value is 0 which is no output.
#' @return the objective function (e.g., the negative MPSE criterion) for given choice of model parameters or `NULL` upon errors
objFunFactory <- function(x, y = NULL,
                          distribution = c("exponential", "weibull"), twoPhase = FALSE, bind = NULL,
                          method = c('MPSE', 'MLE0'), ties = c('density', 'equidist', 'random', 'error'),
                          verbose = 0L) {

  # setup ----
  stopifnot( is.numeric(x), length(x) > 0L, is.null(y) || is.numeric(y) && length(y) > 0L )
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  stopifnot( is.null(bind) || is.character(bind) && length(bind) >= 1L )
  ties <- match.arg(ties)

  twoPhase <- isTRUE(twoPhase)

  #standard ('original') names of distribution
  oNames <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, bind = NULL, transformed = FALSE)
  #bind: intersect with oNames enforces the canonical order of dist-parameters and drops unused parameters and empty strings
  bind <- intersect(oNames, bind)


  # data preparation ----

  # @param obs: data vector of one group
  # @return sorted, cleaned up data vector or NULL in case of trouble
  preprocess <- function(obs) {

    if ( is.null(obs) || ! is.numeric(obs)) return(NULL)

    ind_neg <- which(obs < 0L)
    if (length(ind_neg)){
      warning("Negative values in data", deparse(substitute(obs)), "! These are dropped.", call. = FALSE)
      obs <- obs[-ind_neg]
    }# fi

    # drop NA and +/-Inf & sort
    obs <- sort(obs[is.finite(obs)])

    if (!length(obs)) {
      warning("No valid data! Only non-negative and finite real values are valid.", call. = FALSE)
      return(invisible(NULL))
    }# fi

    if (is.null(obs)) return(NULL)

    # tie break
    # || ties == 'groupedML') # groupedML not implemented yet
    if ( startsWith(method, 'MLE') || ties == 'density' ) return(obs)

    diffobs <- diff(obs)
    stopifnot( all(diffobs >= 0L) ) # i.e. sorted obs

    tiesDiffInd <- which(diffobs == 0L) # < .Machine$double.xmin

    if (length(tiesDiffInd)){
      #rl <- rle(diff(tiesDiffInd))
      if (verbose > 0L){
        #length(which(rl$values > 1L))+1L,
        cat(glue('{length(tiesDiffInd) + sum(diff(tiesDiffInd)>1L) + 1L} tied observations ',
                 'in {sum(diff(tiesDiffInd)>1L) + 1L} group(s) within data vector.\n'))
      }

      if ( ties == 'error' ) stop('Ties within data are not allowed!', call. = FALSE)

      roundOffPrecision <- estimRoundingError(obs, maxObs = 1000L)
      if (verbose > 0L){
        cat(glue("Round-off error has magnitude {roundOffPrecision}\n"))
      }


      # rounding radius can't be wider than smallest observed diff.
      # plogis to mitigate the effect of sample size: the larger the sample the more we can 'trust' the observed minimal diff
      # obs[1L] = min(obs) = diff of minimal obs with 0
      rr <- .5 * min(stats::plogis(q = length(obs), scale = 11) * diffobs[which(diffobs > 0L)],
                     # rounding precision here
                     roundOffPrecision, obs[[1L]], na.rm = TRUE)

      ## modify tied observations per group of ties
      startInd <- endInd <- 1L
      repeat {
        #proceed to end of tie-group
        while (endInd < length(tiesDiffInd) && tiesDiffInd[endInd+1L] == tiesDiffInd[endInd] + 1L) {endInd <- endInd+1L}
        #include adjacent index to complete tie-group
        obsInd <- c(tiesDiffInd[startInd:endInd], tiesDiffInd[endInd]+1L)
        stopifnot( stats::sd(obs[obsInd]) == 0L ) #check: tie-group
        obs[obsInd] <- obs[obsInd] + if (ties == 'random') {
          # sort ensures that data after tie-break is still sorted from small to large
          sort(stats::runif(n = length(obsInd), min = -rr, max = +rr)) } else {
            stopifnot( ties == 'equidist' )
            # use evenly spaced observations to break tie as proposed by Cheng (1989) on Moran test statistic
            #+They first use the ties = 'density' approach for initial estimation of parameters for Moran's statistic
            seq.int(from = -rr, to = +rr, length.out = length(obsInd))
          }
        startInd <- endInd <- endInd+1L
        if ( startInd > length(tiesDiffInd) ) break
      } #repeat

      if (verbose > 1L && length(obs) < 50L ){
        cat(glue("New data: {paste(obs, collapse = ', ')}\n"))
      }
    } #fi tiesdiff

    # we have broken all ties
    stopifnot( !any(diff(obs)==0L) )

    obs
  } #fn preprocess

  # overwrite the data vectors with pre-processed data
  if (is.null({x <- preprocess(obs = x)})) return(invisible(NULL))
  y <- preprocess(obs = y)


  # do we have two groups after pre-processing?
  twoGroup <- isTRUE(!is.null(y) && is.numeric(y) && length(y))
  stopifnot( ! twoPhase || (distribution == 'exponential') ) #XXX twoPhase for Weibull not implemented yet!


  # checks ------------------------------------------------------------------

  if (!twoGroup && !is.null(bind) && length(bind)) {
    bind <- NULL
    warning("bind= has a given non-null argument but it is ignored as we have only a single group!",
            call. = FALSE)
  }

  if (startsWith(method, 'MLE') && twoGroup ) { #XXX not implemented yet!
    warning('MLE fitting is currently only supported for single group setting!', call. = FALSE)
    return(invisible(NULL))
  } # MLE


  # check that there is enough data (here we also look at bind= if twoGroup)
  if (!twoGroup && length(x) < length(oNames) || twoGroup && length(x) + length(y) < 2L * length(oNames) - length(bind) && min(length(x), length(y)) < length(oNames) - length(bind)) {
    warning('Too few valid observations provided!', call. = FALSE)
    return(invisible(NULL))
  }


  # optimization arguments -----

  # parameters as used inside the optimization function
  optpar_names <- getDist(distribution = distribution, type = "param", transformed = TRUE,
                          twoPhase = twoPhase, twoGroup = twoGroup, bind = bind)


  # get optimization start values and upper limits based on observations from a single group
  # for exponential distribution, the function respects the twoPhase-flag
  # @return list with transformed par for single group and upper limits for delay parameters, in canonical order (bind has no effect here!)
  getParSetting.gr <- function(obs){
    # contract: obs is sorted!

    DELAY_MIN <- 1e-9

    parV <- switch (EXPR = distribution,
                    # min(obs) = obs[1L]
                    exponential = {
                      parV0 <- c( max(DELAY_MIN, obs[[1L]] - 2/length(obs)),
                                  mean(obs - obs[[1L]] + 2/length(obs))**-1L )

                      # two extra parameters when exponential with *two* phases
                      if (twoPhase) parV0 <- c(parV0, obs[[floor(.5 + length(obs)/2L)]], parV0[[2L]])

                      parV0 <- purrr::set_names(parV0, nm = oNames)
                      # transformation of parameters into optfun-parametrization
                      extractPars(parV = parV0, distribution = "exponential", transform = TRUE)
                    },
                    weibull = {
                      stopifnot(! twoPhase) #XXXX 2-phase not supported yet for Weibull!
                      # start values from 'Weibull plot'
                      #+using the empirical distribution function
                      ## in MASS::fitdistr they simplify:
                      # lx <- log(x)
                      # m <- mean(lx)
                      # v <- var(lx)
                      # shape <- 1.2/sqrt(v)
                      # scale <- exp(m + 0.572/shape)
                      # take out extreme values for robustness (against possible outliers)
                      #+ when at least 40 observations
                      obs_f <- obs[max(1L, floor(length(obs)*.02)):ceiling(length(obs)*.975)] #assume sorted data
                      start_y <- log(-log(1-(seq_along(obs_f)-.3)/(length(obs_f)+.4)))
                      # cf. lm.fit(x = cbind(1, log(obs)), y = start_y))$coefficients
                      start_shape <- stats::cor(log(obs_f), start_y) * stats::sd(start_y) / stats::sd(log(obs_f))
                      start_scale <- exp(mean(log(obs_f) - mean(start_y) / start_shape))


                      parV0 <- c( max(DELAY_MIN, obs[[1L]] - 2/length(obs)),
                         start_shape,
                         start_scale )

                      parV0 <- purrr::set_names(parV0, nm = oNames)
                      # transformation of parameters into optfun-parametrization
                      extractPars(parV = parV0, distribution = "weibull", transform = TRUE)
                    },
                    # default:
                    stop(glue("Provided distribution {sQuote(distribution)} is not implemented!"), call. = FALSE)
    )

    list(
      par = parV,
      delay1_upper = max(DELAY_MIN, obs[[1L]] - .01/length(obs), obs[[1L]]*.9999),
      delay2_upper = log(max(DELAY_MIN, obs[[length(obs)]] - .02/length(obs), obs[[length(obs)]]*.999))
    )
  }# fn getParSetting.gr


  # parameter bounds: set lower & upper bounds
  lowerB <- upperB <- purrr::set_names(rep_len(NA_real_, length(optpar_names)),
                                       nm = optpar_names)

  PAR_BOUNDS <- list(delay1 = c(lower = 0, upper = NA_real_),
                     delay2 = c(lower = -Inf, upper = NA_real_),
                     rate  = c(lower = -Inf, upper = +Inf),
                     shape = c(lower = sqrt(.Machine$double.eps), upper = +Inf), ##XXXX for the time being, as we are currently untransformed in Weibull ## on log would be: -Inf
                     scale = c(lower = sqrt(.Machine$double.eps), upper = +Inf)) ##XXXX scale = -Inf


  # alas, purrr::iwalk did not work for me here
  for (nam in names(PAR_BOUNDS)) {
    idx <- startsWith(optpar_names, prefix = nam)
    if (any(idx)) {
      lowerB[idx] <- purrr::chuck(PAR_BOUNDS, nam, 'lower')
      upperB[idx] <- purrr::chuck(PAR_BOUNDS, nam, 'upper')
    } #fi
  } #rof



  par0_x <- getParSetting.gr(x)
  # set parameter vector for group 1 and finish upper bound: match delay1 and delay2
  parV <-
    if (! twoGroup) {
      upperB[['delay1_tr']]  <- par0_x[['delay1_upper']]
      if (twoPhase) upperB[['delay2_tr']] <- par0_x[['delay2_upper']]

      par0_x[['par']]

    } else { #twoGroup

      # all parameters are bound
      if ( length(bind) == length(oNames) ) {

        # treat x and y as a single group for upper limit & start value heuristic
        par0_xy <- getParSetting.gr(c(x,y))

        upperB['delay1_tr'] <- par0_xy[['delay1_upper']]
        if (twoPhase) upperB[['delay2_tr']] <- par0_xy[['delay2_upper']]

        par0_xy[['par']]

      } else { #twoGroup, not all params bound!

        par0_y <- getParSetting.gr(y)

        start_x <- par0_x[['par']]
        start_y <- par0_y[['par']]

        # set upper bound for delay parameter(s)!
        if ('delay1' %in% bind) {
          upperB['delay1_tr'] <- min(par0_x[['delay1_upper']], par0_y[['delay1_upper']])
        } else {
          upperB['delay1_tr.x'] <- par0_x[['delay1_upper']]
          upperB['delay1_tr.y'] <- par0_y[['delay1_upper']]
        } # fi

        if (twoPhase) {
          if ('delay2' %in% bind){
            upperB[['delay2_tr']] <- max(par0_x[['delay2_upper']], par0_y[['delay2_upper']])
          } else {
            upperB['delay2_tr.x'] <- par0_x[['delay2_upper']]
            upperB['delay2_tr.y'] <- par0_y[['delay2_upper']]
          }
        }

        # return start value
        if ( is.null(bind) ){ # two groups unbound
          c(start_x, start_y)
        } else if (distribution == 'exponential') {

          switch( EXPR = paste(bind, collapse = '+'),
                  delay1 = {
                    c(min(start_x[[1L]], start_y[[1L]]), start_x[-1L], start_y[-1L])
                  },
                  rate1 = {
                    # arithmetic mean of log-rates
                    c((start_x[[2L]] + start_y[[2L]]) / 2L,
                      start_x[-2L], start_y[-2L])
                  },
                  delay2 = {
                    stopifnot(twoPhase)
                    # arithmetic mean of log-scaled delay2-start values
                    c((start_x[[3L]] + start_y[[3L]]) / 2L,
                      start_x[-3L], start_y[-3L])
                  },
                  rate2 = {
                    stopifnot(twoPhase)
                    # arithmetic mean of log-rates
                    c((start_x[[4L]] + start_y[[4L]]) / 2L,
                      start_x[-4L], start_y[-4L])

                  },
                  `delay1+rate1` = {
                    stopifnot(twoPhase) #otherwise, we would have bound all parameters which is already handled above!
                    # min for delay1, mean for rate1 on log-scale
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[2L]] + start_y[[2L]]) / 2L,
                      start_x[-c(1L, 2L)], start_y[-c(1L, 2L)])
                  },
                  `delay1+delay2` = {
                    stopifnot(twoPhase)
                    # min for delay1, mean for log-scaled delay2
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[3L]] + start_y[[3L]])/2L,
                      start_x[-c(1L, 3L)], start_y[-c(1L, 3L)])
                  },
                  `delay1+rate2` = {
                    stopifnot(twoPhase)
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[4L]] + start_y[[4L]]) / 2L,
                      start_x[-c(1L, 4L)], start_y[-c(1L, 4L)])
                  },
                  `rate1+delay2` = {
                    stopifnot(twoPhase)
                    c((start_x[[2L]]+start_y[[2L]]) / 2L, (start_x[[3L]]+start_y[[3L]]) / 2L,
                      start_x[-c(2L, 3L)], start_y[-c(2L, 3L)])
                  },
                  `rate1+rate2` = {
                    stopifnot(twoPhase)
                    c((start_x[[2L]]+start_y[[2L]]) / 2L, (start_x[[4L]]+start_y[[4L]]) / 2L,
                      start_x[-c(2L, 4L)], start_y[-c(2L, 4L)])
                  },
                  `delay2+rate2` = {
                    stopifnot(twoPhase)
                    c((start_x[[3L]]+start_y[[3L]]) / 2L, (start_x[[4L]]+start_y[[4L]]) / 2L,
                      start_x[-c(3L, 4L)], start_y[-c(3L, 4L)])
                  },
                  `delay1+rate1+delay2` = {
                    stopifnot(twoPhase)
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[2L]]+start_y[[2L]]) / 2L, (start_x[[3L]]+start_y[[3L]]) / 2L,
                      start_x[[4L]], start_y[[4L]])
                  },
                  `delay1+rate1+rate2` = {
                    stopifnot(twoPhase)
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[2L]]+start_y[[2L]]) / 2L, (start_x[[4L]]+start_y[[4L]]) / 2L,
                      start_x[[3L]], start_y[[3L]])
                  },
                  `delay1+delay2+rate2` = {
                    stopifnot(twoPhase)
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[3L]]+start_y[[3L]]) / 2L, (start_x[[4L]]+start_y[[4L]]) / 2L,
                      start_x[[2L]], start_y[[2L]])
                  },
                  `rate1+delay2+rate2` = {
                    stopifnot(twoPhase)
                    c((start_x[[2L]]+start_y[[2L]]) / 2L, (start_x[[3L]]+start_y[[3L]]) / 2L, (start_x[[4L]]+start_y[[4L]]) / 2L,
                      start_x[[1L]], start_y[[1L]])
                  },
                  # default
                  stop("This bind-configuration is not implemented!", call. = FALSE)
          )

        } else {
          stopifnot( distribution == 'weibull' )

          switch( EXPR = paste(bind, collapse = '+'),
                  delay1 = {
                    c(min(start_x[[1L]], start_y[[1L]]), start_x[-1L], start_y[-1L])
                  },
                  shape1 = {
                    # geometric mean of shapes
                    c(sqrt(start_x[[2L]] * start_y[[2L]]), start_x[-2L], start_y[-2L])
                  },
                  scale1 = {
                    # arithmetic mean of scales (corresponds to harmonic mean of rates)
                    c((start_x[[3L]] + start_y[[3L]])/2L, start_x[-3L], start_y[-3L])
                  },
                  `delay1+shape1` = {
                    c(min(start_x[[1L]], start_y[[1L]]), sqrt(start_x[[2L]] * start_y[[2L]]), start_x[[3]], start_y[[3]])
                  },
                  `delay1+scale1` = {
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[3L]] + start_y[[3L]])/2L, start_x[[2]], start_y[[2]])
                  },
                  `shape1+scale1` = {
                    c(sqrt(start_x[[2L]] * start_y[[2L]]), (start_x[[3L]] + start_y[[3L]])/2L, start_x[[1L]], start_y[[1L]])
                  },
                  # default
                  stop("This bind-configuration is not implemented!", call. = FALSE)
          )

        } # weibull
      } #twoGrp, not all params bound!
    } # twoGrp

  # ensure we have names of transformed parameters
  parV <- purrr::set_names(parV, nm = optpar_names)

  stopifnot( ! any(is.na(lowerB), is.na(upperB)) )
  # clean up env # or use local() XXX
  remove(list=c("PAR_BOUNDS", "par0_x"))


  optim_args <- list(
    par = parV,
    method = "L-BFGS-B",
    lower = lowerB,
    upper = upperB,
    # Most parameters are on log-scale.
    control = list(parscale = scalePars(parV, lowerB = 1e-3, upperB = 1e3))
  )




  # objective function ----

  # log spacings:
  # calculate the differences in EDF (for given parameters in group) of adjacent observations on log scale
  # @return n+1 cumulative diffs on log-scale (or single negative number in twoPhase when delay2 <= delay in quick fix)
  getCumDiffs <- function(pars, group) {
    # get() would (by default) start looking in the execution environment of getCumDiffs and would find the object in its parent
    # or: env = rlang::env_parent() # parent of caller env is (normally!?) the function-env
    # or: env = rlang::fn_env(getCumDiffs) # referring directly to the function environment (but requires function obj)
    obs <- rlang::env_get(env = rlang::env_parent(rlang::current_env(), n=1L), nm = group, inherit = FALSE)
    # back-transform parameters to original scale (for CDF)
    pars.gr <- extractPars(pars, distribution = distribution, group = group, transform = TRUE)

    if (verbose > 1L){
      cat(glue("Parameter vector for group {group} after back-transformation: ",
               "{paste(names(pars.gr), round(pars.gr, 2), sep = ': ', collapse = ', ')}"), "\n")
    }

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
      stopifnot( ties == 'density' ) # other tie-strategies have already dealt with ties in *preprocess*
      # increase index by 1 to get from diff(obs)-indices to cumDiffs-indices
      cumDiffs[1L+ind_t] <- purrr::exec(getDist(distribution, type = "dens"), !!! c(list(x = obs[ind_t]), pars.gr))
    } #fi

    # respect the machine's numerical lower limit
    cumDiffs[which(cumDiffs < .Machine$double.xmin)] <- .Machine$double.xmin

    log(cumDiffs)

  }# fn getCumDiffs


  # Objective function like negative mean log-spacings for MPSE or negative log-likelihood for MSE0
  # Estimate parameters by minimizing this function.
  # param `pars` the parameter vector.
  # param `aggregated` a logical flag. For two group case, `FALSE` returns individual mean log cum-diffs per group
  objFun <- function(pars, aggregated = TRUE) {
    stopifnot( length(pars) == length(optpar_names) )
    pars <- purrr::set_names(pars, optpar_names) # extractPars in getCumDiffs needs names!

    switch(method,
           MPSE = {
             - if (! twoGroup) {
               mean(getCumDiffs(pars, group = "x"))
             } else {
               #twoGroup:
               #the approach to first merge x and y and then do the cumDiffs, log and mean does *not* work out
               #because the parameters should be optimized within group.
               #merged data lead to frequent non-convergence or visually bad fits
               res <- c(mean(getCumDiffs(pars, group = "x")), mean(getCumDiffs(pars, group = "y")))

               if (aggregated) stats::weighted.mean(res, w = c(length(x), length(y))) else res
             }
           },
           MLE0 = {
             # MLE0 is currently only implemented for single group situation
             stopifnot( ! twoGroup )
             stopifnot( ! twoPhase ) #XXX not implemented yet!
             nObs <- length(x)
             xc <- x - pars[['delay1']]
             - if (distribution == 'exponential')
               nObs * (log(pars[['rate1']]) - pars[['rate1']] * mean(xc)) else
                 nObs * (log(pars[['shape1']]) - pars[['shape1']] * log(pars[['scale1']]) + (pars[['shape1']]-1L) * mean(log(xc)) - mean(xc**pars[['shape1']])/pars[['scale1']]**pars[['shape1']])
           },
           stop(glue('Objective function for method {method} is not implemented!'), call. = FALSE)
    )
  }

  # attach analytical solution for MLE
  if ( method == 'MLE0' && ! twoGroup && ! twoPhase && distribution == 'exponential' ){
    attr(objFun, which = "opt") <- list(par = extractPars(parV = c(delay1 = x[[1L]], rate1 = 1L/(mean(x) - x[[1L]])),
                                                          distribution = 'exponential', transform = TRUE), #expect transformed parameters
                                        value = length(x) * ( log(mean(x) - x[[1L]]) + 1L ),
                                        convergence = 0L,
                                        message = "analytic solution for standard MLE ('MLE0')",
                                        counts = 1L)
  }

  objFun
}



#' Fit optimal parameters according to the objective function (either MPSE or MLE0).
#'
#' The objective function carries the given data in its environment and it is to be minimized.
#' R's standard routine `stats::optim` does the numerical optimization, using numerical derivatives.
#' or the analytical solution is returned directly if available.
#' @param objFun objective function to be minimized
#' @param optim_args list of own arguments for optimization. If `NULL` it uses the default optim arguments associated to the objective function.
#' @param verbose integer that indicates the level of verboseness. Default 0 is quiet.
#' @return optimization object including a named parameter vector or `NULL` in case of errors during optimization
delay_fit <- function(objFun, optim_args = NULL, verbose = 0) {

  if (is.null(objFun)) return(invisible(NULL))
  stopifnot( is.function(objFun) )
  objFunEnv <- rlang::fn_env(objFun)

  # check if there is already a solution provided by the objective function
  optObj <- attr(objFun, which = "opt", exact = TRUE)

  if ( is.list(optObj) && all( c('par', 'value', 'convergence') %in% names(optObj)) ){
    if (verbose > 0L) cat("Using provided (analytical) solution to objective function.\n")
  } else {
    optObj <- NULL #start from scratch
    # numeric optimization
    if (verbose > 0L) message("Start with numeric optimziation of objective function.")

    if (is.null(optim_args)) optim_args <- rlang::env_get(objFunEnv, nm = "optim_args")
    stopifnot( is.list(optim_args), 'par' %in% names(optim_args) )
    # set objective function (overwrite entry 'fn' if it is already present)
    optim_args[["fn"]] <- objFun

    # optim: 1st attempt
    try({optObj <- purrr::exec(stats::optim, !!! optim_args)}, silent = TRUE)

    if (is.null(optObj)){
      if (verbose > 0L) warning(glue("{rlang::env_get(env = objFunEnv, nm = 'method')}-optimization failed during model fit!"),
                                call. = FALSE)
    } else if ( isTRUE(optObj$convergence > 0L) ){
      # do a 2nd attempt of optim in case it did not converge in the first place
      if (verbose > 1L) message("No proper convergence during 1st optimization in delay fit. Re-try with different parameter scaling.")

      # Use parameter values of non-converged fit as new start values (and adapt parscale accordingly)
      #+The objFun is to be minimized,  smaller is better!
      if ( isTRUE(is.numeric(optObj$par) && all(is.finite(optObj$par)) && optObj$value < objFun(optim_args$par)) ){
        optim_args[['par']] <- optObj$par  # purrr::assign_in(where = "par", value = optObj$par)

        if ( isTRUE("parscale" %in% names(optim_args[["control"]])) ){
          optim_args[['control']][['parscale']] <- scalePars(optim_args[['par']])
        }

        # optim: 2nd attempt
        optObj <- NULL
        try({optObj <- purrr::exec(stats::optim, !!! optim_args)}, silent = TRUE)

        if ( is.null(optObj) || isTRUE(optObj$convergence > 0L && verbose > 0L) ) warning("No proper convergence after re-try.",
                                                                                          call. = FALSE)
      }## fi rescaling for 2nd attempt
    }## fi 2nd attempt necessary?

    # set names to parameter vector
    if (! is.null(optObj)){
      stopifnot( 'par' %in% names(optObj) )
      # set canonical names for parameters
      optObj$par <- purrr::set_names(optObj$par, rlang::env_get(env = objFunEnv, nm = "optpar_names"))
      # save optim_args in optimization object (but w/o objective function)
      optim_args$fn <- NULL
      optObj <- append(optObj, values = list(optim_args = optim_args))
    }
  } #esle numeric optimization

  optObj
}



#' Fit a delayed Exponential or Weibull model to one or two given sample(s).
#'
#' Maximum product of spacings estimation is used by default to fit the parameters. Estimation via standard maximum likelihood (`method = 'MLE0`) is available, too,
#' but MLE yields biased estimates.
#'
#' Numerical optimization is done by `stats::optim`.
#' @param x numeric. observations of 1st group. Can also be a list of data from two groups.
#' @param y numeric. observations from 2nd group
#' @param distribution character. Which delayed distribution is assumed? Exponential or Weibull.
#' @param method character. Which method to fit the model? 'MPSE' = maximum product of spacings estimation *or* 'MLE0' = standard maximum likelihood estimation
#' @param bind character. parameter names that are bind together in 2-group situation.
#' @param ties character. How to handle ties.
#' @param optim_args list. optimization arguments to use. Use `NULL` to use the data-dependent default values.
#' @param verbose integer. level of verboseness. Default 0 is quiet.
#' @return `incubate_fit` the delay-model fit object. Or `NULL` if optimization failed (e.g. too few observations).
#' @export
delay_model <- function(x = stop('Specify observations for at least one group x=!', call. = FALSE), y = NULL,
                        distribution = c("exponential", "weibull"), twoPhase = FALSE,
                        bind = NULL, ties = c('density', 'equidist', 'random', 'error'),
                        method = c('MPSE', 'MLE0'),
                        optim_args = NULL, verbose = 0) {


  # setup -------------------------------------------------------------------

  if (is.logical(verbose)) verbose <- as.numeric(verbose)
  if ( is.null(verbose) || ! is.numeric(verbose) || ! is.finite(verbose) ) verbose <- 0L
  verbose <- verbose[[1L]]


  # unpack x if it is a list of two vectors
  if (is.list(x)){
    stopifnot( length(x) == 2L )
    y <- x[[2L]]
    x <- x[[1L]]
  }

  # enforce that the first argument x= is properly instantiated
  stopifnot( !is.null(x) && is.numeric(x) && length(x) )

  distribution <- match.arg(distribution)

  method <- toupper(method)
  if (length(method) == 1L && method == 'MSE') {
    message("The method name 'MPSE' is prefered over the previously used name 'MSE'!")
    method <- 'MPSE'
  }
  method <- match.arg(method)
  ties <- match.arg(ties)

  if (is.character(bind)){
    if (any(grepl(pattern = "_tr", bind, fixed = TRUE))){
      stop("Parameter names to bind= refer to the distribution parameters and not to the transformed parameters of the objective function.", call. = FALSE)
    }

    # translate convenience names (for single phase) to canonical names
    unNmbrdIdx <- !grepl(pattern = "[12]", bind, fixed = FALSE)
    if (any(unNmbrdIdx)){
      bind[unNmbrdIdx] <- paste0(bind[unNmbrdIdx], "1") #interpret un-numbered parameters as referring to phase 1
      if (verbose > 0L) cat("The unnumbered parameter names in bind= are taken to refer to initial phase and are translated to canonical parameter names.\n")
    }

  }#fi bind=


  # objective function ------------------------------------------------------

  objFun <- objFunFactory(x = x, y = y, method = method, distribution = distribution,
                          twoPhase = twoPhase, bind = bind, ties = ties, verbose = verbose)
  objFunEnv <- rlang::fn_env(objFun)

  # optimise objective function
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  twoGroup <- rlang::env_get(env = objFunEnv, nm = "twoGroup")
  # overwrite data with  pre-processed data
  x <- rlang::env_get(env = objFunEnv, nm = "x")
  y <- rlang::env_get(env = objFunEnv, nm = "y", default = NULL)



  structure(
    list(
      data = if (twoGroup) list(x = x, y = y) else x,
      distribution = distribution,
      twoPhase = twoPhase,
      twoGroup = twoGroup,
      method = method,
      bind = bind,
      ties = ties,
      objFun = objFun,
      par = extractPars(optObj$par, distribution = distribution, transform = TRUE), # /!\ keep in sync with update()!
      val = optObj$value,
      optimizer = purrr::compact(c(list(parOpt = optObj$par), optObj[c('convergence', 'message', 'counts', 'optim_args')]))),
    class = "incubate_fit")
}

#' @export
print.incubate_fit <- function(x, ...){
  coe <- coef(x)
  cat(glue::glue_data(x, .sep = "\n",
                      "Fit a delayed {distribution} {c('', 'with two delay phases')[[1L+twoPhase]]} through {c('Maximum Product of Spacings Estimation (MPSE)', 'standard Maximum Likelihood Estimation (MLE0)')[[1L+(method=='MLE0')]]} for {c('a single group', 'two independent groups')[[1L+twoGroup]]}.",
                      "Data: {if (twoGroup) paste(lengths(data), collapse = ' and ') else length(data)} observations, ranging from {paste(signif(range(data), 4), collapse = ' to ')}",
                      "Fitted coefficients: {paste(paste('\n  ', names(coe)), signif(coe,5L), sep = ': ', collapse = ' ')}\n\n")
  )
}

#' Coefficients of a delay-model fit.
#' @param object object that is a `incubate_fit`
#' @param transformed flag. Do we request the transformed parameters as used within the optimization?
#' @param group character string to request the canonical parameter for one group
#' @param ... further arguments, currently not used.
#' @return named coefficient vector
#' @export
coef.incubate_fit <- function(object, transformed = FALSE, group = NULL, ...){
  #stopifnot( inherits(object, "incubate_fit") )
  transformed <- isTRUE(transformed)
  stopifnot( object[["distribution"]] == 'exponential' || ! transformed ) #XXX currently transformation only implemented for 'expon'

  extractPars(parV = if (transformed) purrr::chuck(object, "optimizer", "parOpt") else object[["par"]],
              distribution = object[["distribution"]], group = group, transform = FALSE)
}

#' @export
summary.incubate_fit <- function(object, ...){
  print(object)
}

#' Refit an `incubate_fit`-object with specified optimization arguments.
#' If more things need to be changed go back to `delay_model` and start from scratch.
#' @param object `incubate_fit`-object
#' @param optim_args optimization arguments
#' @param verbose integer flag. Requested verbosity during `delay_fit`
#' @param ... further arguments, currently not used.
#' @return The updated fitted object of class `incubate_fit`
#' @export
update.incubate_fit <- function(object, optim_args, verbose = 0, ...){

  stopifnot( all(c("objFun", "twoPhase", "twoGroup", "data", "par", "val", "optimizer") %in% names(object)) )

  ## fit model with given optim_args
  optObj <- delay_fit(object[["objFun"]], optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  # update all relevant fields in the list
  object[c("par", "val", "optimizer")] <- list(par = extractPars(optObj$par, distribution = object[["distribution"]], transform = TRUE),
                                               val = optObj$value,
                                               # drop NULLs from list (e.g. if optim_args is not present)
                                               optimizer = purrr::compact(c(
                                                 list(parOpt = optObj$par),
                                                 optObj[c("convergence", "message", "counts", "optim_args")])))

  object
}

#' @export
plot.incubate_fit <- function(x, y, title, subtitle, ...){
  stopifnot( inherits(x, "incubate_fit") )
  # parameter y comes from the plot-generic. y is not used here.

  rlang::check_installed(pkg = 'ggplot2', reason = 'to get plots', version = '3.3')

  cumFun <- getDist(x[["distribution"]], type = "cdf")

  p <- grNames <- NULL

  # catch the one-group case!
  if ( isTRUE(x[["twoGroup"]]) ){
    stopifnot( is.list(x[["data"]]) )

    grNames <- names(x[["data"]])

    p <- ggplot2::ggplot(data = tibble::enframe(unlist(x[["data"]]), name = "group"),
                         mapping = ggplot2::aes_(x = ~value, col = ~substr(x = group, 1L, 1L))) +
      # add estimated delay model(s)
      purrr::map(.x = grNames,
                 .f = ~ ggplot2::geom_function(mapping = ggplot2::aes(col = .x), inherit.aes = FALSE,
                                               fun = cumFun,
                                               args = coef(x, group = .x), linetype = "dashed"))
  } else {
    grNames <- "x"
    p <- ggplot2::ggplot(data = tibble::tibble(value=x[["data"]]),
                         mapping = ggplot2::aes_(x = ~value)) +
      # add estimated delay model
      ggplot2::geom_function(inherit.aes = FALSE,
                             fun = cumFun,
                             args = coef(x, group = grNames), linetype = "dashed")
  }


  if (missing(title)) title <- glue::glue_data(x,
                                               "Fitted {distribution} delay {c('model ', 'models ')[[1L+twoGroup]]}",
                                               "{c('', 'with two delay phases')[[1L+twoPhase]]}")
  if (missing(subtitle)) subtitle <- paste(purrr::map_chr(grNames,
                                                          ~ paste(names(coef(x, group = .)), signif(coef(x, group = .), 4),
                                                                  sep = ": ", collapse = ", ")),
                                           collapse = " - ")


  p +
    ggplot2::stat_ecdf(pad=TRUE) +
    ggplot2::xlim(0L, NA) +
    ggplot2::labs(x = 'Time', y = 'Cumulative prop. of events',
                  col = 'Group',
                  title = title, subtitle = subtitle) +
    ggplot2::scale_y_reverse()
}


#' @export
simulate.incubate_fit <- function(object, nsim = 1, seed = NULL, ...){
  stopifnot(inherits(object, 'incubate_fit'))

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

#' Generate bootstrap distribution of model parameters to fitted incubate model.
#'
#' Bootstrap data are here estimated coefficients from models fitted to bootstrap samples.
#' The bootstrap data is used to make bootstrap inference in the second step.
#' It is an internal function, the main entry point is [confint.incubate_fit()].
#' @param object an `incubate_fit`-object
#' @param bs_data character. Which type of bootstrap method to generate data?
#' @param R integer. Number of bootstrapped model coefficient estimates
#' @param useBoot flag. Do you want to use the boot-package? Default value is `FALSE`.
#' @param smd_factor numeric. smooth-delay factor: influence the amount of smoothing. 0 means no smoothing at all. Default is 0.25 (as was optimal in simulation for log-quantile together with log-delay-shift = 5)
#' @return bootstrap data, either as matrix or of class `boot` (depending on the `useBoot`-flag)
bsDataStep <- function(object, bs_data = c('parametric', 'ordinary'), R, useBoot = FALSE, smd_factor = 0.25) {
  bs_data <- match.arg(bs_data)
  stopifnot(is.numeric(R), length(R) == 1L, R > 1L)
  R <- ceiling(R)
  useBoot <- isTRUE(useBoot)
  ranFun <- getDist(object$distribution, type = "r")
  dFun <- getDist(object$distribution, type = "d")
  twoGroup <- isTRUE(object$twoGroup)
  nObs <- if (twoGroup) lengths(object$data) else length(object$data)
  # full untransformed parameter vector
  coefVect <- coef.incubate_fit(object, group = NULL, transformed = FALSE)
  del1_ind <- grep('delay1', names(coefVect)) # indices of coefficients that involve delay1, e.g. 'delay1' or 'delay1.y'
  ncoef <- length(coefVect)

  stopifnot( ncoef > 0L )

  stopifnot( is.numeric(smd_factor), length(smd_factor) == 1L, smd_factor >= 0L )
  smoothDelay <- isTRUE(smd_factor > 0L)

  if (smoothDelay && bs_data != 'parametric') {
    smoothDelay <- FALSE
    smd_factor <- 0L
    # how could smooth_delay work also for ordinary bootstrap?!
    warning('Smoothing of delay is only implemented for parametric bootstrap!', call. = FALSE)
  }


  # smooth first delay: sample delay values according to objective function (where delay is varied and other parameters are kept fixed) in the vicinity of the estimated first delay
  # This reflects the certainty we have in the delay estimation.
  # Low variability in event time data (or high sample size) will lead to a quickly deteriorating objective function.
  # return vector of length R with candidate values for first delay
  getSMDCandidates <- function(group = 'x'){
    obs <- if (twoGroup) object$data[[group]] else object$data
    obs1 <- obs[[1L]]
    del_coef <- coef.incubate_fit(object, transformed = FALSE, group = group)[['delay1']]

    # avoid smoothing if 1st observation or estimated delay is too close to zer0
    if ( min(obs1, del_coef) < sqrt(.Machine$double.eps) ) return(rep_len(del_coef, length.out = R))

    stopifnot( is.function(object$objFun) )

    groupIdx <- 1L + (twoGroup && group == 'y')
    # in case of a delay per group ('delay.x' and 'delay.y') use the right one
    if (length(del1_ind) > 1L) del1_ind <- del1_ind[[groupIdx]]


    # look at differences of first observations
    obs_d <- diff(obs[seq_len(min(23L, nObs[[groupIdx]]))])
    obs_d <- obs_d[is.finite(obs_d) & obs_d > 0L] #get rid of ties
    obs_d <- if (! length(obs_d)) .0001 else min(obs_d)

    # candidate region for delay parameters
    #+min(..) ensures that we are not too close at obs1, otherwise for MLE we have only a single point
    #+ del_coef - (obs1 - del_coef) = 2 * del_coef - obs1
    del_interv <- c(low = max(0L, min(del_coef - (obs1 - del_coef), del_coef - obs_d,
                                      obs1 - .0001, obs1 * .9999, na.rm = TRUE)),
                    high = obs1)

    #+areas for delay with high values of objective function are more likely to be sampled
    #+candidate region: symmetric around coef_del as midpoint, up to smallest observed value
    #+candidate region becomes finer sampled the broader the interval is
    #+point estimate for delay is part of sample (if lower bound is not cut to be 0, via max in from= argument)
    delayCandDF <- tibble(
      delay1 = seq.int(from = del_interv[['low']], to = del_interv[['high']],
                       # uneven number of grid points (hence, MPSE-estimate for delay will be one of the grid points)
                       # grid step width at most 0.005
                       length.out = max(997L, 2L * min(ceiling(R/2), 100L*ceiling(diff(del_interv)))+1L)),
      # fixing all parameter estimates other than delay1
      objVal = purrr::map_dbl(.x = .data[["delay1"]],
                              # objective function expects transformed parameters for all involved groups
                              # delay1 is normally transformed via ident (for both exponential and Weibull),
                              #+hence we could spare us the round trip to canonical parameters and to replace there and to transform back for optimization
                              #+but this round-trip approach is more safe (in case we ever introduce a non-ident transformation for delay1)
                              .f = ~ object$objFun(pars = extractPars(parV = replace(coefVect, del1_ind, .x),
                                                                      distribution = object$distribution, group = NULL, transform = TRUE),
                                                   aggregated = FALSE)[[groupIdx]])
    )

    # we like to drop last entry (delay = 1st observation) as objective function tends to explode
    # but we have to keep last entry if it corresponds to the delay estimate (e.g., as is the case for MLE0-fitting)
    if (delayCandDF$delay1[NROW(delayCandDF)] > del_coef) delayCandDF <- delayCandDF[-NROW(delayCandDF),, drop = FALSE]
    # relative change to optimal value, will be negative as objective function is minimized
    delayCandDF$objValInv <- (object$val - delayCandDF$objVal) / (object$val+.01)
    # shift upwards into non-negative area
    delayCandDF$objValInv <- delayCandDF$objValInv - min(delayCandDF$objValInv, na.rm = TRUE)
    # scale to be between 0 and 1.
    # small smd_factor => high exponent => peaked distribution
    delayCandDF$objValInv <- (delayCandDF$objValInv / (max(delayCandDF$objValInv, na.rm = TRUE) + .001))**(1L/(smd_factor+.01))
    delayCandDF$cumSum0 <- cumsum(delayCandDF$objValInv)
    # scale cumSum0 to 1.
    delayCandDF$cumSum <- delayCandDF$cumSum0 / max(delayCandDF$cumSum0)
    # lag-1: have it start with 0 and end with a single 1 (the last cumSum is most often 0 as largest delay value has typically objValInv = 0)
    delayCandDF$cumSum <- c(0L, delayCandDF$cumSum[-NROW(delayCandDF)])

    # rightmost.closed = TRUE for the unlikely/impossible?! case that we draw a 1 by runif
    delayCandDF$delay1[findInterval(x = stats::runif(R), vec = delayCandDF$cumSum, rightmost.closed = TRUE)]
  }

  delayCandX <- if (smoothDelay) getSMDCandidates(group = 'x')
  delayCandY <- if (smoothDelay && twoGroup) getSMDCandidates(group = 'y')

  if (useBoot) {
    stopifnot(!twoGroup) # for the time being only single group calls are supported!
    boot::boot(data = object$data,
               statistic = function(d, i) coef(delay_model(x=d[i], distribution = object$distribution, twoPhase = object$twoPhase,
                                                           ties = object$ties,
                                                           method = object$method, bind = object$bind), transformed = FALSE),
               sim = bs_data, mle = coef(object), R = R,
               ran.gen = function(d, coe){ # ran.gen function is only used for parametric bootstrap
                 if (smoothDelay){
                   coe[['delay1']] <- delayCandX[sample.int(n = R, size = 1L)]
                 }
                 rlang::exec(ranFun, !!! as.list(c(n=nObs[[1L]], coe)))
               })

  } else { # no boot-library
    # own implementation: we inline data generation (simulate) and model fitting in one function
    # get coefficients from bootstrapped data
    #+(either by ordinary bootstrap of data or by parametric bootstrap)
    coefBSFun <- switch(bs_data,
                        ordinary = function(dummy) {
                          # draw bootstrap samples from the data
                          x <- (if (twoGroup) object$data$x else object$data)[sample.int(n = nObs[[1L]], replace = TRUE)]
                          y <- if (twoGroup) object$data$y[sample.int(n = nObs[[2L]], replace = TRUE)]

                          coef(delay_model(x=x, y=y, distribution = object$distribution, twoPhase = object$twoPhase,
                                           ties = object$ties,
                                           method = object$method, bind = object$bind), transformed = FALSE)
                        },
                        parametric = {
                          # generate data from the fitted model
                          # for performance reasons, we 'inline' the simulate code, cf. test_diff

                          # arguments to the random function generation
                          ranFunArgsX <- as.list(c(n=nObs[[1L]], coef.incubate_fit(object, transformed = FALSE, group = "x")))
                          ranFunArgsY <- if (twoGroup) as.list(c(n=nObs[[2L]], coef.incubate_fit(object, transformed = FALSE, group = "y")))

                          function(ind) {
                            if (smoothDelay){
                              #+smooth delay according to how sure are we about the delay-estimate:
                              #+the more sure the smaller is the smoothing
                              ranFunArgsX[['delay1']] <- delayCandX[ind]
                              if (twoGroup) ranFunArgsY[['delay1']] <- delayCandY[ind]
                            }

                            # cf simulate (but inlined here for performance reasons)
                            x <- rlang::exec(ranFun, !!! ranFunArgsX)
                            y <- if (twoGroup) rlang::exec(ranFun, !!! ranFunArgsY)

                            dm <- suppressWarnings(delay_model(x=x, y=y, distribution = object$distribution, twoPhase = object$twoPhase,
                                                               ties = object$ties,
                                                               method = object$method, bind = object$bind))
                            retVec <- rep.int(NA_real_, times = ncoef)
                            if (! is.null(dm) && inherits(dm, "incubate_fit")) retVec <- coef.incubate_fit(dm, transformed = FALSE)

                            retVec
                          }
                        },
                        stop('Unkown bootstrap data generation type!')
    )

    # add originally fitted coefficients as first column!
    retM <- cbind(coef(object),
                  future.apply::future_vapply(X = seq_len(R), FUN.VALUE = numeric(length = ncoef),
                                              FUN = coefBSFun, future.seed = TRUE))

    # drop columns that contain NA-values (as bootstrap coefficient estimates)
    retM <- retM[,!.colSums(!is.finite(retM), m = ncoef, n = R+1L)]

    # return at most R columns
    retM[, seq_len(min(R, NCOL(retM)))]

    # more clear and shorter but less efficient!
    # future.apply::future_vapply(simulate(object, nsim = R), FUN.VALUE = numeric(length(cf)),
    #  FUN = \(d) coef(delay_model(x=d, distribution = object$distribution, ties = object$ties, method = object$method, bind = object$bind)))

  }
}

#' Confidence intervals for parameters of incubate-model fits.
#'
#' Bias-corrected bootstrap confidence limits (either quantile-based or normal-approximation based) are generated.
#' Optionally, there are also variants that use a log-transformation first.
#' At least R=1000 bootstrap replications are recommended. Default are quantile-based confidence intervals that internally use a log-transformation.
#' @param object object of class `incubate_fit`
#' @param parm character. Which parameters to get confidence interval for?
#' @param level numeric. Which is the requested confidence level for the interval? Default value is 0.95
#' @param R number of bootstrap replications. Used only if not `bs_data`-object is provided.
#' @param bs_data character or bootstrap data object. If character, it specifies which type of bootstrap is requested and the bootstrap data will be generated. Data can also be provided here directly. If missing it uses parametric bootstrap.
#' @param bs_infer character. Which type of bootstrap inference is requested to generate the confidence interval?
#' @param useBoot logical. Delegate bootstrap confint calculation to the `boot`-package?
#' @param ... further arguments, currently not used.
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.incubate_fit <- function(object, parm, level = 0.95, R = 199L,
                                 bs_data, bs_infer = c('logquantile', 'lognormal', 'quantile', 'quantile0', 'normal', 'normal0'),
                                 useBoot=FALSE, ...){
  stopifnot(inherits(object, 'incubate_fit'))
  stopifnot(is.numeric(level), length(level) == 1L, level < 1L, level > 0L)
  stopifnot(is.numeric(R), length(R) == 1L, R > 0)
  if (missing(bs_data)) bs_data <- 'parametric'
  if (is.vector(bs_data) && is.character(bs_data)) bs_data <- match.arg(bs_data[[1L]], choices = c('parametric', 'ordinary'))
  bs_infer <- match.arg(bs_infer)
  logTransform <- isTRUE(startsWith(bs_infer, 'log'))

  twoGroup <- isTRUE(object$twoGroup)
  nObs <- if (twoGroup) lengths(object$data) else length(object$data)

  useBoot <- isTRUE(useBoot) || inherits(bs_data, 'boot')

  genBootstrapData <- is.character(bs_data) && length(bs_data == 1L) && ! is.na(bs_data) && nzchar(bs_data)
  stopifnot( genBootstrapData || useBoot && inherits(bs_data, 'boot') || is.matrix(bs_data) )


  # check if we can really use boot
  if ( useBoot &&
       (! requireNamespace("boot", quietly = TRUE) || twoGroup || ! bs_infer %in% c('normal', 'lognormal', 'quantile', 'logquantile', 'quantile0')) ) {
    warning('Using own implementation as package', sQuote('boot'), 'is not available or scenario not implemented.',
            call. = FALSE)
    useBoot <- FALSE
  }

  cf <- coef(object)
  pnames <- names(cf)
  stopifnot( is.numeric(cf), is.character(pnames), nzchar(pnames), length(cf) == length(pnames) )

  if (missing(parm)) parm <- pnames else
    if (is.numeric(parm)) parm <- pnames[parm]
  parm <- intersect(pnames, parm) # in any case

  if (is.null(parm) || ! length(parm) || any(! nzchar(parm))) {
    warning('Invalid parameter name given in argument parm=', call. = FALSE)
    return(invisible(NULL))
  }

  stopifnot( is.character(parm), length(parm) >= 1L )

  a <- (1L - level) / 2L
  a <- c(a, 1L - a)

  # if not already provided get bootstrap data (i.e. coefficients) from fitted model to bootstrapped observations
  if (genBootstrapData) {
    bs_data <- bsDataStep(object = object, bs_data = bs_data, R = R, useBoot = useBoot)
  }
  stopifnot( ! is.vector(bs_data) && ! is.character(bs_data) )
  # set R according to the provided bs_data (in particular important when both R & bs_data object are given)
  R <- if (useBoot) bs_data[['R']] else NCOL(bs_data)
  if (R < 999) warning(glue('Be cautious with the confidence interval(s) because the number of effective bootstrap samples R = {R} is rather low (R<999).'),
                       call. = FALSE)

  # logShift: needed only when log-transformation is requested. Start with a small standard value for all parameters
  logshift <- purrr::set_names(rep_len(.0001, length.out=length(pnames)), nm = pnames)
  # for delay, the transformation needs to be independent of the scale of delay, so we subtract the minimum and add a shift
  #+use fixed logshift_delay = 5 (which performed well in simulation at single group, exponential distribution, together with smd=0.25)
  if (logTransform){
    LOGSHIFT_DELAY <- 5
    for (i in which(startsWith(pnames, 'delay'))){
      logshift[i] <- -min(if (useBoot) bs_data$t[,i] else bs_data[i,], na.rm = TRUE) + LOGSHIFT_DELAY
      # using low quantiles would make it less dependent on R but then we needed to check that x-logshift remains positive (for log)
      #stats::quantile(..i.., probs = c(0, 0.001), na.rm = TRUE, names = FALSE) # catch when diff() > LOGSHIFT_DELAY
    }#rof
  }#fi

  # do bootstrap inference on bootstrap data
  ci <- if (useBoot) {
    stopifnot( inherits(bs_data, 'boot') )

    # 'perc' just takes the quantiles,
    #+'basic' uses quantiles of the difference to the observed value (bias-correction)
    ci_type <- switch(bs_infer,
                      quantile0 = 'perc',
                      quantile =,
                      logquantile = 'basic',
                      normal =,
                      lognormal = 'norm',
                      stop('This boot.ci-type is not supported!'))

    matrix(unlist(
      purrr::map(seq_len(length.out = length(coef(object))), .f = ~ {
        # the output of boot.ci can have different CIs as named matrix list entries
        ci_bo <- {if (logTransform)
          boot::boot.ci(bs_data, index = ., conf = level, type = ci_type,
                        h = function(t) log(t + logshift[[.]]), hdot = function(t) 1/(t + logshift[[.]]),
                        hinv = function(t) exp(t) - logshift[[.]]) else
                          boot::boot.ci(bs_data, index = ., conf = level, type = ci_type)}[[switch(ci_type,
                                                                                                   norm = 'normal',
                                                                                                   perc = 'percent',
                                                                                                   ci_type)]]
        # depending on the CI-type: normal yields 3 columns, perc and others give 5 columns
        stopifnot( is.matrix(ci_bo), NCOL(ci_bo) > 2L )
        # the last two columns are always the lower and upper bound
        ci_bo[, c(NCOL(ci_bo)-1L, NCOL(ci_bo))] })),
      ncol = 2L, byrow = TRUE)
  } else {

    stopifnot( is.matrix(bs_data) )

    # bootstrapped confidence limits
    # bias-correction for parametric bootstrap only!?
    #delayH_mle_bias <- mean(delay_mle_bs) - delayH_mle
    switch(bs_infer,
           quantile0 = {
             t(apply(bs_data, 1L, stats::quantile, probs = a, na.rm = TRUE))
           },
           quantile = {
             # bias-corrected quantile-based CI
             # see Davison, p28
             # vector - matrix: vector is expanded column-wise, and the row-dimension fits (=number of coefs)
             2L * cf - t(apply(bs_data, 1L, stats::quantile, probs = rev(a), na.rm = TRUE))

           },
           logquantile = local({
             # #bs_min <- apply(bs_data, 1L, min) - .15
             # bs_min <- purrr::set_names(rep.int(-.001, length(cf)), nm = names(cf))
             # # for delay, the transformation should be independent of the scale of delay
             # if ('delay' %in% names(bs_min)) bs_min['delay'] <- min(bs_data['delay',], na.rm = TRUE) - .1

             ## bias-corrected normal-based CI after log-transformation
             -logshift + exp(
               2L * log(cf + logshift) - log(t(apply(bs_data, 1L, stats::quantile, probs = rev(a), na.rm = TRUE))+logshift)
             )
           }),
           normal0 = {
             t(c(1L, 1L) %o% .rowMeans(bs_data, m = length(cf), n = R) + stats::qnorm(a) %o% apply(bs_data, 1L, stats::sd))
           },
           normal = {
             ## bias-corrected normal-based CI
             ## ci_delay_mle <- delayH_mle - delayH_mle_bias + c(-1, 1) * qnorm(.975) * delayH_mle_sd
             t(c(1L, 1L) %o% (2L * cf - .rowMeans(bs_data, m = length(cf), n = R)) + stats::qnorm(a) %o% apply(bs_data, 1L, stats::sd))
           },
           lognormal = local({
             # #bs_min <- apply(bs_data, 1L, min) - .15
             # bs_min <- purrr::set_names(rep.int(-.001, length(cf)), nm = names(cf))
             # # for delay, the transformation should be independent of the scale of delay
             # if ('delay' %in% names(bs_min)) bs_min['delay'] <- min(bs_data['delay',], na.rm = TRUE) - .1

             bs_data_h <- log(bs_data + logshift)
             ## bias-corrected normal-based CI after log-transformation
             -logshift + exp(
               t(c(1L, 1L) %o% (2L * log(cf + logshift) - .rowMeans(bs_data_h, m = length(cf), n = R)) + stats::qnorm(a) %o% apply(bs_data_h, 1L, stats::sd)))
           }),
           stop('This type of bootstrap confidence interval is not supported!')
    )
  } #esle useBoot

  # ensure formatted row and column names
  rownames(ci) <- pnames
  colnames(ci) <- paste0(format(a*100, trim = TRUE, nsmall = 1L), '%')

  # enforce parameter bounds also for CI
  # all parameters are non-negative!
  ci[which(ci<0L)] <- 0L


  ci[parm, , drop = FALSE]
}

#' Transform observed data to unit interval
#'
#' The transformation is the probability integral transform. It uses the cumulative distribution function with the estimated parameters of the model fit.
#' All available data in the model fit is transformed.
#'
#' @note
#' This S3-method implementation is quite different from its default method that allows for non-standard evaluation on data frames, primarily intended for interactive use.
#' But the name `transform` fits so nicely to the intended purpose that it is re-used for the probability integral transform, here.
#'
#' @param _data a fitted model object of class `incubate_fit`
#' @param ... currently ignored
#' @return The transformed data, either a vector (for single group) or a list with entries x and y (in two group scenario)
#' @export
transform.incubate_fit <- function(`_data`, ...){
  stopifnot(inherits(`_data`, "incubate_fit"))

  cdfFun <- getDist(`_data`$distribution, type = "cdf")

  twoGroup <- `_data`$twoGroup
  x <- if (twoGroup) `_data`$data$x else `_data`$data

  tr <- purrr::exec(cdfFun, !!! c(list(q=x), coef(`_data`, group = 'x')))
  if (twoGroup) tr <- list(x = tr, y = purrr::exec(cdfFun, !!! c(list(q=`_data`$data$y), coef(`_data`, group = 'y'))))

  tr
}
