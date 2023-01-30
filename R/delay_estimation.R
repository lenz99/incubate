

#' Factory method for objective function, either according to maximum product of spacings estimation ('MPSE')
#' or according to some flavour of maximum likelihood estimation (e.g., naive 'MLEn' or corrected 'MLEc').
#'
#' Given the observed data this factory method produces an objective function
#' which is either the negative of the MPSE-criterion H or the negative log-likelihood for MLE.
#'
#' The objective function takes a vector of model parameters as argument.
#'
#' @details
#' From the observations, negative or infinite values are discarded during pre-processing.
#' In any case, the objective function is to be **minimized**.
#'
#' @param x numeric. observations
#' @param y numeric. observations in second group.
#' @param distribution character(1). delayed distribution family
#' @param twoPhase logical flag. Do we allow for two delay phases where event rate may change? Default is `FALSE`, i.e., a single delay phase.
#' @param bind character. parameter names that are bind together (i.e. equated) between both groups
#' @param method character(1). Specifies the method for which to build the objective function. Default value is `MPSE`. `MLEn` is the naive MLE-method, calculating the likelihood function as the product of density values. `MLEc` is the modified MLE.
#' @param profiled logical. Should scale parameter be profiled out prior to optimization?
#' @param ties character. How to handle ties within data of a group.
#' @param verbose integer flag. How much verbosity in output? The higher the more output. Default value is 0 which is no output.
#' @return the objective function (e.g., the negative MPSE criterion) for given choice of model parameters or `NULL` upon errors
objFunFactory <- function(x, y = NULL,
                          distribution = c("exponential", "weibull"), twoPhase = FALSE, bind = NULL,
                          method = c('MPSE', 'MLEn', 'MLEc', 'MLEw'), profiled = FALSE, ties = c('density', 'equidist', 'random', 'error'),
                          verbose = 0) {

  # setup ----
  stopifnot( is.numeric(x), length(x) > 0, is.null(y) || is.numeric(y) && length(y) > 0 )
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  stopifnot( is.null(bind) || is.character(bind) && length(bind) >= 1 )
  ties <- match.arg(ties)

  stopifnot( is.logical(twoPhase), length(twoPhase) == 1L)
  stopifnot( is.logical(profiled), length(profiled) == 1L)
  # enforce either TRUE or FALSE
  twoPhase <- isTRUE(twoPhase)
  profiled <- isTRUE(profiled)


  # original names: standard names of distribution (say, for a single group)
  oNames <- getDist(distribution, type = "param", twoPhase = twoPhase,
                    twoGroup = FALSE, bind = NULL, transformed = FALSE, profiled = FALSE)



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

      roundOffPrecision <- estimRoundingError(obs, maxObs = 1001L)
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

      if (verbose > 1L && length(obs) <= 50L ){
        cat(glue("New data: {paste(obs, collapse = ', ')}\n"))
      }
    } #fi tiesdiff

    # we have broken all ties
    stopifnot( !any(diff(obs) == 0) )

    obs
  } #fn preprocess

  # overwrite the data vectors with pre-processed data
  if (is.null({x <- preprocess(obs = x)})) return(invisible(NULL))
  y <- preprocess(obs = y)


  # do we have two groups after pre-processing?
  twoGroup <- isTRUE(!is.null(y) && is.numeric(y) && length(y))


  # adjust bind:
  #+enforce the canonical order of dist-parameters and drop unused parameters and empty strings
  #+set to NULL if not effectively two group setting
  bind <- intersect(oNames, bind)

  if (!twoGroup && !is.null(bind) && length(bind)) {
    bind <- NULL
    warning("'bind=' was specified in vain as we have only a single group!", call. = FALSE)
  }#fi


  # adjust profiled:
  #+profiling is only possible if rate1/scale1 is not bound and single phase
  #XXX profiling is not implemented for all cases! we sometimes reverse it to FALSE and just issue a warning
  profiled0 <- profiled
  profiled <- profiled && (! any(c("rate1", "scale1") %in% bind) || length(bind) == length(oNames)) && ! twoPhase
  #&& method %in% c("MLEn", "MLEc", "MLEw") #&& distribution == 'weibull' &&

  if (xor(profiled0, profiled)){
    warning(glue("Option `profiled={profiled0}` was reversed!"), call. = FALSE)
  }
  rm("profiled0")


  # set some coefficient names:
  # coefficient names (now that we have settled the profiling flag)
  # transformed (within optimization function)
  trNames <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = FALSE, bind = NULL,
                     profiled = profiled, transformed = TRUE)
  # full parameter names (for the whole parameter vector spanning all groups)
  # original (including parameters that are profiled out in optimization)
  oNamesFull <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = twoGroup, bind = bind,
                        profiled = FALSE, transformed = FALSE) # profiled = FALSE because we consider here original parameters
  # transformed
  trNamesFull <- getDist(distribution, type = "param", twoPhase = twoPhase, twoGroup = twoGroup, bind = bind,
                         profiled = profiled, transformed = TRUE) # profiled as requested because we consider optimization parameters



  # checks ------------------------------------------------------------------

  # MLEw works only with profiling
  stopifnot( method != 'MLEw' || profiled )

  # check that there is enough data (here we also look at bind= if twoGroup)
  if (!twoGroup && length(x) < length(oNames) ||
      twoGroup && length(x) + length(y) < 2L * length(oNames) - length(bind) && min(length(x), length(y)) < length(oNames) - length(bind)) {
    warning("Too few valid observations provided!", call. = FALSE)
    return(invisible(NULL))
  }


  # parameter handling ----

  # MLEw's W1 as median (for later reference to get scale parameter)
  W1 <- if (method == 'MLEw')
    c(x=(1-1/(9*length(x)))^3, y = if (! is.null(y)) (1-1/(9*length(y)))^3 else 1) else
      c(x=1, y=1)

  stopifnot( ! twoPhase ) #XXX not implemented yet!!


  # provide indices for x and for y
  # where to find the parameters per group in the parameter vector of the objective function
  extractParOptInd <- if (! twoGroup) {
    # single group!
    list(x = seq_along(trNames)) ## Cave: trNames reacts to twoPhase-setting (which I've not thought through, yet)
  } else {
    # two group!
    #XXX exponential && profiled: indices are not correct for two groups, yet!!
    #XXX continue here!! (this would allow to run simul_test.R!)
    if (is.null(bind)) {
      if (distribution == 'exponential') {
        if (profiled) list(x = c(1L), y = c(2L)) else list(x = c(1L, 2L), y = c(3L, 4L))
      } else {
        #weibull:
        if (profiled) list(x = c(1L, 2L), y = c(3L, 4L)) else
          list(x = c(1L, 2L, 3L), y = c(4L, 5L, 6L))
      }
    } else if (length(oNames) == length(bind)) {
      # twoGroup, but all parameters are bound!
      if (distribution == 'exponential') {
        # profiled can actually be true (as each group leads to own scale/rate
        #+but it will be averaged (see mergePars!!)
        if (profiled) {
          #warning("Did not expect `profiled=TRUE` and full bind on all parameters!", call. = FALSE)
          list(x = c(1L), y = c(1L))
        } else list(x = c(1L, 2L), y = c(1L, 2L))
      } else {
        #weibull:
        if (profiled) {
          #warning("Did not expect `profiled=TRUE` and full bind on all parameters!", call. = FALSE)
          list(x = c(1L, 2L), y = c(1L, 2L))
        } else list(x = c(1L, 2L, 3L), y = c(1L, 2L, 3L))
      }
    } else {
      # twoGroups & non-trivial bind

      local({
        # getDist with profiled=TRUE & transformed = FALSE removes the profile parameters (although it is original scale)
        #+ as we need the original parameter names without those of profiling
        oNamesFullProf <- if (!profiled) oNamesFull else
          getDist(distribution = distribution, type = "param", twoPhase = twoPhase, bind = bind,
                  twoGroup = TRUE, profiled = TRUE, transformed = FALSE)

        # locally, drop "rate1/scale1" from oNames when in profiling mode
        # Cave: not robust! Think about (e.g.) twoPhase when profiling! (currently profiling is switched off when twoPhase)
        if (profiled) oNames <- setdiff(oNames, c("rate1", "scale1")) # only *local* temporary change
        nonbind <- setdiff(oNames, bind)
        list(
          x = as.vector(rlang::set_names(charmatch(c(bind, paste0(nonbind, ".x")), oNamesFullProf), nm = c(bind, nonbind))[oNames]),
          y = as.vector(rlang::set_names(charmatch(c(bind, paste0(nonbind, ".y")), oNamesFullProf), nm = c(bind, nonbind))[oNames])
        )
      })
    }
  }

  # provide indices for x and for y
  # where to find the parameters per group in the common parameter vector
  extractParInd <- if (!profiled) {
    # = optimization indices when no profiling
    extractParOptInd
  } else {
    # profiled!
    if (! twoGroup) {
      # profiled single group!
      list(x = seq_along(oNames)) ## Cave: oNames reacts to twoPhase-setting (which I've not thought through, yet)
      #if (distribution == 'exponential') list(x = c(1L, 2L)) else list(x = c(1L, 2L, 3L))
    } else {
      # twoGroup && profiled
      if (is.null(bind)){
        if (distribution == 'exponential') list(x = c(1L, 2L), y = c(3L, 4L)) else
          # weibull
          list(x = c(1L, 2L, 3L), y = c(4L, 5L, 6L))
      } else {
        if ( length(oNames) == length(bind)) {
          # profiled can actually be true (as each group leads to own scale/rate
          #+but it will be averaged (see mergePars!!)
          #warning("Did not expect `profiled=TRUE` and full bind on all parameters!", call. = FALSE)
          if (distribution == 'exponential') list(x = c(1L, 2L), y = c(1L, 2L)) else
            # weibull
            list(x = c(1L, 2L, 3L), y = c(1L, 2L, 3L))
        } else {
          # twoGroup & non-trivial bind
          local({
            # we consider the parameter names on original scale, with profiled parameters also back in!
            nonbind <- setdiff(oNames, bind)
            list(
              x = as.vector(rlang::set_names(charmatch(c(bind, paste0(nonbind, ".x")), oNamesFull), nm = c(bind, nonbind))[oNames]),
              y = as.vector(rlang::set_names(charmatch(c(bind, paste0(nonbind, ".y")), oNamesFull), nm = c(bind, nonbind))[oNames])
            )
          })
        }
      }
    }
  } #esle !profiled

  # parameter transformation matrices
  paramTransf <- list(
    M = switch(distribution,
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
    ),
    Minv = switch(distribution,
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
    ),
    F = list(exponential = c(identity, log, log, log),
             weibull = c(identity, log, #log1p, #identity, #=shape1
                         log, log, log, log))[[distribution]],
    Finv = list(exponential = c(identity, exp, exp, exp),
                weibull = c(identity, exp, #expm1, #identity, #=shape1
                            exp, exp, exp, exp))[[distribution]]
  )

  # transform parameter vector for a single group. Does not use parameter names.
  # transformed parameters are used within optimization. The transformation helps to ensure side-conditions (e.g. log-transformation ensures non-negativitiy of original parameter)
  # @param parV1 parameter vector for a single group
  # @param inverse logical. If `inverse=TRUE` does inverse transformation, from optimization parameters back to original parameters
  # @return transformed parameter vector, unnamed!
  transformPars1 <- function(parV1, inverse = FALSE){

    if (inverse) {
      # b = Ainv %*% Finv(b')
      as.numeric(paramTransf[["Minv"]][seq_along(parV1), seq_along(parV1)] %*%
                   as.numeric(.mapply(FUN = function(f, x) f(x),
                                      dots = list(paramTransf[["Finv"]][seq_along(parV1)], parV1),
                                      MoreArgs = NULL)))
    } else {
      # b' = F(A %*% b)
      as.numeric(.mapply(FUN = function(f, x) f(x),
                         dots = list(paramTransf[["F"]][seq_along(parV1)],
                                     as.numeric(paramTransf[["M"]][seq_along(parV1), seq_along(parV1)] %*% parV1)),
                         MoreArgs = NULL))
    }
  }# fn transformPars1

  # merge two parameter vectors
  # @param isOpt flag: are the parameters on optimization scale?
  # @return merged parameter vector
  mergePars <- function(parx, pary, isOpt){
    exParInd <- if (isOpt) extractParOptInd else extractParInd
    # aggregate parameters (via mean if isOpt or geometric mean if on original scale).
    #+this is necessary for merging start vector for parameter-optimization
    res <- as.vector(tapply(X = c(parx, pary),
                            INDEX = unlist(exParInd),
                            # arithmetic or geometric mean
                            FUN = function(x) {
                              stopifnot(length(x) <= 2L)
                              if (length(x) <= 1L) x else
                                if (isOpt) (x[[1L]] + x[[2L]])/2L else #mean(x)
                                  sqrt(x[[1L]] * x[[2L]]) #prod(x)**(1/length(x))
                              },
                            simplify = TRUE))
    # .. only for delay1 we use minimum as aggregation function (in this case first entry in exParInd$x and exParInd$y is 1!)
    if (exParInd$x[[1L]] + exParInd$y[[1L]] == 2){
      res[[1L]] <- min(parx[[1L]], pary[[1L]])
    }

    res
  }

  # extract parameter vector for a specified group
  # if parameters are for optimization and transformation is requested, unprofiling is done (if relevant)
  # @param group character. Extract parameters for the given group. If NULL, keep all parameters.
  # @param isOpt logical. Are the given parameters on optimization function scale?
  # @param named logical. Extract parameters as named vector?
  # @return parameter vector
  extractPars <- function(parV, group = NULL, isOpt = TRUE, transform = FALSE, named = FALSE){
    # result is on optimization scale?
    resIsOpt <- xor(isOpt, transform)

    # basically, ignore group= when single group: use always canonical "x" then
    if (!twoGroup) {
      group <- "x"
    }

    if (is.null(group)){
      return(local({

        # recursive calls for the individual groups
        parx <- extractPars(parV, group = "x", isOpt = isOpt, transform = transform, named = FALSE)
        pary <- extractPars(parV, group = "y", isOpt = isOpt, transform = transform, named = FALSE)

        # merge the two parameter vectors back together (after a potential transformation)
        res0 <- mergePars(parx = parx, pary = pary, isOpt = resIsOpt)
        if (named){
          res0 <- rlang::set_names(res0, nm = if (resIsOpt) trNamesFull else oNamesFull)
        }
        res0
      }))
    } #fi is.null(group)

    # index vector for specified group
    ind <- if (isOpt) extractParOptInd[[group]] else extractParInd[[group]]

    if (is.null(ind)) return(NULL)

    res <- if (!transform) {
      parV[ind]
    } else {
      # do transform
      local({
        res0 <- transformPars1(parV[ind], inverse = isOpt)

        if (profiled) {
          # un-profile (when going from profiled par_opt to par_orig)
          if (isOpt) {
            # access observations for specified group
            obs <- if (group == "y") y else x
            # add scale parameter at the end of parameter vector
            k <- if (distribution == 'weibull') res0[[2L]] else 1L
            scale0 <- (1/W1[[group]]*mean((obs-res0[[1L]])^k))^(1/k)
            res0 <- c(res0, if (distribution == 'exponential') 1/scale0 else scale0)
          } else {
            # extract only remaining parameters
            res0 <- res0[extractParOptInd[[group]]]
          }
        }
        res0
      })
    }


    # single group names
    if (named) {
      rlang::set_names(res, nm = if (resIsOpt) trNames else oNames)
    } else {
      as.vector(res)
    }
  }


  # optimization arguments -----


  # get optimization start values and upper limits based on observations from a single group
  # for `twoPhase=TRUE` there will be more parameters
  # with profiling no scale parameter is returned (as it is not optimized)
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

                      #parV0 <- rlang::set_names(parV0, nm = oNames)
                      # transform start-parameters for optfun-parametrization
                      parV0 <- transformPars1(parV0, inverse = FALSE)

                      # drop scale if profiling
                      if (profiled) {
                        stopifnot(! twoPhase) #XXX not implemented, yet!
                        parV0 <- parV0[1L] # drop "rate1"
                      }
                      parV0
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
                      # use median rank approximation for empirical Weibull CDF: F(i,n) = (i - 0.3) / (n + 0.4)
                      # and then ordinate is log(1/(1-F)) = -log(1-F) on log-scale
                      start_y <- log(-log(1-stats::ppoints(obs, a=.3)))
                      # cf. lm.fit(x = cbind(1, log(obs)), y = start_y)$coefficients
                      # weighted version with more weight in the middle:
                      # w <- seq_along(obs); w <- w * (max(w)+1-w) #or use plogis-weights to downweight the early obs
                      # lm.wfit(x = cbind(1, log(obs)), y = start_y, w = plogis(-2:(length(obs)-3)))$coefficients
                      start_shape <- stats::cor(log(obs), start_y) * stats::sd(start_y) / stats::sd(log(obs))
                      start_scale <- exp(mean(log(obs)) - mean(start_y) / start_shape) # scale from intercept


                      parV0 <- c( max(DELAY_MIN, obs[[1L]] - 2/(length(obs)+3)),
                                  start_shape,
                                  start_scale )

                      # support 2-phase with additional start parameters
                      if (twoPhase) parV0 <- c(parV0, obs[[floor(.5 + length(obs)/2L)]], parV0[-1L])

                      #parV0 <- rlang::set_names(parV0, nm = oNames)
                      # transform start-parameters for optfun-parametrization
                      parV0 <- transformPars1(parV0, inverse = FALSE)

                      # drop scale if profiling
                      if (profiled) {
                        stopifnot(! twoPhase) #XXX not implemented, yet!
                        parV0 <- parV0[c(1L, 2L)] # drop "scale1"
                      }
                      parV0

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

  # profile likelihood: maximize profiled log-lik f directly
  # if FALSE, go indirectly: consider min(f'^2). As necessary condition, f'^2 == 0
  profiled_llik_directly <- TRUE

  # parameter bounds: set lower & upper bounds
  lowerB <- upperB <- rlang::set_names(rep_len(NA_real_, length(trNamesFull)),
                                       nm = trNamesFull)

  #XXX #QQQ Should this go up to extractPars-function where the transformations are defined???
  PAR_BOUNDS <- list(delay1 = c(lower = 0, upper = NA_real_),
                     delay2 = c(lower = -Inf, upper = NA_real_),
                     rate  = c(lower = -Inf, upper = +Inf),
                     # shape lower bound for MLEnp (actually for shape1)
                     shape = c(lower = if (profiled && method == 'MLEn' && !profiled_llik_directly) 1.49e-8 else -Inf, upper = +Inf),
                     scale = c(lower = -Inf, upper = +Inf))


  # alas, purrr::iwalk did not work for me here
  for (nam in names(PAR_BOUNDS)) {
    idx <- startsWith(trNamesFull, prefix = nam)
    if (any(idx)) {
      lowerB[idx] <- purrr::chuck(PAR_BOUNDS, nam, 'lower')
      upperB[idx] <- purrr::chuck(PAR_BOUNDS, nam, 'upper')
    } #fi
  } #rof



  par0_x <- getParSetting.gr(x)
  parV <-
    if (! twoGroup) {
      # set parameter vector for group 1 and finish upper bound: match delay1 and delay2
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
        } else {

          mergePars(parx = start_x, pary = start_y, isOpt = TRUE)
        }
      } #twoGrp, not all params bound!
    } # twoGrp

  # ensure we have names of transformed parameters
  parV <- rlang::set_names(parV, nm = trNamesFull)

  stopifnot( ! any(is.na(lowerB), is.na(upperB)) )
  # clean up env. # use local() more???
  remove(list=c("PAR_BOUNDS", "par0_x"))


  optim_args <- list(
    par = parV,
    method = "L-BFGS-B",
    lower = lowerB,
    upper = upperB,
    # most parameters are on log-scale.
    control = list(parscale = scalePars(parV, lowerB = 1e-3, upperB = 1e3))
  )




  # objective function ----

  # calculate the log-likelihood, either naive, weighted or in corrected form.
  # What precisely is calculated depends on method but also on the profiled-flag.
  # @param criterion logical. if criterion = TRUE, then pars are on original scale and the proper log-likelihood is returned
  getLogLik <- function(pars, group, criterion = FALSE) {

    # access observations of group
    obs <- if (group == "y") y else x #direct access by name
    #rlang::env_get(env = rlang::env_parent(rlang::current_env(), n=1L), nm = group, inherit = FALSE)

    # extract parameters for specified group on original scale (for CDF)
    pars.gr <- extractPars(pars, group = group, isOpt = !criterion, transform = !criterion)


    # for Weibull, do we want to penalize high shape values in MLE?
    penalize_shape <- distribution == 'weibull' && FALSE #currently turned off (could become an option)

    if (criterion){
      # criterion = log-likelihood
      sum(rlang::exec(getDist(distribution, type = "density"), !!! c(list(x=obs, log=TRUE), pars.gr)))
    } else {
      # calculate the objective function which depends on
      #+ method
      #+ profiled

      n <- length(obs)
      # shape parameter (candidate)
      k <- if (distribution == 'weibull') pars.gr[[2L]] else 1L

      switch(method,
             # MLEc =,
             MLEn = {
               if (profiled && distribution == 'weibull'){
                 obs_c <- obs - pars.gr[[1L]]
                 #cat("\nDelay a: ", pars.gr[["delay1"]], "Shape k: ", k, " (", pars[2], ")\n") #DDD debug

                 # objective function to maximize:
                 if (profiled_llik_directly){
                   # we use 1st derivative to profile out scale parameter but use log-likelihood function directly otherwise
                   # 2nd & 3rd summand could also be: -log(mean(obs_c**k)) + log(k)
                   n * ((k-1) * mean(log(obs_c)) - log(sum(obs_c**k)) + log(n*k) - 1) - penalize_shape*log(k+1)
                 } else {
                   - (1/k + mean(log(obs_c)) - sum(log(obs_c) * obs_c**k) / sum(obs_c**k))**2 -
                     # 1st factor is inverse of harmonic mean
                     (mean(1/obs_c) * sum(obs_c**k)/sum(obs_c**(k-1)) - k/(k-1))**2 -
                     # optional penalization term
                     penalize_shape*log(k+1)
                 }

               } else {
                 # log-likelihood with all parameters
                 sum(rlang::exec(getDist(distribution, type = "density"), !!! c(list(x=obs, log=TRUE), pars.gr)))
               }
             },
             MLEw = {
               stopifnot(profiled)

               obs_c <- obs - pars.gr[[1L]]

               # z is an estimate for the ordered ((x_(i) - a)/gamma)^k ~ Exp(1)
               z <- -log(1-stats::ppoints(n=n, a=.3)) ## = median rank estimates for -log(1-F_i)

               # last term: approximation for -log(GM_n Z) = - AM_n(logZ)
               W2 <- sum(z * log(z)) / (n * W1[[group]]) - log(log(2) - 0.1316 * (1 - 1/n))
               W3 <- if (k==1) {
                 W1[[group]] * mean(1/z)
               } else {
                 W1[[group]] * mean(z^(-1/k)) / mean(z^((k-1)/k))
               }

               if (verbose > 1L){
                 cat(glue("Weights: W2 = {W2} and W3 = {W3} for {group}.",
                          "Candidate values: delay {pars.gr[[1L]]} and shape {k}."), "\n")
               }

               # objective function to maximize
               - (W2/k + mean(log(obs_c)) - sum(log(obs_c) * obs_c**k)/sum(obs_c**k))**2 -
                 # 1st factor is inverse of harmonic mean
                 (mean(1/obs_c) * sum(obs_c**k)/sum(obs_c**(k-1)) - W3)**2 -
                 # optional penalization term
                 penalize_shape*log(k+1)

             },
             MLEc = {
               stopifnot( n >= 2L )
               # objective function to maximize
               log(diff(rlang::exec(getDist(distribution, type = "cdf"), !!! c(list(q=obs[1:2]), pars.gr)))) +
                 sum(rlang::exec(getDist(distribution, type = "density"), !!! c(list(x=obs[-1L], log=TRUE), pars.gr))) -
                 # optional penalization term
                 penalize_shape * log(k+1)
             },

             stop("This method is not handled here!", call. = FALSE))
    }
  }

  # log spacings:
  # calculate the differences in EDF (for given parameters in group) of adjacent observations on log scale
  # The criterion is always the negative mean of these log-spacings.
  # @param pars vector of parameters (by default, on transformed scale, i.e. when criterion = FALSE)
  # @param criterion logical. When `criterion = TRUE`, then pars are on original scale.
  # @return n+1 cumulative diffs on log-scale (or single negative number in twoPhase when delay2 <= delay in quick fix)
  getCumDiffs <- function(pars, group, criterion = FALSE) {

    # access observations of group
    obs <- if (group == "y") y else x #direct access by name
    # get() would (by default) start looking in the execution environment of getCumDiffs and would find the object in its parent
    # or: env = rlang::env_parent() # parent of caller env is (normally!?) the function-env
    # or: env = rlang::fn_env(getCumDiffs) # referring directly to the function environment (but requires function obj)
    #obs <- rlang::env_get(env = rlang::env_parent(rlang::current_env(), n=1L), nm = group, inherit = FALSE)

    # extract parameters for specified group on original scale (for CDF)
    pars.gr <- extractPars(pars, group = group, isOpt = !criterion, transform = !criterion)

    if (verbose > 1L){
      cat(glue("Parameter vector for group {group} after back-transformation: ",
               "{paste(round(pars.gr, 2), collapse = ', ')}"), "\n")
    }

    # calculate spacings
    # contract: data is sorted!
    cumDiffs <- diff(c(0L,
                       rlang::exec(getDist(distribution, type = "cdf"), !!! c(list(q=obs), pars.gr)),
                       1L))


    # use densFun for ties
    # we check difference of obs directly (not cumDiffs)
    #+because cumDiffs can be 0 even if obs are different, in particular for non-suitable parameters!
    ind_t <- which(diff(obs) == 0L)
    if ( length(ind_t) ){
      stopifnot( ties == 'density' ) # other tie-strategies have already dealt with ties in *preprocess*
      # increase index by 1 to get from diff(obs)-indices to cumDiffs-indices
      cumDiffs[1L+ind_t] <- rlang::exec(getDist(distribution, type = "dens"), !!! c(list(x = obs[ind_t]), pars.gr))
    } #fi

    # respect the machine's numerical lower limit
    cumDiffs[which(cumDiffs < .Machine$double.xmin)] <- .Machine$double.xmin

    log(cumDiffs)

  }# fn getCumDiffs


  # Objective function like negative mean log-spacings for MPSE or negative log-likelihood for MSE0
  # Estimate parameters by minimizing this function.
  # param `pars` the vector of parameters. transformed when criterion=FALSE and not transformed when criterion=TRUE
  # param `criterion` a logical flag. If requested, give the original criterion to minimize. Then, the parameters are on original scale
  # param `aggregated` a logical flag. For two group case, `aggregated=FALSE` returns values per group, like mean log cum-diffs per group.
  objFun <- function(pars, criterion = FALSE, aggregated = TRUE) {

    switch(method,
           MPSE = {
             - if (! twoGroup) {
               mean(getCumDiffs(pars, group = "x", criterion = criterion))
             } else {
               #twoGroup:
               #the approach to first merge x and y and then do the cumDiffs, log and mean does *not* work out
               #because the parameters should be optimized within group.
               #merged data lead to frequent non-convergence or visually bad fits
               res <- c(mean(getCumDiffs(pars, group = "x", criterion = criterion)), mean(getCumDiffs(pars, group = "y", criterion = criterion)))

               if (aggregated) stats::weighted.mean(res, w = c(length(x), length(y))) else res
             }
           },
           MLEn = ,
           MLEw = ,
           MLEc = {
             stopifnot( ! twoPhase ) #XXX not implemented yet!

             - if (! twoGroup) getLogLik(pars, group = "x", criterion = criterion) else {
               res <- c(getLogLik(pars, group = "x", criterion = criterion), getLogLik(pars, group = "y", criterion = criterion))

               if (aggregated) sum(res) else res
             }

           },
           stop(glue('Objective function for method {method} is not implemented!'), call. = FALSE)
    )
  } #fn objFun

  # attach analytical solution for MLE
  if ( method == 'MLEn' && ! twoGroup && ! twoPhase && distribution == 'exponential' ){
    attr(objFun, which = "opt") <- local({
      par_analytic <- c(delay1 = x[[1L]], rate1 = 1L/(mean(x) - x[[1L]]))
      list(par_orig = par_analytic,
           #transformed parameters
           par = extractPars(par_analytic, isOpt = FALSE, transform = TRUE, named = TRUE),
           value = length(x) * ( log(mean(x) - x[[1L]]) + 1L ),
           methodOpt = "analytic",
           convergence = 0L,
           message = "analytic solution for naive MLE ('MLEn')",
           counts = 0L)
    })
  }

  objFun
}




#' Fit optimal parameters according to the objective function (either MPSE or MLE-based).
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

  # gather information from objective function environment
  objFunObjs <- rlang::env_get_list(env = objFunEnv,
                                    nms = c("bind", "distribution", "method", "optim_args" ,"trNamesFull", "profiled", "twoGroup", "twoPhase", "x", "y", "extractPars", "oNames"))

  # check if there is already a solution provided by the objective function
  optObj <- attr(objFun, which = "opt", exact = TRUE)

  if ( is.list(optObj) && all( c('par', "par_orig", 'value', 'convergence') %in% names(optObj)) ){
    if (verbose > 0L) cat("Using provided (analytical) solution to objective function.\n")
  } else {
    optObj <- NULL #start from scratch
    # numeric optimization
    if (verbose > 0L) message("Start with numeric optimiziation of objective function.")

    if (is.null(optim_args)) optim_args <- objFunObjs[["optim_args"]]
    stopifnot( is.list(optim_args), 'par' %in% names(optim_args),
               is.numeric(optim_args$par), length(optim_args$par) == length(objFunObjs$trNamesFull) )
    if (is.null(names(optim_args$par))) optim_args$par <- rlang::set_names(optim_args$par, objFunObjs$trNamesFull)
    stopifnot( identical(names(optim_args$par), objFunObjs$trNamesFull))
    # set objective function (overwrite entry 'fn' if it is already present)
    optim_args[["fn"]] <- objFun


    # optim: first attempts ----

    try({
      optObj <- purrr::exec(stats::optim, !!! optim_args)
      optObj$methodOpt <- optim_args$method
    }, silent = TRUE)


    if (is.null(optObj)){
      if (verbose > 0L) warning(glue("{objFunObjs$method}-optimization failed during model fit!"),
                                call. = FALSE)
    } else if ( isTRUE(optObj$convergence > 0L) ){
      # do a 2nd attempt of optim in case it did not converge in the first place
      if (verbose > 1L) message("No proper convergence during 1st optimization in delay fit. Re-try with different parameter scaling.")

      # Use parameter values of non-converged fit as new start values (and adapt parscale accordingly)
      #+The objFun is to be minimized, smaller is better!
      if ( isTRUE(is.numeric(optObj$par) && all(is.finite(optObj$par)) && optObj$value < objFun(optim_args$par)) ){
        optim_args[["par"]] <- optObj$par  # purrr::assign_in(where = "par", value = optObj$par)

        if ( "parscale" %in% names(optim_args[["control"]]) ){
          optim_args[['control']][['parscale']] <- scalePars(optim_args[['par']])
        }

        # optim: 2nd attempt --
        optObj <- NULL
        if (verbose > 1L) message("Do 2nd attempt with renewed start values and parameter scaling")
        try({
          optObj <- purrr::exec(stats::optim, !!! optim_args)
          optObj$methodOpt <- optim_args$method
        }, silent = TRUE)

        if ( is.null(optObj) || isTRUE(optObj$convergence > 0L && verbose > 0L) ) warning("No proper convergence after re-try.", call. = FALSE)
      }## fi rescaling for 2nd attempt
    }## fi 2nd attempt necessary?


    # nlminb (PORT): last attempt ----

    if (is.null(optObj) || optObj$convergence > 0L){
      if (verbose > 1L) message("Do another final attempt with PORT-optimizer.")

      optObj <- minObjFunPORT(objFun = objFun, start = optim_args$par, lower = optim_args$lower, upper = optim_args$upper, verbose = verbose)
    }


    # post-process optObj -----

    #XXX for MLE with indirect profiling (currently MLEw):
    #+check that we have indeed an local maximum for the log-likelihood (as we have only found candidate values by looking for roots of f')

    # set names to parameter vector
    if (! is.null(optObj)){
      stopifnot( 'par' %in% names(optObj) )
      stopifnot( identical(names(optObj$par), objFunObjs$trNamesFull))
      # # set canonical names for parameters
      # optObj$par <- rlang::set_names(optObj$par, objFunObjs$trNamesFull)
      # save optim_args in optimization object (but w/o objective function)
      optim_args$fn <- NULL
      optObj <- append(optObj, values = list(optim_args = optim_args))
    }

    # add par_orig
    optObj <- append(optObj, values = list(par_orig = objFunObjs$extractPars(parV = optObj$par, group = NULL, isOpt = TRUE, transform = TRUE, named = TRUE)))
  } #esle numeric optimization

  optObj
}



#' Fit a delayed Exponential or Weibull model to one or two given sample(s).
#'
#' Maximum product of spacings estimation is used by default to fit the parameters. Estimation via naive maximum likelihood (`method = 'MLEn`) is available, too,
#' but MLEn yields biased estimates. MLEc is a corrected version of MLE due to Cheng.
#'
#' Numerical optimization is done by `stats::optim`.
#' @param x numeric. observations of 1st group. Can also be a list of data from two groups.
#' @param y numeric. observations from 2nd group
#' @param distribution character. Which delayed distribution is assumed? Exponential or Weibull.
#' @param twoPhase logical. Allow for two phases?
#' @param bind character. parameter names that are bind together in 2-group situation.
#' @param ties character. How to handle ties.
#' @param method character. Which method to fit the model? 'MPSE' = maximum product of spacings estimation *or* 'MLEn' = naive maximum likelihood estimation *or* 'MLEw' = weighted MLE' *or* MLEc' = corrected MLE
#' @param profiled logical. Profile out scale from log-likelihood if possible.
#' @param optim_args list. optimization arguments to use. Use `NULL` to use the data-dependent default values.
#' @param verbose integer. level of verboseness. Default 0 is quiet.
#' @return `incubate_fit` the delay-model fit object. Or `NULL` if optimization failed (e.g. too few observations).
#' @export
delay_model <- function(x = stop('Specify observations for at least one group x=!', call. = FALSE), y = NULL,
                        distribution = c("exponential", "weibull"), twoPhase = FALSE,
                        bind = NULL, ties = c('density', 'equidist', 'random', 'error'),
                        method = c('MPSE', 'MLEn', 'MLEw', 'MLEc'), profiled = method == 'MLEw',
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

  method <- if (length(method) == 1L && toupper(method) == 'MSE') {
    message("The method name 'MPSE' is prefered over the previously used name 'MSE'!")
    'MPSE'
  } else method[1L]
  method <- match.arg(method)
  ties <- match.arg(ties)

  if (is.character(bind)){
    if (any(endsWith(bind, suffix = "_tr"))){
      stop("Parameter names to bind= refer to the distribution parameters and not to the transformed parameters of the objective function.", call. = FALSE)
    }

    # translate convenience names (for single phase) to canonical names
    unNmbrdIdx <- !grepl(pattern = "[12]", bind, fixed = FALSE)
    if (any(unNmbrdIdx)){
      bind[unNmbrdIdx] <- paste0(bind[unNmbrdIdx], "1") #interpret un-numbered parameters as referring to phase 1
      if (verbose > 0L) cat("The unnumbered parameter names in bind= are translated to canonical parameter names (=phase 1).\n")
    }
  }#fi bind=


  # objective function ------------------------------------------------------

  objFun <- objFunFactory(x = x, y = y, method = method, profiled = profiled, distribution = distribution,
                          twoPhase = twoPhase, bind = bind, ties = ties, verbose = verbose)
  objFunEnv <- rlang::fn_env(objFun)

  # optimise objective function
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  # return -----
  twoGroup <- rlang::env_get(env = objFunEnv, nm = "twoGroup")
  # overwrite data with  pre-processed data
  x <- rlang::env_get(env = objFunEnv, nm = "x")
  y <- rlang::env_get(env = objFunEnv, nm = "y", default = NULL)

  # /!\ keep in sync with update()!
  structure(
    list(
      data = if (twoGroup) list(x = x, y = y) else x,
      distribution = distribution,
      twoPhase = twoPhase,
      twoGroup = twoGroup,
      method = method,
      bind = rlang::env_get(env = objFunEnv, nm = "bind"),
      ties = ties,
      objFun = objFun,
      par = optObj$par_orig,
      criterion = objFun(pars = optObj$par_orig, criterion = TRUE, aggregated = TRUE),
      optimizer = purrr::compact(c(list(parOpt = optObj$par, valOpt = optObj$value, profiled = profiled),
                                   optObj[c("methodOpt", 'convergence', 'message', 'counts', 'optim_args')]))),
    class = "incubate_fit")
}

#' @export
print.incubate_fit <- function(x, ...){
  coe <- coef(x)
  cat(glue::glue_data(x, .sep = "\n",
                      "Fit a delayed {distribution}{c('', ' with two delay phases')[[1L+twoPhase]]} through{c('', ' profiled')[[1L+optimizer$profiled]]} {switch(method,
                      MPSE = 'Maximum Product of Spacings Estimation (MPSE)', MLEn = 'naive Maximum Likelihood Estimation (MLEn)',
                      MLEw = 'weighted Maximum Likelihood Estimation (MLEw)',
                      MLEc = 'corrected Maximum Likelihood Estimation (MLEc)', '???')} for {c('a single group', 'two independent groups')[[1L+twoGroup]]}.",
                      "Data: {if (twoGroup) paste(lengths(data), collapse = ' and ') else length(data)} observations, ranging from {paste(signif(range(data), 4), collapse = ' to ')}",
                      "Criterion: {signif(criterion,3)}",
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
  stopifnot( inherits(object, "incubate_fit") )
  transformed <- isTRUE(transformed)

  rlang::env_get(rlang::fn_env(object$objFun), nm = "extractPars")(purrr::chuck(object, !!! if (transformed) list("optimizer", "parOpt") else "par"),
                                                                   group = group, isOpt = transformed, transform = FALSE, named = TRUE)
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
update.incubate_fit <- function(object, optim_args = NULL, verbose = 0, ...){

  stopifnot( all(c("data", "distribution", "method", "objFun", "twoPhase", "twoGroup", "par", "criterion", "optimizer") %in% names(object)) )

  ## fit model with given optim_args
  objFun <- object[["objFun"]]
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))


  x <- if (object[["twoGroup"]]) object[["data"]]["x"] else object[["data"]]
  y <- if (object[["twoGroup"]]) object[["data"]]["y"] else NULL

  # update all relevant fields in the list
  # /!\ keep in sync with delay_model() /!\
  object[c("par", "criterion", "optimizer")] <- list(par = optObj$par_orig,
                                                     criterion = objFun(pars = optObj$par_orig, criterion = TRUE, aggregated = TRUE),
                                                     # drop NULLs from list (e.g. if optim_args is not present)
                                                     optimizer = purrr::compact(c(
                                                       list(parOpt = optObj$par, valOpt = optObj$value, profiled = object$optimizer$profiled),
                                                       optObj[c("methodOpt", "convergence", "message", "counts", "optim_args")])))

  object
}

#' @export
plot.incubate_fit <- function(x, y, title, subtitle, ...){
  # parameter y comes from the plot-generic. y is not used here.
  stopifnot( inherits(x, "incubate_fit") )

  rlang::check_installed(pkg = 'ggplot2', reason = 'to draw plots', version = '3.3')

  cumFun <- getDist(x[["distribution"]], type = "cdf")

  p <- grNames <- NULL

  # catch the one-group case!
  if ( x[["twoGroup"]] ){
    stopifnot( is.list(x[["data"]]) )

    grNames <- names(x[["data"]])

    p <- ggplot2::ggplot(data = tibble::enframe(unlist(x[["data"]]), name = "group"),
                         mapping = ggplot2::aes(x = value, col = substr(x = group, 1L, 1L))) +
      # add estimated delay models
      ggplot2::geom_function(mapping = ggplot2::aes(col = "x"), inherit.aes = FALSE,
                             fun = cumFun, args = coef(x, group = "x"), linetype = "dashed") +
      ggplot2::geom_function(mapping = ggplot2::aes(col = "y"), inherit.aes = FALSE,
                             fun = cumFun, args = coef(x, group = "y"), linetype = "dashed")

  } else {
    grNames <- "x"
    p <- ggplot2::ggplot(data = tibble::tibble(value=x[["data"]]),
                         mapping = ggplot2::aes(x = value)) +
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
    ggplot2::coord_trans(y = "reverse") + # transforms "after_stat" which matters for stat_ecdf
    ggplot2::labs(x = 'Time', y = 'Cumulative prop. of events',
                  col = 'Group',
                  title = title, subtitle = subtitle)
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
                              # objective function with delay1-entries a little bit varied
                              # use `criterion = TRUE` to operate directly on the original parameters
                              # del1_ind: delay1-index within group
                              .f = ~ object$objFun(pars = replace(coefVect, del1_ind, .x), criterion = TRUE, aggregated = FALSE)[[groupIdx]])
    )

    # we like to drop last entry (delay = 1st observation) as objective function tends to explode
    # but we have to keep last entry if it corresponds to the delay estimate (e.g., as is the case for MLEn-fitting)
    if (delayCandDF$delay1[NROW(delayCandDF)] > del_coef) delayCandDF <- delayCandDF[-NROW(delayCandDF),, drop = FALSE]
    # relative change to optimal value, will be negative as objective function is minimized
    delayCandDF$objValInv <- (object$criterion - delayCandDF$objVal) / (object$criterion+.01)
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
  logshift <- rlang::set_names(rep_len(.0001, length.out=length(pnames)), nm = pnames)
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
             # bs_min <- rlang::set_names(rep.int(-.001, length(cf)), nm = names(cf))
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
             # bs_min <- rlang::set_names(rep.int(-.001, length(cf)), nm = names(cf))
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
