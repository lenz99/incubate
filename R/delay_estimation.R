

#' Factory method for objective function
#'
#' Given the observed data this factory method produces an objective function
#' which is either the negative of the MPSE-criterion H or some flavour of the negative log-likelihood for MLE.
#' Implemented variants of MLE-objective functions are naive MLE (`'MLEn'`), corrected MLE (`'MLEc'`) or weighted MLE (`'MLEw'`).
#' In any case, the objective function is to be **minimized**.
#'
#' @details
#' The objective function takes a vector of model parameters as argument.
#' From the observations, negative or infinite values are discarded during pre-processing.
#'
#' @param x numeric. observations
#' @param y numeric. observations in second group.
#' @param distO distribution object
#' @param twoPhase logical flag. Do we allow for two delay phases where event rate may change? Default is `FALSE`, i.e., a single delay phase.
#' @param bind character. parameter names that are bind together (i.e. equated) between both groups
#' @param method character(1). Specifies the method for which to build the objective function. Default value is `MPSE`. `MLEn` is the naive MLE-method, calculating the likelihood function as the product of density values. `MLEc` is the modified MLE.
#' @param profiled logical. Should scale parameter be profiled out prior to optimization?
#' @param ties character. How to handle ties within data of a group.
#' @param verbose integer flag. How much verbosity in output? The higher the more output. Default value is 0 which is no output.
#' @return the objective function (e.g., the negative MPSE criterion) for given choice of model parameters or `NULL` upon errors
objFunFactory <- function(x, y = NULL, distO,
                          twoPhase = FALSE, bind = NULL,
                          method = c('MPSE', 'MLEn', 'MLEc', 'MLEw'), profiled = FALSE, ties = "density",
                          verbose = 0) {

  # setup ----
  stopifnot(is.numeric(x), length(x) > 0, is.null(y) || is.numeric(y) && length(y) > 0)
  method <- match.arg(method)
  stopifnot(!missing(distO), is.list(distO))
  stopifnot(is.null(bind) || is.character(bind) && length(bind) >= 1)


  stopifnot(is.logical(twoPhase), length(twoPhase) == 1L)
  stopifnot(is.logical(profiled), length(profiled) == 1L)
  # enforce either TRUE or FALSE
  twoPhase <- isTRUE(twoPhase) && distO$twoPhaseAllowed
  profiled <- isTRUE(profiled)


  # original names: standard names of distribution (say, for a single group)
  oNames <- distO$param(twoPhase = twoPhase, twoGroup = FALSE, bind = NULL, transformed = FALSE, profiled = FALSE)



  # data preparation ----

  # unify Surv-type but keep numeric if no censoring
  respL <- prepResponseVar(x0 = x, y0 = y, simplify = TRUE)
  stopifnot(is.list(respL), identical(names(respL), c("x", "y")))
  x <- respL[["x"]]
  y <- respL[["y"]]
  rm(list = "respL")

  # flag if we have Surv-data or not
  isSurv <- inherits(x, what = "Surv")

  # Data preprocessing per group:
  # negative, NA and infinite values are dropped. Data gets sorted.
  # @param obs: data vector of one group
  # @return sorted, cleaned up data vector or NULL in case of trouble
  preprocessF <- function(obs) {

    if (is.null(obs) || ! is.numeric(obs)) return(NULL)


    if (!isSurv) {
      # numeric response, non-Surv

      # fix numeric instabilities to have proper ties (when observations are pretty close)
      obs <- survival::aeqSurv(Surv(obs), tolerance = TOL_NUM)[, 1L, drop = TRUE]

      ind_neg <- which(obs < 0L)
      if (length(ind_neg) && !distO$negAllowed) {
        warning("Negative values in data", deparse(substitute(obs)), "! These are dropped.", call. = FALSE)
        obs <- obs[-ind_neg]
      }# fi
      # drop NA and +/-Inf & sort #XXX sort.int?
      obs <- sort(obs[is.finite(obs)])

      if (!length(obs)) {
        warning("Insufficient data! Only ", if (!distO$negAllowed) "non-negative and ", "finite real values are valid.", call. = FALSE)
        return(invisible(NULL))
      }# fi

      # check spread in data
      if (obs[[length(obs)]] < obs[[1L]] + 3L*TOL_NUM) { # && method %in% c("MPSE", "MLEc")) {
        warning("Too small spread in data for this estimation method!", call. = FALSE)
        return(invisible(NULL))
      }

    } else {
      # Surv response
      survType <- attr(obs, which = "type", exact = TRUE)
      # for MPSE: check we only have right-censoring
      if (method == 'MPSE' && survType != 'right') {
        warning("MPSE-fitting supports only right censored observations currently.", call. = FALSE)
        return(invisible(NULL))
      }

      # fix numeric instabilities to have proper ties
      obs <- survival::aeqSurv(obs, tolerance = TOL_NUM)

      # drop negative times
      ind_neg <- which(obs[,1L] < 0L)
      if (length(ind_neg) && !distO$negAllowed) {
        warning("Negative values in data", deparse(substitute(obs)), "! These are dropped.", call. = FALSE)
        obs <- obs[-ind_neg, , drop=FALSE]
      }
      # check finite for right-, left-censored or interval-censored Surv-times
      if (survType %in% c("right", "left", "interval")) {
        obs <- obs[which(is.finite(obs[,1L])), , drop=FALSE]
      }
      # sort by time (first column)
      obs <- sort(obs)

      if (!length(obs) || length(which(obs[, "status"] == 1L)) < 1L + method %in% c("MPSE", "MLEc") ) {
        warning(glue("Insufficient data! Only ", if (!distO$negAllowed) "non-negative and ", "finite real values are valid, ",
                     "at least {c('one observed event time', 'two observed event times')[[1L + method %in% c('MPSE', 'MLEc')]]}",
                     "required for estimation method {method}."), call. = FALSE)
        return(invisible(NULL))
      }# fi

      # check spread in data (any, observed or censored times)
      if (obs[length(obs), 1L] < obs[1L, 1L] + 3L*TOL_NUM ) { #&& method %in% c("MPSE", "MLEc")) {
        warning("Too small spread in data for this estimation method!", call. = FALSE)
        return(invisible(NULL))
      }
    } #esle isSurv

    obs
  } #fn preprocessF

  # overwrite the data vectors with pre-processed data
  if (is.null({x <- preprocessF(obs = x)})) return(invisible(NULL))
  y <- preprocessF(obs = y)


  # do we have two groups after pre-processing?
  twoGroup <- isTRUE(!is.null(y) && is.numeric(y) && length(y))


  # tie-informations for a group
  # @param obs data from a single group.
  # @return list with tie indices or NULL if no data is given
  tieInformationF <- function(obs) {
    ##obs <- graphite[12:20][-4] ##test case: twins + triplicate
    if (is.null(obs)) return(invisible(NULL))

    roundOffPrecision <- estimRoundingError(obs, n_obs = 1001L)
    if (verbose > 0L) {
      cat(glue("Round-off error has magnitude {roundOffPrecision}."), "\n")
    }

    if (!length(obs)) return(invisible(NULL))

    # for Surv, we only consider duplicated observed event times, here.
    # as we only need to fix ties within observed times, ties in censored times use interpolation!?
    dupInd <- if (isSurv) {
      which(duplicated(obs) & obs[, "status"] == 1)
    } else {
      which(duplicated(obs))
    }

    # # reduce obs to only the observed event times for the remainder of the function
    # obs <- obs[which(obs[, "status"] == 1), 1L]

    tieGrp <- matrix(NA_real_, nrow=0, ncol = 0)
    # index vector for cumDiff where spacing will be zer0 due to ties
    cumDiffInd <- integer(0L)

    # rounding radius:
    # used to break ties later on when evaluating MPSE-criterion
    # it can't be wider than smallest observed diff.
    # plogis to mitigate the effect of sample size: the larger the sample the more we can 'trust' the observed minimal diff
    diffObs <- if (isSurv) {
      outInd <- union(dupInd, which(obs[, "status"] != 1))
      diff(obs[if (length(outInd)) -outInd else TRUE, 1L])
    } else {
      diff(obs[if (length(dupInd)) -dupInd else TRUE]) #use unique???
    }

    rRad <- TOL_NUM + .5 * min(roundOffPrecision,
                               # obs[1L] = min(obs) = diff of minimal obs with 0
                               abs(if (isSurv) obs[which.max(obs[, "status"] == 1), 1L] else obs[[1L]]), # very first time obs[[1L]] should be non-negative, anyways.
                               #stats::plogis(q = .1+length(diffObs), scale = 17) * diffObs,
                               diffObs,
                               na.rm = TRUE)

    if (length(dupInd)) {
      stopifnot(dupInd[[1]] > 1L) # duplicated entries start at least 2
      if (ties == "error") stop("Ties within data are not allowed (ties == 'error')!", call. = FALSE)

      gapsInDupInds <- c(1L, which(diff(dupInd)>1)+1) # +1 to be on dupInd-scale
      nbrTieGroups <- length(gapsInDupInds)
      if (verbose > 0L) {
        cat(glue('{nbrTieGroups + length(dupInd)} tied observations ',
                 'in {nbrTieGroups} group(s) within data vector.\n'))
      }
      # start one position before duplicated-indices
      startInd <- dupInd[gapsInDupInds] - 1L
      # tabulate instead of for-loop
      len <- tabulate(bin = if (isSurv) factor(obs[c(startInd, dupInd), 1L], levels = obs[startInd, 1L]) else
        factor(obs[c(startInd, dupInd)], levels = obs[startInd]),
        nbins = length(startInd))
      # len <- rep_len(-1L, length.out = length(startInd))
      # for (i in seq_along(len)) {
      #   j <- 1
      #   while (startInd[[i]]+j <= length(obs) && ! obs[[startInd[[i]]+j]] > obs[[startInd[i]]]) { j <- j+1 }
      #   len[[i]] <- j
      # } #rof
      stopifnot(all(len >= 2)) # at least 2 observations per tie-group
      tieGrp <- cbind(startInd, len) # 2-column matrix
      # index vector for cumDiff where spacing will be zer0 due to ties
      # unlist could become purrr::list_c (req v1.0.0)
      cumDiffInd <- rep.int(startInd, times = len-1L) +
        unlist(purrr::map(.x = len-1, .f = function(.x) seq_len(length.out = .x)))
    }#fi dupInd

    # return tie information (per group)
    list(tieGrp = tieGrp,
         cumDiffInd = cumDiffInd,
         numPrecision = c(precEstim = roundOffPrecision, rRad = rRad))
  } #tieInformationF

  # store tie information per group
  tieInfo <- purrr::compact(list(strategy = ties,
                                 x = tieInformationF(obs = x),
                                 y = tieInformationF(obs = y)))


  # adjust bind:
  #+enforce the canonical order of dist-parameters and drop unused parameters and empty strings
  #+set to NULL if not effectively two group setting
  bind <- intersect(oNames, bind)

  if (!twoGroup && !is.null(bind) && length(bind)) {
    bind <- NULL
    warning("'bind=' was specified in vain as we have only a single group!", call. = FALSE)
  }#fi


  # adjust profiled:
  #+profiling is not implemented for some cases! Then, we reverse it to FALSE and just issue a warning
  #+profiling is only possible if rate1/scale1 is not bound and single phase
  profiled0 <- profiled
  profiled <- profiled && (! any(c("rate1", "scale1") %in% bind) || length(bind) == length(oNames)) && ! twoPhase
  #&& method %in% c("MLEn", "MLEc", "MLEw") #&& distribution == 'weibull' &&

  if (xor(profiled0, profiled)) {
    warning(glue("Option `profiled={profiled0}` was reversed to profiled={profiled}!"), call. = FALSE)
  }
  rm("profiled0")

  # KM fit (used for plotting later)
  survDat <- tibble(time = if (isSurv) c(x, y) else Surv(c(x,y)),
                    groupVar = rep.int(c("x", "y"), times = c(length(x), length(y))))
  kmFit <- survival::survfit(time ~ groupVar, data = survDat,
                             start.time = 0, se.fit = FALSE, conf.type = "none")

  cens <- local({
    # little helper function to count the censored observed by type (right, left, interval)
    censDescF <- function(.x, what = c("n", "ind", "rcens")) {
      what <- match.arg(what)

      # mock survfit-object
      rcensDummy <- list(surv=rlang::rep_along(kmFit$surv, 1))

      if (!isSurv) {
        return(list(n = c(right = 0L, left = 0L, interval = 0L, any = 0L),
                    ind = list(right = integer(0), left = integer(0), interval = integer(0), obs = seq_along(.x)),
                    rcens = rcensDummy)[[what]])
      }

      switch(what,
             n = {
               nvctr <- switch(attr(.x, which = "type", exact = TRUE),
                               right = c(sum(.x[, "status"] == 0), 0, 0),
                               left = c(0, sum(.x[, "status"] == 0), 0),
                               interval = tabulate(.x[, "status"]+1L, nbins = 4L)[-2L],
                               stop("This type of censoring is not supported!", call. = FALSE))

               rlang::set_names(append(nvctr, sum(nvctr)),
                                nm = c("right", "left", "interval", "any"))
             },
             ind = {
               switch(attr(.x, which = "type", exact = TRUE),
                      right = list(right = which(.x[, "status"] == 0), left = integer(0L), interval = integer(0L), obs = which(.x[, "status"] == 1)),
                      left = list(right = integer(0L), left = which(.x[, "status"] == 0), interval = integer(0L), obs = which(.x[, "status"] == 1)),
                      interval = list(right = which(.x[, "status"] == 0), left = which(.x[, "status"] == 2), interval = which(.x[, "status"] == 3),
                                      obs = which(.x[, "status"] == 1)),
                      stop("This type of censoring is not supported!", call. = FALSE) )
             },
             rcens = {
               stopifnot(is.numeric(.x), length(.x) == 1L, .x >= 0L) #.x is nbr of right censorings in the data
               if (.x > 0L) {
                 # treat right-censorings as events and rest as censoring
                 survival::survfit(Surv(time[, 1L], event = !time[, "status"], type = "right") ~ groupVar,
                                   data = survDat, conf.type = "none", se.fit = FALSE)
               } else {
                 rcensDummy
               }
             },
             stop("This request ", sQuote(what, q = FALSE), " is not supported here!", call. = FALSE)
      )
    } #fn censDescF

    retL <- list(isSurv = isSurv,
                 n = purrr::compact(list(x=censDescF(x, what = "n"),
                                         y = if (twoGroup) censDescF(y, what = "n"))),
                 ind = purrr::compact(list(x=censDescF(x, what = "ind"),
                                           y = if (twoGroup) censDescF(y, what = "ind"))))
    # add KM estimator for right-censorings (combines all data, even when two groups, in one object)
    retL[["rcens"]] <- censDescF(retL$n$x[["right"]] + if (twoGroup) retL$n$y[["right"]] else 0,
                                 what = "rcens")

    retL
  })

  # kmFit and rcens have same number of rows
  stopifnot(length(kmFit$surv) == length(cens$rcens$surv))

  # indices of first two relevant observations (for MLEc)
  indForefront <- if (method != "MLEc") NULL else local({
    # little helper function to get the indices for the first two smallest observed values (non-censorings)
    #+in a sorted vector of observations
    forefrontIndF <- function(group) {
      stopifnot(!missing(group), is.character(group), length(group) == 1L)
      obs <- if (group == "y") y else x

      ind_obs1 <- ind_next <- integer()

      if (isSurv) {
        # Surv-response
        cind_gr <- cens$ind[[group]]
        cindo_gr <- cens$ind[[group]]$obs
        stopifnot(length(cindo_gr) >= 2L)

        # check for easy case: no tie at first two observed event times
        if (obs[cindo_gr[2L], 1L] > obs[cindo_gr[1L], 1L] + TOL_NUM) {
          ind_obs1 <- cindo_gr[1L]
          ind_next <- cindo_gr[2L]
        } else {
          # walk down the observed event times
          i1 <- 2L
          while (obs[cindo_gr[i1], 1L] == obs[cindo_gr[1L], 1L]) {
            i1 <- i1 + 1L
          }
          ind_obs1 <- cindo_gr[seq_len(i1-1L)]
          if (length(cindo_gr) >= i1) ind_next <- cindo_gr[i1]
        }

      } else {
        # numeric response, non-Surv
        stopifnot(length(obs) >= 2L)
        # check for easy case: no tie at beginning
        if (obs[[2L]] > obs[[1L]] + TOL_NUM) {
          ind_obs1 <- 1L
          ind_next <- 2L
        } else {
          # get indices for 1st and 2nd observation. Try with few first observations first (for better performance)
          for (l in sort.int(unique(c(5, 10, 50, 100, 500, 1000, length(obs))))) {
            if (l > length(obs)) break
            obs_r <- rank(obs[seq_len(l)], ties.method = "min", na.last = TRUE)
            #which.max(obs_r > 1) # 1st index of 2nd obs
            firstTwoRanks <- unique(obs_r)[c(1L, 2L)]
            # check that there are two distinct values
            if (any(is.na(firstTwoRanks))) {
              if (l < length(obs)) next
              if (method %in% c("MPSE", "MLEc")) stop("At least two different distinct observation values per group required!", call. = FALSE)
              #else warning("Only a single unique distinct observation value in a group.", call. = FALSE)
            } #fi

            ind_obs1 <- which(obs_r == firstTwoRanks[1L])
            ind_next <- which(obs_r == firstTwoRanks[2L])
            if (length(ind_next)) ind_next <- ind_next[1L]
          } #rof
        } #esle
      } #esle (non-Surv)

      list(inds_obs1 = ind_obs1, ind_next = ind_next)
    }

    purrr::compact(list(x = forefrontIndF(group = "x"), y = if (twoGroup) forefrontIndF(group = "y")))
  })

  # set some coefficient names:
  # coefficient names (now that we have settled the profiling flag)
  # transformed (within optimization function)
  trNames <- distO$param(twoPhase = twoPhase, twoGroup = FALSE, bind = NULL, profiled = profiled, transformed = TRUE)
  # full parameter names (for the whole parameter vector spanning all groups)
  # original (including parameters that are profiled out in optimization)
  oNamesFull <- distO$param(twoPhase = twoPhase, twoGroup = twoGroup, bind = bind,
                            profiled = FALSE, transformed = FALSE) # profiled = FALSE because we consider here original parameters
  # transformed
  trNamesFull <- distO$param(twoPhase = twoPhase, twoGroup = twoGroup, bind = bind,
                             profiled = profiled, transformed = TRUE) # profiled as requested because we consider optimization parameters



  # checks ------------------------------------------------------------------

  # MLEw works only with profiling
  stopifnot(method != 'MLEw' || profiled)

  # check that there is enough data (here we also look at bind= if twoGroup)
  if (!twoGroup && length(x) < length(oNames) ||
      twoGroup && length(x) + length(y) < 2L * length(oNames) - length(bind) && min(length(x), length(y)) < length(oNames) - length(bind)) {
    warning("Too few valid observations provided!", call. = FALSE)
    return(invisible(NULL))
  }


  # parameter handling ----

  weights <- if (method != "MLEw") list(W1 = c(x=1, y=1)) else {
    local({

      # Little helper to get the so-called z-values z_i := -log(1-F_i) = log(1/(1-F_i))
      #+which define the weights W1-W3
      # Median rank is a general way to estimate F_i (using a binomial model)
      # Benard's approximation estimates F_i as (i - a) / (N + 1 - 2*i) for some a
      #+a=.3 is recommended by Fothergill (1990) ***
      #+a=.3175 due to Filliben, "The probability plot.." (1975)
      # Assuming 3-parameter Weibull holds, we have z_i = ((x_(i) - a)/gamma)^k ~ Exp(1).
      # Exact values for 1st and last (=nth) entry are known (see "A reliable algorithm..", Jacquelin, 1993)
      # Cousineau uses MC-simulation, drawing from Exp(1). He does not use the observed data to derive F_i.
      # He chooses weights W1-W3 as median of their sampling distribution in MC irrespective of the concrete sample.
      # @param group Specifies for which group to estimate the z's
      # @param method How to estimate the z's. mr = median rank method to estimate F_i. mr_exact for short data sample
      # @param propagateTies logical. Should ties in the observations lead to ties in the z's as well?
      # @return numeric vector of ordered z's, same length as number of observed event time values in group
      zF <- function(group = "x", method = c("mr_exact", "mr_benard"), a = 0.3, propagateTies = FALSE) {
        method <- match.arg(method)
        nObs <- if (group == "y") length(y) else length(x)

        if (!isSurv) {
          # numeric response, non-Surv

          if (propagateTies) {
            nObs0 <- nObs # save original length just to double-check
            obs <- if (group == "y") y else x
            ind_doz <- which(diff(obs) == 0)
            nObs <- nObs - length(ind_doz)
          }

          z0 <- if (nObs < 2) {
            .5
          } else if (method == "mr_exact" && nObs < 89L) { #use exact median rank values if not too many observations
            stats::qbeta(p=.5, shape1 = seq_len(nObs), shape2 = rev(seq_len(nObs)))
          } else {
            # Benard-style approximation for long observation vectors
            # 1st and last entry are still exact median rank values
            z_n <- .5^(1/nObs)
            c(1-z_n, stats::ppoints(n = nObs, a = a)[1L+seq_len(nObs-2L)], z_n)
          }

          if (propagateTies && length(ind_doz)) {
            ind_rept <- rep_len(1L, length.out = length(z0))
            # index to update ind_rept
            iupd <- ind_doz[1L]
            idoz <- 1L

            while (idoz <= length(ind_doz)) {
              tie_cnt <- 1
              # count ties in group
              while(idoz + tie_cnt <= length(ind_doz) && ind_doz[idoz+tie_cnt] == ind_doz[idoz + tie_cnt-1] + 1) {
                tie_cnt <- tie_cnt + 1
              } #elihw

              # update tie count (for rep-times)
              ind_rept[iupd] <- ind_rept[iupd] + tie_cnt
              # update indices
              idoz <- idoz + tie_cnt - 1 # to end of tie group
              if (idoz < length(ind_doz)) {
                iupd <- iupd + ind_doz[idoz+1]-ind_doz[idoz]-1 # move iupd for next update
              }
              # move on
              idoz <- idoz + 1
            }# elihw

            z0 <- rep.int(z0, times = ind_rept)
            stopifnot(length(z0) == nObs0)
          } #fi

          return(-log(1-z0))
        }#fi !isSurv

        stopifnot(isSurv)
        switch(EXPR = attr(x, which = "type", exact = TRUE),
               right = {
                 # nbr of events observed
                 n_ev <- nObs - cens$n[[group]][["right"]]

                 # Cousineau estimates the weights from Monte-Carlo simulation study (MCSS) irrespective of the concrete sample
                 # But here we have censorings which also effect F_i and we hence chose to estimate F_i from the concrete sample with its censoring scheme
                 # Estimate F_i via Kaplan-Meier (copes with censorings) with Benard-style median rank estimation to avoid 0 and 1
                 # z is an estimate for the ordered z_i = -log(1-F_i) = ((x_(i) - a)/gamma)^k ~ Exp(1)
                 # unique event times (in all available groups)
                 ind_evKM <- which(kmFit$n.event > 0.99) #at least one event (type=left/interval makes that we get fractional numbers here [but 0 is 0 also for interval!?])
                 # get the right subset of indices for specified group (when having two groups)
                 if (twoGroup) {
                   # Cave: works only for two groups (x or y) as I only use the strata[[1L]] as cutpoint
                   ind_evKM <- if (group == "x") ind_evKM[ind_evKM <= kmFit$strata[[1L]]] else ind_evKM[ind_evKM > kmFit$strata[[1L]]]
                 }
                 # n.event is generally not integer for type=interval/left. It is increased by a fraction (depending on number of events) and sums to nbr of events+1 (per group)
                 # floor(n.event + n.censor) = n
                 stopifnot(sum(as.integer(kmFit$n.event[ind_evKM]),
                               if (twoGroup) kmFit$n.censor[(if (group == "x") 1 else -1) * seq_len(kmFit$strata[[1L]])] else kmFit$n.censor) == kmFit$n[[if (group == "x") 1L else 2L]])

                 # estimated survival probabilities for event times, replicated
                 # n.event is not always integer for Surv-type=interval/left. rep.int truncates floats & it should always work.
                 kmSurvProb <- rep.int(kmFit$surv[ind_evKM], times = kmFit$n.event[ind_evKM])
                 stopifnot(length(kmSurvProb) == n_ev)

                 # Benard-style median-rank estimation (avoid 0 and 1)
                 -log(1-((1-kmSurvProb) * n_ev - a) / (n_ev + 1 - 2*a))
               },
               stop("This Surv-type is not handled here!", call. = FALSE)
        )
      } #fn zF

      z_x <- zF(group = "x", propagateTies = TRUE)
      z_y <- if (twoGroup) zF(group = "y", propagateTies = TRUE)

      # how to calculate the weights W1-W3?
      method_w1 <- if (isSurv) "sample" else "sdist_median"
      method_w2 <- if (isSurv) "sample" else "sdist_median"
      method_w3 <- if (isSurv) "sample" else "sdist_median"

      # little helper function to calculate W1-weight (as function of n)
      # W1 = mean(z_i) follows a gamma-dist with parameters shape=n and scale=1/n and we estimate W1 as median of it.
      # W1 is also used to get scale parameter during un-profiling.
      # We count all events because it is used to get scale parameter (and in this formula we already correct for censorings),
      #+e.g., nObs = length(x), even when there is cens$n$x[["any"]]
      # @param method By which method to calculate weights W1? 'sample' will use the mean of the provided sample of z-values, sdist_median uses the median of the sampling distribution (MC-sim)
      w1F <- function(nObs, z, method = c("sample", "sdist_median", "hybrid")) {
        method <- match.arg(method)

        #Using a MC-simulation
        #+cf. Cousineau's simulation results for median of W1's sampling distribution
        #+"Nearly unbiased estimators.." (2009), Table 2, column J_1

        switch(EXPR = method,
               sample = {
                 if (missing(z) || !is.numeric(z) || length(z) == 0L) {
                   stop("Please provide the vector of z's to estimate W1!", call. = FALSE)
                 }
                 mean(z)
               },
               sdist_median = {
                 # use W1-function
                 .MLEw_approx$fun$w1F(nObs)
               },
               hybrid = {
                 (.MLEw_approx$fun$w1F(nObs) + mean(z)) / 2L
               },
               stop("This method for estimating W1 is not handled here!", call. = FALSE)
        )
      } #w1F

      # W1 weights: use full length (even when censored obs are present?!)
      W1_x <- w1F(nObs = length(x), z = z_x, method = method_w1) #or # nObs = length(x) - cens$n$x[["any"]]),
      W1_y <- if (twoGroup) w1F(nObs = length(y), z = z_y, method = method_w1) else 1 #length(y) - cens$n$y[["any"]])


      w2F <- function(z, method = c("sample", "sdist_median", "hybrid")) {
        method <- match.arg(method)
        nz <- length(z)

        # MC-simulation on W2 for n=1..16
        # cf. Cousineau's simulation results for median of W2's sampling distribution
        #+"Nearly unbiased estimators.." (2009), Table 3, column J_2
        # W2-approximation via asymptotic regression model SSasymp on log(n):
        # We hence model: W2 = 1 + (R0 - 1) * n^(-r)

        switch(EXPR = method,
               sample = {
                 sum(z * log(z)) / sum(z) - mean(log(z))
               },
               sdist_median = {
                 .MLEw_approx$fun$w2F(nz)
               },
               hybrid = {
                 # mean betw med-approx and sample estimate
                 (.MLEw_approx$fun$w2F(nz) + sum(z * log(z)) / sum(z) - mean(log(z))) / 2L
               },
               stop("This method for W2-estimation is not handled here!", call. = FALSE)
        )
      } #nf w2F


      # w3 function. Not vectorized in argument k
      w3FF <- function(group = "x", method = c("sample", "sdist_median", "hybrid")) {
        method <- match.arg(method)

        if (!twoGroup || group != "y") {
          obs <- x
          z <- z_x
          W1 <- W1_x
        } else {
          obs <- y
          z <- z_y
          W1 <- W1_y
        }
        nObs <- length(obs)

        # catch all for n = 1
        if (nObs < 2L) return(function(k) 1)

        # fn of shape k
        switch(EXPR = method,
               sdist_median = {
                 .MLEw_approx$fun$w3FF(nObs)
               },
               sample = function(k) {
                 W1 * if (log(k) < -5) 1 else if (k==1) mean(1/z) else sum(1/z^(1/k)) / sum(z^((k-1)/k))
               },
               stop("This method for W3 approximation is not handled here!", call. = FALSE))
      }#nf w3FF

      # return list of weights
      list(W1 = c(x = W1_x, y = W1_y),
           W2 = c(x = w2F(z = z_x, method = method_w2),
                  y = if (twoGroup) w2F(z = z_y, method = method_w2)),
           #function(k) W1_x * if (log(k) < -5) 1 else if (k==1) mean(1/z_x) else sum(1/z_x^(1/k)) / sum(z_x^((k-1)/k)),
           W3 = purrr::compact(list(x = w3FF(group = "x", method = method_w3),
                                    y = if (twoGroup) w3FF(group = "y", method = method_w3))))
    })
  } #esle weights


  stopifnot(!twoPhase) #XXX not implemented yet!!


  # provide indices for x and for y
  # where to find the parameters per group in the parameter vector of the objective function
  extractParOptInd <- if (!twoGroup) {
    # single group!
    list(x = seq_along(trNames)) ## Cave: trNames reacts to twoPhase-setting (which I've not thought through, yet)
  } else {
    # two group!
    #XXX exponential && profiled: indices are not correct for two groups, yet!!
    #+(this would allow to run simul_test.R!) #YYY already done?!
    if (is.null(bind)) {
      switch(distO$dist,
             exponential = {
               if (profiled) list(x = c(1L), y = c(2L)) else list(x = c(1L, 2L), y = c(3L, 4L))
             },
             weibull = {
               if (profiled) list(x = c(1L, 2L), y = c(3L, 4L)) else
                 list(x = c(1L, 2L, 3L), y = c(4L, 5L, 6L))
             },
             normal = {
               stopifnot(!profiled)
               list(x = c(1L, 2L), y = c(3L, 4L))
             },
             stop(glue("Unsupported distribution {distO$dist}!"), call. = FALSE)
      )
    } else if (length(oNames) == length(bind)) {
      # twoGroup, but all parameters are bound!
      switch(distO$dist,
             exponential = {
               # profiled can actually be true (as each group leads to own scale/rate
               #+but it will be averaged (see mergePars!!)
               if (profiled) {
                 #warning("Did not expect `profiled=TRUE` and full bind on all parameters!", call. = FALSE)
                 list(x = c(1L), y = c(1L))
               } else list(x = c(1L, 2L), y = c(1L, 2L))
             },
             weibull = {
               if (profiled) {
                 #warning("Did not expect `profiled=TRUE` and full bind on all parameters!", call. = FALSE)
                 list(x = c(1L, 2L), y = c(1L, 2L))
               } else list(x = c(1L, 2L, 3L), y = c(1L, 2L, 3L))
             },
             normal = {
               stopifnot(!profiled)
               list(x = c(1L, 2L), y = c(1L, 2L))
             },
             stop(glue("Unsupported distribution {distO$dist}!"), call. = FALSE)
      )
    } else {
      # twoGroups & non-trivial bind

      local({
        # param with profiled=TRUE & transformed = FALSE removes the profile parameters (although it is original scale)
        #+ as we need the original parameter names without those of profiling
        oNamesFullProf <- if (!profiled) oNamesFull else
          distO$param(twoPhase = twoPhase, bind = bind, twoGroup = TRUE, profiled = TRUE, transformed = FALSE)

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
  }#esle twoGroup

  # provide indices for x and for y
  # where to find the parameters per group in the common parameter vector
  extractParInd <- if (!profiled) {
    # = optimization indices when no profiling
    extractParOptInd
  } else {
    # profiled!
    if (!twoGroup) {
      # profiled single group!
      list(x = seq_along(oNames)) ## Cave: oNames reacts to twoPhase-setting (which I've not thought through, yet)
      #if (distO$dist == 'exponential') list(x = c(1L, 2L)) else list(x = c(1L, 2L, 3L))
    } else {
      # twoGroup && profiled
      stopifnot(distO$dist != "normal")
      if (is.null(bind)) {
        if (distO$dist == 'exponential') list(x = c(1L, 2L), y = c(3L, 4L)) else
          # weibull
          list(x = c(1L, 2L, 3L), y = c(4L, 5L, 6L))
      } else {
        if (length(oNames) == length(bind)) {
          # profiled can actually be true (as each group leads to own scale/rate
          #+but it will be averaged (see mergePars!!)
          #warning("Did not expect `profiled=TRUE` and full bind on all parameters!", call. = FALSE)
          if (distO$dist == 'exponential') list(x = c(1L, 2L), y = c(1L, 2L)) else
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

  # parameter transformation matrices (for single group)
  paramTransf <- list(
    M = switch(distO$dist,
               exponential = matrix(c( 1, 0, 0, 0,
                                       0, 1, 0, 0,
                                       -1, 0, 1, 0,
                                       0, 0, 0, 1), nrow = 4L, byrow = TRUE,
                                    dimnames = list(c("delay1_tr", "rate1_tr",
                                                      "delay2_tr", "rate2_tr"))),
               weibull = matrix(c( 1, 0, 0, 0, 0, 0,
                                   0, 1, 0, 0, 0, 0,
                                   0, 0, 1, 0, 0, 0,
                                   -1, 0, 0, 1, 0, 0,
                                   0, 0, 0, 0, 1, 0,
                                   0, 0, 0, 0, 0, 1), nrow = 6L, byrow = TRUE,
                                dimnames = list(c("delay1_tr", "shape1_tr", "scale1_tr",
                                                  "delay2_tr", "shape2_tr", "scale2_tr"))),
               normal = matrix(c( 1, 0,
                                  0, 1), nrow = 2, byrow = TRUE,
                               dimnames = list(c("mean_tr", "sd_tr"))),
               stop("Unknown distribution!", call. = FALSE)
    ),
    Minv = switch(distO$dist,
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
                  normal = matrix(c( 1, 0,
                                     0, 1), nrow = 2, byrow = TRUE,
                                  dimnames = list(c("mean", "sd"))),
                  stop("Unknown distribution", call. = FALSE)
    ),
    F = list(exponential = c(identity, log, log, log),
             weibull = c(identity, log, #log1p, #identity, #=shape1
                         log, log, log, log),
             normal = c(identity, identity))[[distO$dist]],
    Finv = list(exponential = c(identity, exp, exp, exp),
                weibull = c(identity, exp, #expm1, #identity, #=shape1
                            exp, exp, exp, exp),
                normal = c(identity, identity))[[distO$dist]]
  )

  # transform parameter vector for a single group. Does not use parameter names.
  # transformed parameters are used within optimization.
  # The transformation helps to ensure side-conditions (e.g. log-transformation ensures non-negativity of original parameter)
  # @param parV1 parameter vector for a single group
  # @param inverse logical. `inverse=TRUE` takes optimization parameters back to original parameters
  # @return transformed parameter vector, unnamed!
  transformPars1 <- function(parV1, inverse = FALSE) {

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
  # @param isOpt logical. Are the parameters on optimization scale?
  # @return merged parameter vector
  mergePars <- function(parx, pary, isOpt) {
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
                                  sqrt(x[[1L]] * x[[2L]]) #prod(x)^(1/length(x))
                            }, simplify = TRUE))
    # .. only for delay1 we use minimum as aggregation function
    #+(in this case first entry in exParInd$x and exParInd$y is 1!)
    if (distO$hasDelay && exParInd$x[[1L]] + exParInd$y[[1L]] == 2) {
      res[[1L]] <- min(parx[[1L]], pary[[1L]])
    }

    res
  }

  # Extract parameter vector for a specified group
  # if parameters are for optimization and transformation is requested, profiling is undone (if relevant)
  # @param group character. Extract parameters for the given group. If NULL, keep all parameters.
  # @param isOpt logical. Are the given parameters on optimization function scale?
  # @param transform logical. Transform parameters?
  # @param named logical. Extract parameters as named vector?
  # @return parameter vector
  extractPars <- function(parV, group = NULL, isOpt = TRUE, transform = FALSE, named = FALSE) {
    if (is.null(parV)) return(NULL)
    # result is on optimization scale?
    resIsOpt <- xor(isOpt, transform) #TRUE if different

    # basically, ignore group= when single group: use always canonical "x" then
    if (!twoGroup) group <- "x"

    if (is.null(group)) {
      return(local({

        # recursive calls for the individual groups
        parx <- extractPars(parV, group = "x", isOpt = isOpt, transform = transform, named = FALSE)
        pary <- extractPars(parV, group = "y", isOpt = isOpt, transform = transform, named = FALSE)

        # merge the two parameter vectors back together (after a potential transformation)
        res0 <- mergePars(parx = parx, pary = pary, isOpt = resIsOpt)
        if (named) {
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
          stopifnot(distO$dist != 'normal')
          # un-profile (when going from profiled par_opt to par_orig)
          if (isOpt) {
            # access observations for specified group
            obs <- if (group == "y") y else x
            k <- if (distO$dist == 'weibull') res0[[2L]] else 1L
            # calculate scale parameter
            scale0 <- if (!isSurv) {
              (mean((obs-res0[[1L]])^k) / weights$W1[[group]])^(1/k)
            } else {
              # Surv: only right-censored observations currently implemented!
              stopifnot(attr(obs, which = "type", exact = TRUE) == 'right')
              # we used to consider all times (including censorings) but would divide (for mean) only by the number of events
              # this seems wrong and could fail when first obs is right-censored prior to delay candidate
              #(mean((obs[,1L]-res0[[1L]])^k) * length(obs)/(length(obs) - cens$n[[group]][["right"]]) / weights$W1[[group]])^(1/k)

              # consider only event times: we take them from the KM-fit
              (mean((summary(kmFit)$time - res0[[1L]])^k) / weights$W1[[group]])^(1/k)
            }
            # add scale/rate parameter at the end of parameter vector
            res0 <- append(res0, values = if (distO$dist == 'exponential') 1/scale0 else scale0)
          } else {
            # extract only remaining parameters
            res0 <- res0[extractParOptInd[[group]]]
          }
        }# profiled

        res0
      })
    }


    # single group names
    if (named) {
      rlang::set_names(res, nm = if (resIsOpt) trNames else oNames)
    } else {
      as.vector(res)
    }
  } #fn extractPars


  # optimization arguments -----

  # Get optimization start values and upper limits based on observations from a single group
  # for `twoPhase=TRUE` there will be more parameters
  # with profiling no scale parameter is returned (as it is not optimized)
  # @param obs vector of observations from single group
  # @return list with transformed par for single group and upper limits for delay parameters, in canonical order (bind has no effect here!)
  getParSetting.gr <- function(obs) {
    # contract: obs is sorted!
    DELAY_MIN <- .Machine$double.xmin ##1e-9

    if (isSurv && (cens$n$x[["left"]] %||% 0) + (cens$n$y[["left"]] %||% 0) > 0) {
      stop("Left-censoring is not supported here!", call. = FALSE)
    }
    # extract first event time (we assume there is no left-censoring!)
    firstEvTime <- if (isSurv) obs[which(obs[, "status"] == 1)[[1L]], 1L] else obs[[1L]]

    # Surv: convert to numeric, quick fix, use only event times as numeric vector that are observed or right censored
    # XXX improve here?, e.g., use flatten_surv from lme4cens?! # could use cens-list here
    if (isSurv) {
      obs <- obs[, 1L, drop=TRUE][obs[, "status", drop = TRUE] <= 1]
    }

    parV <- switch(EXPR = distO$dist,
                   # min(obs) = obs[1L]
                   exponential = {
                     parV0 <- c(max(DELAY_MIN, obs[[1L]] - 3 / (length(obs)+1)),
                                mean(obs - obs[[1L]] + 2 / length(obs))^-1L)

                     # two extra parameters when exponential with *two* phases
                     if (twoPhase) parV0 <- c(parV0, obs[[floor(.5 + length(obs)/2L)]], parV0[[2L]])

                     #parV0 <- rlang::set_names(parV0, nm = oNames)
                     # transform start-parameters for optfun-parametrization
                     parV0 <- transformPars1(parV0, inverse = FALSE)

                     # drop scale if profiling
                     if (profiled) {
                       stopifnot(!twoPhase) #XXX not implemented, yet!
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

                     parV0 <- local({
                       start_delay <- max(DELAY_MIN, obs[[1L]] - 3 / (length(obs)+1))

                       # use median rank approximation for empirical Weibull CDF: F(i,n) = (i - 0.3) / (n + 0.4)
                       # and then ordinate is log(1/(1-F)) = -log(1-F) on log-scale
                       start_y <- log(-log(1-stats::ppoints(n = length(obs), a=.3)))

                       # log of centred observations
                       # avoid negative values (as DELAY_MIN is positive)
                       lobs0 <- log(pmax.int(DELAY_MIN, obs-start_delay))

                       # simple linear regression of Y=start_y vs X=log-obs
                       # cf. lm.fit(x = cbind(1, log(obs)), y = start_y)$coefficients
                       # weighted version with more weight in the middle:
                       # w <- seq_along(obs); w <- w * (max(w)+1-w) #or use plogis-weights to downweight the early obs
                       # lm.wfit(x = cbind(1, log(obs)), y = start_y, w = plogis(-2:(length(obs)-3)))$coefficients
                       start_shape <- stats::cor(lobs0, start_y) * stats::sd(start_y) / stats::sd(lobs0)
                       start_scale <- exp(mean(lobs0) - mean(start_y) / start_shape) # scale from intercept

                       c(start_delay, start_shape, start_scale)
                     })

                     if (verbose > 3) {
                       cat("Start values 1st phase (a single group): ",
                           paste(c("delay", "shape", "scale"),
                                 round(parV0, 2), sep = ": ", collapse = ", "),
                           "\n")
                     } #fi verbose

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
                   normal = {
                     # robust start values
                     # IQR in normal is 1.349 times the std. deviation
                     c(stats::median(obs), stats::IQR(obs)/1.349)
                   },
                   # default:
                   stop(glue("Distribution {sQuote(distO$dist)} is not implemented!"), call. = FALSE)
    )

    list(
      par = parV,
      delay1_upper = max(DELAY_MIN, firstEvTime - .01/length(obs), firstEvTime * .9999),
      delay2_upper = log(max(DELAY_MIN, obs[[length(obs)]] - .02/length(obs), obs[[length(obs)]]*.999))
    )
  }# fn getParSetting.gr

  # profile likelihood: maximize profiled log-lik f directly
  # if FALSE, go indirectly: consider min(f'^2) to hunt for *local* extremum as these local extrema have f'^2 == 0 as necessary condition
  #profiled_llik_directly <- TRUE

  # parameter bounds: set lower & upper bounds
  lowerB <- upperB <- rlang::rep_named(names = trNamesFull, x = NA_real_)


  #XXX #QQQ Should this go up to extractPars-function where the transformations are defined???
  PAR_BOUNDS <- list(delay1 = c(lower = 0, upper = NA_real_),
                     delay2 = c(lower = -Inf, upper = NA_real_),
                     rate  = c(lower = -Inf, upper = +Inf),
                     # shape lower bound for MLEnp (actually for shape1)
                     #shape = c(lower = if (profiled && method == 'MLEn' && !profiled_llik_directly) 1.49e-8 else -Inf, upper = +Inf),
                     shape = c(lower = -Inf, upper = 3.5), # exp(3.5) = 33 is already huge for shape [exp(1.6) = 5, exp(4.5) = 90]
                     scale = c(lower = -Inf, upper = +Inf),
                     mean = c(lower = -Inf, upper = +Inf),
                     sd = c(lower = 0, upper = +Inf))


  # set bounds from lookup table PAR_BOUNDS
  # alas, purrr::iwalk did not work for me here
  for (nam in names(PAR_BOUNDS)) {
    idx <- startsWith(trNamesFull, prefix = nam)
    if (any(idx)) {
      lowerB[idx] <- purrr::chuck(PAR_BOUNDS, nam, 'lower')
      upperB[idx] <- purrr::chuck(PAR_BOUNDS, nam, 'upper')
    } #fi
  } #rof



  par0_x <- getParSetting.gr(x)
  if (verbose > 2) cat("Start parameters for opt, group x: ",
                       paste(round(par0_x$par, 3), collapse = ", "), "\n")
  parV <-
    if (! twoGroup) {
      # set parameter vector for group 1 and finish upper bound: match delay1 and delay2
      upperB[['delay1_tr']]  <- par0_x[['delay1_upper']]
      if (twoPhase) upperB[['delay2_tr']] <- par0_x[['delay2_upper']]

      par0_x[['par']]

    } else { #twoGroup

      # all parameters are bound
      if (length(bind) == length(oNames)) {

        # treat x and y as a single group for upper limit & start value heuristic
        par0_xy <- getParSetting.gr(c(x,y))

        upperB['delay1_tr'] <- par0_xy[['delay1_upper']]
        if (twoPhase) upperB[['delay2_tr']] <- par0_xy[['delay2_upper']]

        par0_xy[['par']]

      } else {
        #twoGroup, not all params bound!
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
          if ('delay2' %in% bind) {
            upperB[['delay2_tr']] <- max(par0_x[['delay2_upper']], par0_y[['delay2_upper']])
          } else {
            upperB['delay2_tr.x'] <- par0_x[['delay2_upper']]
            upperB['delay2_tr.y'] <- par0_y[['delay2_upper']]
          }
        } #fi twoPhase

        # return start value
        if (is.null(bind)) { # two groups unbound
          c(start_x, start_y)
        } else {

          mergePars(parx = start_x, pary = start_y, isOpt = TRUE)
        }
      } #twoGrp, not all params bound!
    } # twoGrp

  # ensure we have names of transformed parameters
  parV <- rlang::set_names(parV, nm = trNamesFull)

  stopifnot(!any(is.na(lowerB), is.na(upperB)))
  # clean up env. # use local() more???
  remove(list=c("PAR_BOUNDS", "par0_x"))

  if (verbose > 1L) {
    cat("Start values for opt: ",
        paste(names(parV), round(parV, 3), sep = "=", collapse = ", "), "\n")
  }

  optim_args <- list(
    par = parV,
    method = "L-BFGS-B",
    lower = lowerB,
    upper = upperB,
    # most parameters are on log-scale.
    control = list(parscale = scalePars(parV, lowerB = 1e-3, upperB = 1e3))
  )


  # Penalization for high values of shape per group
  #
  # For Weibull distribution in MLEw method it penalizes high shape values.
  # The penalization factor increases with sample size as
  # location (median) and spread (mad) of the ML-objective function grow with sample size:
  # most clearly so for MLEn and MLEc. MLEw is less regular (maybe due to bad fits)
  #
  # Older idea was to use start values and estimate lowish average of log-density
  # for the observed values in the group (or: for range of possible values)
  #+But: important is not so much the base level of log-likelihood
  #+but what reduction is possible through optimization of start values, no?!
  #
  # @seealso simulations in `MLEw_shape_penalization.R`
  # @param k candidate value for shape
  # @param nObs number of observations in group
  # @return non-negative penalty value. Big values mean higher penalty (it gets subtracted from the log-likelihood)
  penF <- function(k, nObs = 1) {
    # FALSE #to turn off penalty
    pen_shape <- distO$dist == 'weibull' && method == 'MLEw'

    if (!isTRUE(pen_shape)) return(0)

    pen_shape_shift <- 9.9 #shift parameter of softplus penalty
    pen_shape_steep <- .9 #steepness of softplus penality

    nObs * log1p(exp(pen_shape_steep * (k - pen_shape_shift)) / pen_shape_steep)
  }#fn penF


  # objective function ----

  # calculate the log-likelihood, either naive, weighted or in corrected form.
  # What precisely is calculated depends on its surrounding closure (value of method but also the profiled-flag)
  # @param pars complete vector of parameters (can be refering to two groups)
  # @param criterion logical. If `TRUE`, then pars are on original scale and the proper log-likelihood is returned. This flag currently serves a double purpose! (Disentangle maybe?)
  getLogLik <- function(pars, group, criterion = FALSE) {

    # Old idea was to
    # change signature to be with pars.gr and obs for both getLogLik and getCumDiffs
    #+But what are the benefits?

    # access observations of group
    obs <- if (group == "y") y else x #direct access by name
    #obs <- rlang::env_get(env = rlang::env_parent(rlang::current_env(), n=1L), nm = group, inherit = FALSE)

    # extract parameters for specified group on original scale (for CDF)
    pars.gr <- extractPars(pars, group = group, isOpt = !criterion, transform = !criterion)


    if (criterion) {
      # criterion = proper log-likelihood
      retV <- if (!isSurv) {
        # numeric response, non-Surv
        sum(rlang::exec(distO$pdf, !!! c(list(x=obs, log=TRUE), pars.gr)))
      } else {
        # Surv-response
        #if (verbose > 2) cat(glue("Parameter {paste(pars.gr, collapse = '; ')}"))
        switch (attr(obs, which = "type", exact = TRUE),
                right = {
                  sum(rlang::exec(distO$pdf, !!! c(list(x=obs[cens$ind[[group]]$obs,  1L], log=TRUE), pars.gr)),
                      rlang::exec(distO$cdf, !!! c(list(q=obs[cens$ind[[group]]$right,1L], lower.tail = FALSE, log.p = TRUE), pars.gr)))
                },
                stop("This type of censoring is not supported!", call. = FALSE)
        )
      } #esle
      return(retV)
    } #fi criterion


    # !criterion
    # calculate the objective function which depends on
    #+method
    #+profiled
    nObs <- length(obs)
    stopifnot(!criterion, nObs > 1L)

    # shape parameter (candidate)
    k <- if (distO$dist == 'weibull') pars.gr[[2L]] else 1L

    switch(EXPR = method,
           MLEn = {
             if (profiled && distO$dist == 'weibull') {
               if (isSurv) {
                 switch(attr(obs, which = "type", exact = TRUE),
                        right = {
                          obs_c <- obs[cens$ind[[group]]$obs, 1L] - pars.gr[[1L]]

                          # objective function to maximize:
                          # we use 1st derivative to profile out scale parameter, but otherwise, use log-likelihood function directly
                          (nObs - cens$n[[group]][["right"]]) * ((k-1) * mean(log(obs_c)) - log(mean(obs_c^k)) + log(k) - 1) +
                            # contribution of right censorings
                            sum(rlang::exec(distO$cdf, !!! c(list(q=obs[cens$ind[[group]]$right,1L], lower.tail = FALSE, log.p = TRUE), pars.gr))) +
                            # optional penalty term for large values of shape
                            -penF(k, nObs = nObs)

                        },
                        stop("This Surv-type is not supported!", call. = FALSE)
                 )
               } else {
                 # numeric response, non-Surv
                 obs_c <- obs - pars.gr[[1L]]
                 #cat("\nDelay a: ", pars.gr[["delay1"]], "Shape k: ", k, " (", pars[2], ")\n") #DDD debug

                 # return early when we have too high delay parameter
                 if (obs_c[[1L]] < 0) return(NA_real_)
                 # objective function to maximize:
                 # we use 1st derivative to profile out scale parameter but use log-likelihood function directly otherwise
                 # 2nd & 3rd summand could also be: - log(sum(obs_c^k)) + log(n*k)
                 nObs * ((k-1) * mean(log(obs_c)) - log(mean(obs_c^k)) + log(k) - 1) - penF(k, nObs = nObs)

                 # alternative:
                 #indirect way: ! profiled_llik_directly
                 #consider min(f'^2) to hunt for *local* extremum as these local extrema have f'^2 == 0 as necessary condition
                 #We would need to check that we have indeed an local **maximum** for the log-likelihood (as we have only found candidate values by looking for roots of f')
                 #   - (1/k + mean(log(obs_c)) - sum(log(obs_c) * obs_c^k) / sum(obs_c^k))^2 -
                 #     # 1st factor is inverse of harmonic mean
                 #     (mean(1/obs_c) * sum(obs_c^k)/sum(obs_c^(k-1)) - k/(k-1))^2 -
                 #     # optional penalization term
                 #     penalize_shape*log(k+1)
               } #esle

             } else {
               # log-likelihood with all parameters (scale is not profiled out)
               if (!isSurv) {
                 # numeric, non-Surv
                 sum(rlang::exec(distO$pdf, !!! c(list(x=obs, log=TRUE), pars.gr))) +
                   -penF(k, nObs = nObs)
               } else {
                 # Surv-response
                 switch (attr(obs, which = "type", exact = TRUE),
                         right = {
                           sum(rlang::exec(distO$pdf, !!! c(list(x=obs[cens$ind[[group]]$obs,  1L], log=TRUE), pars.gr)),
                               rlang::exec(distO$cdf, !!! c(list(q=obs[cens$ind[[group]]$right,1L], lower.tail = FALSE, log.p = TRUE), pars.gr))) +
                             # optional penalty term for high shape parameter
                             -penF(k, nObs = nObs)
                         },
                         stop("This Surv-type is not supported!", call. = FALSE))
               } #esle !isSurv
             } #esle
           },

           # weighted MLE
           MLEw = {
             stopifnot(profiled)
             stopifnot(distO$hasDelay, distO$dist != "normal")

             retVal <- if (!isSurv) {
               # numeric response, non-Surv
               obs_c <- obs - pars.gr[[1L]]

               # objective function to maximize
               -(weights$W2[[group]] / k + mean(log(obs_c)) - sum(log(obs_c) * obs_c^k) / sum(obs_c^k))^2 +
                 # 1st factor is inverse of harmonic mean
                 -(mean(1/obs_c) * sum(obs_c^k) / sum(obs_c^(k-1)) - weights$W3[[group]](k))^2 +
                 # optional penalization term
                 -penF(k, nObs = nObs)

             } else {
               # Surv-response
               switch(EXPR = attr(obs, which = "type", exact = TRUE),
                      right = {
                        obs_evc <- obs[cens$ind[[group]]$obs, 1L] - pars.gr[[1L]]

                        # objective function to maximize
                        -(weights$W2[[group]] / k + mean(log(obs_evc)) - sum(log(obs_evc) * obs_evc^k) / sum(obs_evc^k))^2 +
                          # 1st factor is inverse of harmonic mean
                          -(mean(1/obs_evc) * sum(obs_evc^k) / sum(obs_evc^(k-1)) - weights$W3[[group]](k))^2 +
                          # contribution of right-censored obs
                          rlang::exec(distO$cdf,  !!! c(list(q=obs[cens$ind[[group]]$right,1L], lower.tail = FALSE, log.p = TRUE), pars.gr)) +
                          # optional penalization term for big shape
                          -penF(k, nObs = nObs)
                      },
                      stop("This Surv-type is not supported here!", call. = FALSE))
             } #esle !isSurv

             if (verbose > 1L) {
               cat(glue("W1 = {round(weights$W1[[group]],2)}, ",
                        "W2 = {round(weights$W2[[group]],2)}, ",
                        "W3 = {round(weights$W3[[group]](k),4)} for {group}. ",
                        "Candidate values: delay {round(pars.gr[[1L]],3)} shape {round(k,3)} ",
                        "=> LLval: {round(retVal, 3)}"),
                   "\n")
             }

             retVal
           },

           # corrected MLE
           # objective function to maximize
           MLEc = {
             stopifnot(nObs >= 2L)
             # contribution of first observation is corrected for: we take first two different values
             ind12 <- indForefront[[group]]

             if (!isSurv) {
               # numeric response, non-Surv
               sum(length(ind12[["inds_obs1"]]) * log(diff(rlang::exec(distO$cdf, !!! c(list(q=obs[c(1L, ind12[["ind_next"]])]), pars.gr)))),
                   rlang::exec(distO$pdf, !!! c(list(x=obs[-ind12[["inds_obs1"]]], log=TRUE), pars.gr)),
                   # optional penalization term
                   -penF(k, nObs = nObs))

             } else {
               switch (attr(obs, which = "type", exact = TRUE),
                       right = {
                         # we need at least two observed event times
                         stopifnot(nObs - cens$n[[group]]["any"] >= 2L)

                         # first event-time needs correction
                         sum(length(ind12[["inds_obs1"]]) * log(diff(rlang::exec(distO$cdf, !!! c(list(q=obs[c(ind12[["inds_obs1"]][1L], ind12[["ind_next"]]),1L]),
                                                                                                  pars.gr)))),
                             # remaining observed event times
                             rlang::exec(distO$pdf, !!! c(list(x=obs[setdiff(cens$ind[[group]]$obs, ind12[["inds_obs1"]]),1L], log=TRUE), pars.gr)),
                             # right-censored observations do not need correction
                             #+(as they are tail probabilities that do not peak so drastically as densities do)
                             rlang::exec(distO$cdf, !!! c(list(q=obs[cens$ind[[group]]$right,1L], lower.tail = FALSE, log.p = TRUE), pars.gr)),
                             # optional penalization term
                             -penF(k, nObs = nObs))
                       },
                       stop("This type of censoring is not supported!", call. = FALSE))

             } #esle !isSurv
           },
           stop(glue("This method {method} is not handled here!"), call. = FALSE)
    )
  }


  # log spacings:
  # calculate the differences in EDF (for given parameters in group) of adjacent observations on log scale
  # These log-spacings are the heart of the MPSE-criterion which is the negative mean of these log-spacings.
  # Moran's test statistic is the negative sum of these log-spacings.
  # @param pars vector of parameters (by default, on transformed scale, i.e. when criterion = FALSE)
  # @param criterion logical. When `criterion = TRUE`, then pars are on original scale. No other meaning here.
  # @param ties. how to handle ties. By default, use the tie-setting from objective function call.
  # @return n+1 cumulative diffs on log-scale (or single negative number in twoPhase when delay2 <= delay in quick fix)
  getCumDiffs <- function(pars, group, criterion = FALSE, ties. = ties) {

    # access observations of group
    #obs <- rlang::env_get(env = rlang::env_parent(rlang::current_env(), n=1L), nm = group, inherit = FALSE)
    #+or use env = rlang::fn_env(getCumDiffs) # (but requires function obj)
    obs <- if (group == "y") y else x # direct access by name

    # extract parameters for specified group on original scale (for CDF)
    pars.gr <- extractPars(pars, group = group, isOpt = !criterion, transform = !criterion)

    if (verbose > 1L) {
      cat(glue("Parameter vector for group {group} on non-transformed scale: ",
               "{paste(round(pars.gr, 2), collapse = ', ')}"), "\n")
    }



    # calculate spacings
    # contract: data is sorted!
    cumDiffs <- if (!isSurv) {
      # numeric response (non-Surv)
      diff(c(0L, rlang::exec(distO$cdf, !!! c(list(q=obs), pars.gr)), 1L))
    } else {
      # Surv-response
      ind_evKM <- which(kmFit$n.event > 0.99) #at least one event (type="interval" makes that we get fractional numbers here [but 0 is 0 also for interval!?])

      # get the right subset of indices for specified group (when having two groups)
      if (twoGroup) {
        # Cave: works only for two groups (x or y) as I only use the strata[[1L]] as cutpoint
        ind_evKM <- if (group == "x") ind_evKM[ind_evKM <= kmFit$strata[[1L]]] else ind_evKM[ind_evKM > kmFit$strata[[1L]]]
      }
      # h: return object
      h <- rep_len(-1, length.out = length(obs))

      # n.event is generally not integer for type=interval/left. It is increased by a fraction (depending on number of events) and sums to nbr of events+1 (per group)
      # floor(n.event + n.censor) == n
      stopifnot(sum(as.integer(kmFit$n.event[ind_evKM]),
                    if (twoGroup) kmFit$n.censor[c(-1,1)[[1L+(group == "x")]] * seq_len(kmFit$strata[[1L]])] else kmFit$n.censor) == kmFit$n[[if (group == "x") 1L else 2L]])
      # use CDF for all observed event times of right-censored outcome variable
      h[obs[, "status"] == 1] <- rep.int(rlang::exec(distO$cdf, !!! c(list(q=kmFit$time[ind_evKM]), pars.gr)) * (cens$rcens$surv[ind_evKM]) + (1 - cens$rcens$surv[ind_evKM]),
                                         # n.event is not always integer for Surv-type=interval/left. rep.int truncates floats & it should always work.
                                         times = kmFit$n.event[ind_evKM])
      # censored observations get interpolated values
      ind_hrcens <- which(h<0)
      if (length(ind_hrcens)) {
        ind_hobs <- which(h>0)
        # interpolate values for all censored observations
        h[ind_hrcens] <- stats::approx(x = c(0L, ind_hobs, length(obs)+1L), y = c(0L, h[ind_hobs], 1L),
                                       method = "linear", ties = "ordered", # x-values are already ordered!
                                       yleft = NA, yright = NA,
                                       # values where to interpolate
                                       xout = ind_hrcens)$y
      } #fi hrcens

      diff(c(0L, h, 1L))

    }#esle !isSurv

    # check for ties to fix cumDiffs for observed event times
    tig <- tieInfo[[group]] # tie info group (tig)
    nTigs <- NROW(tig[["tieGrp"]])

    if (nTigs) {
      stopifnot(all(cumDiffs[tig[["cumDiffInd"]]] == 0)) # all spacings for tied observed event times are 0

      obsVals <- if (!isSurv) obs[tig$tieGrp[, "startInd"]] else obs[tig$tieGrp[, "startInd"], 1L]

      cumDiffs[tig[["cumDiffInd"]]] <- switch(
        ties.,
        density = {
          # use density instead of diff of CDF for tied observation pairs
          rep.int(
            # take first observation per tie-group
            rlang::exec(distO$pdf, !!! c(list(x = obsVals), pars.gr)),
            #tig$tieGrp[, "len"]-1L # number of repeats per tie group
            times = tig$tieGrp[, "len"]-1L)
        },
        # "equispaced" for CDF-backtransformed using given parameters, then equispaced spacings across tie groups
        # we keep using a standard density strategy for fitting, but can request another tie-strategy for evaluating the MPSE-criterion,
        #e.g., for Moran's test
        equispaced = {
          # per tie group, assume tied observations are maximally spread (within rounding radius).
          # Two reasons why this is leads to bigger cumDiffs (=smaller criterion/Moran's test statistic = conservative)
          # 1/ for adjacent spacings (involving obs directly before and after tie) we have assumed the original tied observation
          # 2/ we use equal spacings in transformed space for all tied observation (within tie group)
          rep.int(diff(rlang::exec(distO$cdf,
                                   !!! c(list(q = rep(obsVals, each = 2L) + c(-1, 1) * tig$numPrecision[["rRad"]]), pars.gr)))[seq.int(from = 1, by = 2, length.out = nTigs)] / (tig$tieGrp[, "len"]-1L),
                  times = tig$tieGrp[, "len"]-1L)

        },
        # for what it's worth: tie-strategy "error" should have already quit
        error = {
          stop("getCumDiffs: ties are not allowed!", call. = FALSE)
        },
        # handle exception
        stop(glue("Unknown strategy {ties.} to handle ties here."), call. = FALSE)
      )
    } #fi

    # respect the machine's numerical lower limit
    cumDiffs[which(cumDiffs < .Machine$double.xmin)] <- .Machine$double.xmin

    log(cumDiffs)

  }# fn getCumDiffs


  # Objective function to be minimized.
  #
  # Depending on method, it is negative mean log-spacings for MPSE or negative log-likelihood for MLEn
  # One can estimate parameters by minimizing this objective function.
  #
  # param `pars` the vector of parameters. transformed when criterion=FALSE and not transformed when criterion=TRUE
  # param `criterion` logical. If `TRUE`, give the original criterion to minimize. In this case, the parameters must be on original scale.
  # param `aggregated` logical. For two group case, `aggregated=FALSE` returns values per group, like mean log cum-diffs per group.
  # param `ties.` how to handle ties for the MPSE-function. Default value is 'density'.
  objFun <- function(pars, criterion = FALSE, aggregated = TRUE, ties. = ties) {

    if (verbose > 1) cat("pars:", pars, "\n")

    retVal <- switch(method,
                     MPSE = {
                       if (!twoGroup) {
                         -mean(getCumDiffs(pars, group = "x", criterion = criterion, ties. = ties.))
                       } else {
                         #twoGroup:
                         #the approach to first merge x and y and then do the cumDiffs, log and mean does *not* work out
                         #because the parameters should be optimized within group.
                         #merged data lead to frequent non-convergence or visually bad fits
                         res <- c(mean(getCumDiffs(pars, group = "x", criterion = criterion, ties. = ties.)),
                                  mean(getCumDiffs(pars, group = "y", criterion = criterion, ties. = ties.)))

                         if (aggregated) -stats::weighted.mean(res, w = c(length(x), length(y))) else -res
                       }
                     },
                     MLEn = ,
                     MLEw = ,
                     MLEc = {
                       stopifnot(!twoPhase) #XXX not implemented yet!

                       if (!twoGroup) -getLogLik(pars, group = "x", criterion = criterion) else {
                         res <- c(getLogLik(pars, group = "x", criterion = criterion), getLogLik(pars, group = "y", criterion = criterion))

                         if (aggregated) -sum(res) else -res
                       }

                     },
                     stop(glue('Objective function for method {method} is not implemented!'), call. = FALSE)
    )

    if (verbose > 2) {
      cat("Objfun value: ", retVal, "\n")
    }
    retVal
  } #fn objFun

  # attach analytical solution for MLE
  if (method == 'MLEn' && !twoGroup && !twoPhase && distO$dist == 'exponential' && !isSurv) {
    attr(objFun, which = "opt") <- local({
      par_analytic <- c(delay1 = x[[1L]], rate1 = 1L/(mean(x) - x[[1L]]))
      list(par_orig = par_analytic,
           #transformed parameters
           par = extractPars(par_analytic, isOpt = FALSE, transform = TRUE, named = TRUE),
           value = length(x) * (log(mean(x) - x[[1L]]) + 1L),
           methodOpt = "analytic",
           convergence = 0L,
           message = "analytic solution for naive MLE ('MLEn')",
           counts = 0L)
    })
  }

  objFun
} #fn objFunFactory




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
  stopifnot(is.function(objFun))
  objFunEnv <- rlang::fn_env(objFun)

  # gather information from objective function environment
  objFunObjs <- rlang::env_get_list(env = objFunEnv,
                                    nms = c("bind", "method", "optim_args" ,"trNamesFull", "profiled", "twoGroup", "twoPhase", "x", "y", "extractPars", "oNames"))

  # check if there is already a solution provided by the objective function
  optObj <- attr(objFun, which = "opt", exact = TRUE)

  if (is.list(optObj) && all(c("par", "par_orig", "value", "convergence") %in% names(optObj))) {
    if (verbose > 0L) message("Using provided (analytical) solution to objective function.")
  } else {
    optObj <- NULL #start from scratch
    # numeric optimization
    if (verbose > 0L) message("Start with numeric optimiziation of objective function.")

    if (is.null(optim_args)) {
      # set standard optim-args
      optim_args <- objFunObjs[["optim_args"]]
    }

    stopifnot(is.list(optim_args), "par" %in% names(optim_args),
              is.numeric(optim_args$par), length(optim_args$par) == length(objFunObjs$trNamesFull))
    # ensure that transformed parameters are named
    if (!rlang::is_named(optim_args$par)) rlang::names2(optim_args$par) <- objFunObjs$trNamesFull
    stopifnot(identical(names(optim_args$par), objFunObjs$trNamesFull))
    # set objective function (overwrite entry 'fn' if it is already present)
    optim_args[["fn"]] <- objFun


    # optim: first attempts ----

    # initial start values for optimization
    par0 <- optim_args$par

    try({
      optObj <- rlang::exec(stats::optim, !!! optim_args)
      optObj$methodOpt <- optim_args$method
    }, silent = TRUE)


    if (is.null(optObj)) {
      if (verbose > 0L) warning(glue("{objFunObjs$method}-optimization failed during model fit!"),
                                call. = FALSE)
    } else if (isTRUE(optObj$convergence > 0L)) {
      # do a 2nd attempt of optim in case it did not converge in the first place
      if (verbose > 1L) message("No proper convergence during 1st optimization in delay fit. Re-try with different parameter scaling.")

      # Use parameter values of non-converged fit as new start values (and adapt parscale accordingly)
      #+The objFun is to be minimized, smaller is better!
      if (isTRUE(is.numeric(optObj$par) && all(is.finite(optObj$par)) && optObj$value < objFun(par0))) {
        if (verbose > 1L) cat("Set new start values for 2nd attempt\n")
        optim_args[["par"]] <- optObj$par  # purrr::assign_in(where = "par", value = optObj$par)

        if ("parscale" %in% names(optim_args[["control"]])) {
          optim_args[["control"]][["parscale"]] <- scalePars(optim_args[["par"]])
        }

        # optim: 2nd attempt --
        optObj <- NULL
        if (verbose > 1L) message("Do 2nd attempt with renewed start values and parameter scaling")
        try({
          optObj <- purrr::exec(stats::optim, !!! optim_args)
          optObj$methodOpt <- optim_args$method
        }, silent = TRUE)

        if (verbose > 0L && (is.null(optObj) || isTRUE(optObj$convergence > 0L))) {
          warning("No proper convergence after re-try.", call. = FALSE)
        }
      }## fi rescaling for 2nd attempt
    }## fi 2nd attempt necessary?


    # nlminb (PORT): last attempt ----

    if (is.null(optObj) || optObj$convergence > 0L) {
      if (verbose > 0L) cat("Do another final attempt with PORT-optimizer.\n")

      # choose best start values for PORT:
      # if there are shape parameters, go for start value that is reasonably small
      par1 <- local({
        shapeInd <- which(startsWith(names(par0), prefix = "shape"))
        keep0 <- length(shapeInd) && sum(pmax.int(par0[shapeInd]-2,0)^2) < sum(pmax.int(optim_args$par[shapeInd]-2,0)^2)
        if (keep0) {
          if (verbose > 1) cat("Keep initial start parameters for final PORT-optimizer attempt.\n")
          par0
        } else {
          if (verbose > 1) cat("Use updated start parameters for final PORT-optimizer attempt.\n")
          optim_args$par
        } #esle
      })

      optim_args$par <- par1 #update optim_args
      optObj <- minObjFunPORT(objFun = objFun, start = optim_args$par,
                              lower = optim_args$lower, upper = optim_args$upper,
                              verbose = verbose)
    } #fi


    # post-process optObj -----

    # set names to parameter vector
    if (! is.null(optObj)) {
      stopifnot("par" %in% names(optObj))
      stopifnot(identical(names(optObj$par), objFunObjs$trNamesFull))
      # # set canonical names for parameters
      # optObj$par <- rlang::set_names(optObj$par, objFunObjs$trNamesFull)
      # save optim_args in optimization object (but w/o objective function)
      optim_args$fn <- NULL
      optObj <- append(optObj, values = list(optim_args = optim_args))
    } #fi

    # add par_orig
    optObj <- append(optObj,
                     values = list(par_orig = objFunObjs$extractPars(parV = optObj$par, group = NULL,
                                                                     isOpt = TRUE, transform = TRUE, named = TRUE))
    )
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
#' @param distribution Which delayed distribution is assumed? Exponential or Weibull. Can be given as character or as distribution object.
#' @param twoPhase logical. Allow for two phases?
#' @param bind character. parameter names that are bind together in 2-group situation.
#' @param ties character. Strategy to handle ties for `method = "MPSE"`.
#' @param method character. Which method to fit the model? 'MPSE' = maximum product of spacings estimation *or* 'MLEn' = naive maximum likelihood estimation *or* 'MLEw' = weighted MLE' *or* MLEc' = corrected MLE
#' @param profiled logical. Profile out scale from log-likelihood if possible.
#' @param optim_args list. optimization arguments to use. Use `NULL` to use the data-dependent default values.
#' @param verbose integer. level of verboseness. Default 0 is quiet.
#' @return `incubate_fit` the delay-model fit object. Or `NULL` if optimization failed (e.g. too few observations).
#' @export
delay_model <- function(x = stop('Specify observations for first group x=!', call. = FALSE), y = NULL,
                        distribution = c('exponential', 'weibull', 'normal'), twoPhase = FALSE,
                        bind = NULL, ties = c('density', 'equispaced', 'error'),
                        method = c('MPSE', 'MLEn', 'MLEw', 'MLEc'), profiled = method == 'MLEw',
                        optim_args = NULL, verbose = 0) {

  # setup -------------------------------------------------------------------

  if (is.logical(verbose)) verbose <- as.numeric(verbose)
  if (is.null(verbose) || !is.numeric(verbose) || !is.finite(verbose) ) verbose <- 0L
  verbose <- verbose[[1L]]


  # unpack x if it is a list of two vectors
  if (is.list(x)) {
    if (length(x) != 2L) {
      stop("If x= is given a list it must be of size 2.", call. = FALSE)
    }
    y <- x[[2L]]
    x <- x[[1L]]
  }

  # enforce that the first argument x= is properly instantiated
  stopifnot(!is.null(x), is.numeric(x), length(x) > 0)

  distO <- if (is.list(distribution)) distribution else buildDist(match.arg(distribution))

  method <- if (length(method) == 1L && toupper(method) == 'MSE') {
    message("The method name 'MPSE' is prefered over the previously used name 'MSE'!")
    "MPSE"
  } else method[1L]
  method <- match.arg(method)
  ties <- match.arg(ties)

  if (is.character(bind)) {
    if (any(endsWith(bind, suffix = "_tr"))) {
      stop("Parameter names to bind= refer to the distribution parameters and not to transformed parameters of the objective function.",
           call. = FALSE)
    }

    # translate convenience names (for single phase) to canonical names
    if (distO$twoPhaseAllowed) {
      unNmbrdIdx <- !grepl(pattern = "[12]", bind, fixed = FALSE)
      if (any(unNmbrdIdx)) {
        bind[unNmbrdIdx] <- paste0(bind[unNmbrdIdx], "1") #interpret un-numbered parameters as referring to phase 1
        if (verbose > 0L) {
          cat("The unnumbered parameter names in bind= are translated to canonical parameter names (=phase 1).\n")
        }
      }
    }#fi twoPhaseAllowed
  }#fi bind=


  # objective function ------------------------------------------------------

  objFun <- objFunFactory(x = x, y = y, method = method, profiled = profiled, distO = distO,
                          twoPhase = twoPhase, bind = bind, ties = ties, verbose = verbose)
  if (is.null(objFun)) return(invisible(NULL))
  objFunEnv <- rlang::fn_env(objFun)

  # optimise objective function
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj) || is.null(optObj$par_orig)) return(invisible(NULL))

  # return -----
  twoGroup <- rlang::env_get(env = objFunEnv, nm = "twoGroup")
  # overwrite data with  pre-processed data
  x <- rlang::env_get(env = objFunEnv, nm = "x")
  y <- rlang::env_get(env = objFunEnv, nm = "y", default = NULL)

  # /!\ keep in sync with update()!
  structure(
    list(
      data = if (twoGroup) list(x = x, y = y) else x,
      nobs = c(x = NROW(x), y = if (twoGroup) NROW(y) else 0),
      distO = distO,
      twoPhase = twoPhase,
      twoGroup = twoGroup,
      method = method,
      bind = rlang::env_get(env = objFunEnv, nm = "bind"),
      ties = ties,
      #isSurv = rlang::env_get(env = objFunEnv, nm = "isSurv"),
      cens = rlang::env_get(env = objFunEnv, nm = "cens", default = 0L), ##if (twoGroup)
      kmFit = rlang::env_get(env = objFunEnv, nm = "kmFit", default = NULL),
      objFun = objFun,
      par = optObj$par_orig,
      criterion = objFun(pars = optObj$par_orig, criterion = TRUE, aggregated = TRUE),
      optimizer = purrr::compact(c(list(parOpt = optObj$par, valOpt = optObj$value, profiled = profiled),
                                   optObj[c("methodOpt", 'convergence', 'message', 'counts', 'optim_args')]))),
    class = "incubate_fit")
}


#' Refit an `incubate_fit`-object with specified optimization arguments.
#' This function is useful when only an optimization argument is to be changed.
#' If more things need to be changed go back to `delay_model` and start from scratch.
#' @param object `incubate_fit`-object
#' @param optim_args optimization arguments
#' @param verbose integer flag. Requested verbosity during `delay_fit`
#' @param ... further arguments, currently not used.
#' @return The updated fitted object of class `incubate_fit` or `NULL` in case of failure.
#' @export
update.incubate_fit <- function(object, optim_args = NULL, verbose = 0, ...) {

  stopifnot(all(c("data", "distO", "method", "objFun", "twoPhase", "twoGroup", "par", "criterion", "optimizer") %in% names(object)))

  ## fit model with given optim_args
  objFun <- object[["objFun"]]
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))


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
simulate.incubate_fit <- function(object, nsim = 1, seed = NULL, ...) {
  stopifnot(inherits(object, "incubate_fit"))

  ranFun <- object$distO$random

  #XXX add option to mirror cens= setting in observed data?
  # arguments to the random function generation
  ranFunArgsX <- as.list(c(n=object$nobs[[1L]], coef(object, group = "x")))
  ranFunArgsY <- if (object$twoGroup) as.list(c(n=object$nobs[[2L]], coef(object, group = "y")))

  simExpr <- if (object$twoGroup) {
    expression(list(x=rlang::exec(ranFun, !!! ranFunArgsX),
                    y=rlang::exec(ranFun, !!! ranFunArgsY)))
  } else {
    expression(rlang::exec(ranFun, !!! ranFunArgsX))
  }

  if (nsim > 1000L) {
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
  ranFun <- object$distO$random
  dFun <- object$distO$pdf
  twoGroup <- isTRUE(object$twoGroup)
  nObs <- object$nobs
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
    if ( min(obs1, del_coef) < TOL_NUM ) return(rep_len(del_coef, length.out = R))

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
    delayCandDF$objValInv <- (delayCandDF$objValInv / (max(delayCandDF$objValInv, na.rm = TRUE) + .001))^(1L/(smd_factor+.01))
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
               statistic = function(d, i) coef(delay_model(x=d[i], distribution = object$distO, twoPhase = object$twoPhase,
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

                          coef(delay_model(x=x, y=y, distribution = object$distO, twoPhase = object$twoPhase,
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

                            dm <- suppressWarnings(delay_model(x=x, y=y, distribution = object$distO, twoPhase = object$twoPhase,
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
    #  FUN = \(d) coef(delay_model(x=d, distribution = object$distO, ties = object$ties, method = object$method, bind = object$bind)))

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
                                 useBoot=FALSE, ...) {
  stopifnot(inherits(object, 'incubate_fit'))
  stopifnot(is.numeric(level), length(level) == 1L, level < 1L, level > 0L)
  stopifnot(is.numeric(R), length(R) == 1L, R > 0)
  if (missing(bs_data)) bs_data <- 'parametric'
  if (is.vector(bs_data) && is.character(bs_data)) bs_data <- match.arg(bs_data[[1L]], choices = c('parametric', 'ordinary'))
  bs_infer <- match.arg(bs_infer)
  logTransform <- isTRUE(startsWith(bs_infer, 'log'))

  twoGroup <- isTRUE(object$twoGroup)
  nObs <- object$nobs

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
  stopifnot(!is.vector(bs_data) && !is.character(bs_data))
  # set R according to the provided bs_data (in particular important when both R & bs_data object are given)
  R <- if (useBoot) bs_data[['R']] else NCOL(bs_data)
  if (R < 999) warning(glue('Be cautious with the confidence interval(s) because the number of effective bootstrap samples R = {R} is rather low (R<999).'),
                       call. = FALSE)

  # logShift: needed only when log-transformation is requested. Start with a small standard value for all parameters
  logshift <- rlang::set_names(rep_len(.0001, length.out=length(pnames)), nm = pnames)
  # for delay, the transformation needs to be independent of the scale of delay, so we subtract the minimum and add a shift
  #+use fixed logshift_delay = 5 (which performed well in simulation at single group, exponential distribution, together with smd=0.25)
  if (logTransform) {
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

    stopifnot(is.matrix(bs_data))

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
