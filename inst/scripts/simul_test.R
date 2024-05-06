#!/usr/bin/env Rscript
# Evaluate test for differences in delayed exponential or Weibull setting
#
# Test delay parameter


# init -----

cat("\nMC-simulations for test for difference in delay parameters.\n")
cat("It is ***", toString(Sys.time()), "***\n")

library("incubate")
# minimal version check:
#+ 0.7.6 for GOF-Pvalues for restricted & unrestricted model: e.g. gof_mo0 (was gof_mo) and gof_mo1 (new)
#+ 0.9.8 for names for P-values have changed: boot => bootstrap, gof_mo0 => moran, etc
#+ 1.1.9.9000 script is developed as part of the incubate package (not separate as part of the MS)
#+ 1.1.9.9014 ties='density' as default now also for tests
#+ 1.1.9.9016 avoid attributes, use transform() for Pearson/AD GOF tests
#+ 1.3.0.9025: rename logrank P-values to logrank and logrank_pp (to avoid confusion with likelihood ratio (=LR) tests)
#+ 1.3.0.9037: allow profiling for MPSE and all MLE-methods, at least with single group..
#+ 1.3.0.9055: support random right-censoring in rexp_delayed() and rweib_delayed()
stopifnot(packageVersion("incubate") >= "1.3.0.9055")
cat('incubate package version: ', toString(packageVersion("incubate")), '\n')

library("dplyr", warn.conflicts = FALSE)
stopifnot(packageVersion("dplyr") > "1.0.10")
library("purrr")
library("tidyr", warn.conflicts = FALSE)
library("tibble")
suppressPackageStartupMessages(library("R.utils"))

TODAY <- Sys.Date()


# command line arguments -----
cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(
                                  # simulation settings
                                  dist='exponential', scenario='MS',
                                  R=150, mcnrep=100,
                                  # technical settings
                                  resultsDir = getwd(),
                                  slice=0, seed=as.integer(TODAY),
                                  chnkSize=0, workers=3))


if (any(c('help', 'h') %in% names(cmdArgs))) {
  cat('Run Monte-Carlo simulations with delayed exponential or Weibull data in a two group setting.\n')
  cat('A test for difference in delay (and sometimes delay+rate) is performed.\n')
  cat('Sample size, delay, scale and scale ratio (between the two groups) and shape use different fixed values (see code in this script).\n')
  cat('Command line parameter options allow to adjust what this script actually does:\n')
  cat('  --help\t print this help\n')
  cat('  --print\t show scenarios to simulate and exit.\n')
  cat('  --resultsDir=\t specify the directory where to put the result files. Defaults to the directory where Rscript is executed.\n')
  cat('  --dist=\t specify distribution that governs the data generation. Default is the exponential distribution.\n')
  cat('  --scenario=\t with respect to the delay in both groups, choose a scenario for the simulation:\n\t\t\tDELAYEQ = no difference in delay,\n\t\t\tDELAYGT = 2nd group y with bigger delay.\n\t\t\tMS = only relevant scenarios shown in manuscript (default)\n\t\t\tALL = all cases\n')
  cat('  --allN\t use different sample sizes in the simulations. Without this option, only a single sample size is used.\n')
  cat('  --scaleSimple\t use only standard value for scale and scale-ratio\n')
  cat('  --includeMLEw\t include also weighted MLE approach\n')
  cat('  --cens\t apply also random right-censoring during the simulation study\n')
  cat('  --slice=\t if given, pick only this number of first scenarios for simulations. If negative, scenarios taken from the tail.\n')
  cat('  --seed=\t if given, set random seed at the start of the script. Default is date-dependent.\n')
  cat('  --chnkSize=\t chunk size to write out results having processed so many scenarios. Default is no chunking (=0).\n')
  cat('  --workers=\t number of parallel computations using `future.callr` and `future.apply`. The only level of parallelization is across the MC-replications for each simulation setting.\n')
  cat('  --R=\t\t number of samples within parametric bootstrap test: it determines the resolution for our P-value, e.g.,\n\t\t R=100 will allow for P-values at per-cent resolution\n')
  cat('  --mcnrep=\t size of Monte-Carlo study: it is the number of replicated bootstrap data sets on which statistical tests are done.\n')
  quit(save = 'no')
}

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir),
           # check read & write permission (first octal information)
           (file.mode(myResultsDir) %>% as.character() %>% substr(1,1) %>% as.octmode() & 6) == '6')

myDist <- cmdArgs[["dist"]]
stopifnot(is.character(myDist), length(myDist) == 1L)
myDist <- match.arg(arg = tolower(myDist), choices = c("exponential", "weibull"))
isExpon <- isTRUE(myDist == "exponential")
stopifnot(isExpon || isTRUE(myDist == "weibull"))

myWorkers <- cmdArgs[["workers"]]
stopifnot(is.numeric(myWorkers), length(myWorkers) == 1L, myWorkers >= 1L)
USE_FUTURE <- myWorkers > 1L

myChnkSize <- cmdArgs[["chnkSize"]]
stopifnot( is.numeric(myChnkSize), length(myChnkSize) == 1L )

myR <- cmdArgs[["R"]]
stopifnot( is.numeric(myR), length(myR) == 1L, myR >= 1L )

myMCNrep <- cmdArgs[["mcnrep"]]
stopifnot(is.numeric(myMCNrep), length(myMCNrep) == 1L, myMCNrep >= 1L)

mySlice <- cmdArgs[["slice"]]
stopifnot(! is.null(mySlice), is.numeric(mySlice), length(mySlice) == 1L)
mySlice <- ceiling(mySlice)

mySeed <- cmdArgs[["seed"]]
stopifnot( is.numeric(mySeed), length(mySeed) == 1L, mySeed >= 0L )

myScenario <- cmdArgs[["scenario"]]
stopifnot(! is.null(myScenario), is.character(myScenario), length(myScenario) == 1L, nzchar(myScenario))
myScenario <- match.arg(arg = toupper(myScenario), choices = c("DELAYEQ", "DELAYGT", "MS", "ALL"))

myPrint <- isTRUE(any(c("print", "p") %in% tolower(names(cmdArgs))))
myAllN <- isTRUE(any(c("alln", "a") %in% tolower(names(cmdArgs))))
myIncludeMLEw <- isTRUE(any("includemlew" %in% tolower(names(cmdArgs))))
myScaleSimple <- isTRUE(any(c("scalesimple", "scale", "scales") %in% tolower(names(cmdArgs))))
myCens <- isTRUE(any(c("cens", "censoring") %in% tolower(names(cmdArgs))))




# set up simulation setting -----

if (mySeed > 0L) set.seed(mySeed)

simSetting <- tidyr::expand_grid(n_x = c(8, 10, 12, 15, 20, 30, 50, 75), #100
                                 delay_x = 5,
                                 delay_y = c(5, 7, 9, 11, 13, 15), #, 20),
                                 scale_x = c(5, 10), #c(1, 2, 5),
                                 scale_ratio = c(2, 1, .5),
                                 # shape values according to distribution
                                 #+effectively filter for distribution
                                 shape = if (isExpon) 1 else c(.5, 2),
                                 cens = c(0, 0.1, 0.2, 0.3))

# avoid duplicates:
# by convention, group y is not less delayed than group x
simSetting <- simSetting %>%
  # symmetry
  dplyr::filter(delay_y >= delay_x) %>%
  # enough expected number of observations
  dplyr::filter(cens >= 0, cens < 1, n_x * (1-cens) > 5) %>%
  # use equally sized groups
  dplyr::mutate(n_y = n_x, .after = n_x)

# default: no censoring (cens = 0)
if (!myCens) {
  simSetting <- simSetting %>%
    dplyr::slice_min(cens)
}

# default is to use only the largest sample size
if (!myAllN) {
  simSetting <- simSetting %>%
    dplyr::slice_max(n_x)
}

if (myScaleSimple) {
  simSetting <- simSetting %>%
    dplyr::filter(dplyr::near(scale_x, 10),
                  dplyr::near(scale_ratio, 1))
}


# filter for target scenario!
# use all capital letters (see definition of myScenario)
simSetting <- switch (myScenario,
                      DELAYEQ = {
                        simSetting %>%
                          dplyr::filter(dplyr::near(delay_x, delay_y))
                        # # for both distributions:
                        # dplyr::near(scale_x, 10), dplyr::near(scale_ratio, 1))
                      },

                      DELAYGT = {
                        simSetting %>%
                          dplyr::filter(delay_y > delay_x + 1e-11)
                      },

                      MS = {

                        #filter only relevant scale_ratio combinations
                        # for exponential:
                        # [G1] scale_x = 10 & scale_ratio = 1  (rate_x = .1, rate_ratio = 1)
                        # [G2] scale_x =  5 & scale_ratio = 1  (rate_x = .2, rate_ratio = 1)
                        # [G3] scale_x =  5 & scale_ratio =  2 [i.e., scale_y = 10]  (rate_x = .2, rate_ratio = .5 [i.e., rate_y = .1])
                        # [G4] scale_x = 10 & scale_ratio = .5 [i.e., scale_y =  5]  (rate_x = .1, rate_ratio = 2  [i.e., rate_y = .2])
                        # for weibull:
                        # dd=0, k=.5|2, scale_x = 10 & scale_ratio = 1
                        # dd=5, k=.5|2, scale_x = 10 & scale_ratio = 1

                        # filter based on scale parameters
                        # this contains the cases which are needed in the manuscript
                        local({
                          simFilterMS <- if (isExpon) {
                            simSetting %>%
                              dplyr::filter(dplyr::near(scale_ratio, 1) |
                                              (dplyr::near(scale_x, 5) & dplyr::near(scale_ratio, 2)) |
                                              (dplyr::near(scale_x, 10) & dplyr::near(scale_ratio, .5)))
                          } else {
                            tibble::tibble(scale_x = 10, scale_ratio = 1)
                          }

                          simSetting %>%
                            dplyr::inner_join(simFilterMS, by = colnames(simFilterMS))
                        })
                      },

                      # ALL = no-op, keep everything! (also unused border cases)
                      ALL = {
                        simSetting
                      },
                      stop("Unknown target!", call. = FALSE)
)



if (!dplyr::near(mySlice, 0)) {
  simSetting <- local({

    sliceF <- if (mySlice > 0) dplyr::slice_head else dplyr::slice_tail

    simSetting %>%
      sliceF(n = abs(mySlice))
  })
}

if (myPrint) {
  print(knitr::kable(simSetting, format = 'pipe', digits = 2))
  cat('\n')
  cat(NROW(simSetting), 'simulation scenarios in total.\n')
  cat('Each scenario is covered by ', myMCNrep, 'MC-data replications.\n')
  cat('Bootstrap tests with R=', myR, 'parametric bootstrap samples (P-value resolution).\n')
  cat('Seed set initially is: ', if (mySeed>0) mySeed else '-not set-', '\n')
  cat('Results directory is set to ', myResultsDir, '\n')

  quit(save = 'no')
}

# set up parallel computing ----
if (USE_FUTURE) {
  library("future.callr")
  library("future.apply")

  future::plan(strategy = future.callr::callr, workers = myWorkers)
  # two level future
  # future::plan(list(
  #   tweak(future.callr::callr, workers = 2L),
  #   tweak(multicore, workers = 4L)
  # ))
}



# functions -----

#' Run Monte-Carlo simulations to test difference in delay using an exponential model for a given simulation setting
#'
#' A fixed set of estimation methods are used.
#' Uses parallel computation (future_replicate) to go through the (=nrep) MC-simulations.
#' Each bootstrap test is also future-aware (and would pick up a nested future-plan setting)
#' @param DGPsetting numeric. a row from `simSetting`. It encodes parameters that specify the data generating process for both groups
#' @param include_MLEw logical. Should we also include MLEw?
#' @return dataframe. P-values in the different Monte-Carlo runs.
doMCSim <- function(DGPsetting, include_MLEw = TRUE) {
  # settings from the environment:
  stopifnot(exists("isExpon"), exists("myMCNrep"), exists("myR"))
  stopifnot(is.numeric(DGPsetting), length(DGPsetting) == 8L)
  stopifnot(is.logical(include_MLEw), length(include_MLEw) == 1L)
  include_MLEw <- isTRUE(include_MLEw)

  n_x <- DGPsetting[[1]]
  n_y <- DGPsetting[[2]]
  delay_x <- DGPsetting[[3]]
  delay_y <- DGPsetting[[4]]

  scale_x <- DGPsetting[[5]]
  scale_ratio <- DGPsetting[[6]]
  shape <- DGPsetting[[7]]
  cens <- DGPsetting[[8]]

  # do we test for parameters combined?
  testParamCombined <- scale_ratio != 1
  scale_y <- scale_x * scale_ratio

  # different estimation methods
  estimMethods <- tidyr::expand_grid(method = c(c("MPSE", "MLEn", "MLEc"), if (include_MLEw) "MLEw"),
                                     profiled = c(FALSE, TRUE),
                                     R = as.integer(myR)) %>%
    # all MLE-methods use only profiled variant, MPSE uses both, profiled & unprofiled
    dplyr::filter(method == 'MPSE' | profiled) %>%
    dplyr::rowwise()

  testDiffList <- future.apply::future_replicate(n = myMCNrep,
                                                 future.packages = c("dplyr", "incubate", if (cens > 0) "survival"),
                                                 future.seed = TRUE,
                                                 expr = {
                                                   # generate data
                                                   x <- y <- 1 #dummy init
                                                   if (isExpon) {
                                                     stopifnot(dplyr::near(shape, 1L))
                                                     x <- rexp_delayed(n = n_x, delay1 = delay_x, rate1 = 1/scale_x, cens = cens)
                                                     y <- rexp_delayed(n = n_y, delay1 = delay_y, rate1 = 1/scale_y, cens = cens)
                                                   } else {
                                                     # weibull
                                                     x <- rweib_delayed(n = n_x, delay1 = delay_x, scale1 = scale_x, shape1 = shape, cens = cens)
                                                     y <- rweib_delayed(n = n_y, delay1 = delay_y, scale1 = scale_y, shape1 = shape, cens = cens)
                                                   }


                                                   estimMethods %>%
                                                     dplyr::mutate(testDiffObj = list({
                                                       te_diff <- NULL
                                                       # test_diff might also use parallel computations depending on future-settings
                                                       try(expr = {
                                                         # test difference in delay1 in exponential model
                                                         te_diff <- test_diff(x = x, y = y, distribution = "expon", param = "delay1",
                                                                              method = method, profiled = profiled, R = R, type = "all",
                                                                              # log-rank test only once
                                                                              doLogrank = method == "MPSE" && !profiled)
                                                         # bootstrap P-value for combined test for difference in parameters delay+rate
                                                         #+only if the scale (=1/rate for exponential) is indeed different betw groups
                                                         if (testParamCombined) {
                                                           te_diff2 <- test_diff(x = x, y = y, distribution = "expon", param = c("delay1", "rate1"),
                                                                                 method = method, profiled = profiled, R = R, type = "bootstrap")
                                                           # store P-value of delay+rate in original test_diff-object
                                                           te_diff$P$bootstrap2 <- purrr::pluck(te_diff2, "P", "bootstrap", .default = NA_real_)
                                                         }#fi
                                                       }, silent = TRUE)

                                                       te_diff })) %>%
                                                     # compact testDiff-list column: drop entries that did not work out!
                                                     dplyr::filter(!is.null(testDiffObj)) %>%
                                                     # extract all P-values/R_eff in long format from each row in estimMethods-df!
                                                     #+dplyr::reframe (beta in v1.1.0) allows to summarize with more than one row
                                                     dplyr::reframe(method, profiled, R,
                                                                    R_eff = length(testDiffObj$testDist),
                                                                    tibble::enframe(unlist(testDiffObj$P),
                                                                                    name = "test", value = "pvalue"))
                                                 }, simplify = FALSE)

  # drop NULLs (just in case)
  testDiffList <- purrr::compact(testDiffList)

  # bind together into a single long tibble
  dplyr::bind_rows(testDiffList, .id = "run")
} #fn doMCSim


#' Run MC-simulations for each scenario sequentially (row-by-row)
#' @param simSetDF dataframe containing simulation scenarios
#' @param ... further arguments passed to `doMCSim` (currently not used!)
#' @returns tibble of simulations settings with results added
applyMCSims <- function(simSetDF, ...) {
  simSetDF %>%
    dplyr::mutate(., results = apply(as.matrix(.), MARGIN = 1L, FUN = doMCSim, includeMLEw = myIncludeMLEw, ...))
}



#' Add meta data to dataframe
#'
#' The special `comment` attribute is used to store the deparsed meta data list.
#' @param da simulation data
#' @param timeTag time stamp to be added in meta data comment
#' @returns simulation data with meta data added as comment
addMetaData <- function(da, timeTag) {
  # add comment as text
  comment(da) <- list(seed = mySeed, R = myR, mcnrep = myMCNrep, workers = myWorkers, chnkSize = myChnkSize,
                      host = Sys.info()[["nodename"]],
                      rversion = R.version.string,
                      incubate = as.character(packageVersion("incubate")),
                      date = TODAY,
                      time = timeTag) %>%
    #paste(names(.), ., sep = '=', collapse = ',')
    deparse()

  da
}



# run & save ----

DATETIME_TAG <- format(Sys.time(), format = "%Y-%m-%d-%Hh%Mm%Ss")
rdsBaseName <- paste0("simRes_test_", DATETIME_TAG)
rdsName <- file.path(myResultsDir, paste0(rdsBaseName, ".rds"))

if (myChnkSize < 1L || NROW(simSetting) <= myChnkSize) {
  # no chunking
  simSetting <- applyMCSims(simSetDF = simSetting) %>%
    addMetaData(timeTag = DATETIME_TAG)

  saveRDS(simSetting, file = rdsName)

} else {

  # work in chunks
  rowIdx <- seq_len(NROW(simSetting))
  # how many chunks?
  chnkNbr <- (length(rowIdx) %/% myChnkSize)+1L
  stopifnot(chnkNbr > 1L, chnkNbr <= 999999L)
  # stripe over the scenarios
  rowIdxLst <- split(rowIdx, f = rep_len(x=seq_len(chnkNbr), length.out = length(rowIdx)))
  stopifnot(length(rowIdxLst) == chnkNbr)

  for (i in seq_along(rowIdxLst)) {
    simSetting_chnk <- dplyr::slice(simSetting, rowIdxLst[[i]]) %>%
      applyMCSims() %>%
      addMetaData(timeTag = DATETIME_TAG)
    #simSetting_chnk <- applyMCSims(simSetDF = simSetting_chnk)

    message("Writing out chunk ", i, " to RDS-file..")
    saveRDS(simSetting_chnk,
            file = file.path(myResultsDir, paste0(rdsBaseName, "_", sprintf("%06d", i), ".rds")))
  }#rof

  # merge chunked output!
  chnkFileNames <- list.files(path = myResultsDir,
                              pattern = paste0('^', rdsBaseName, '_[[:digit:]]+[.]rds$'),
                              full.names = TRUE)
  if (length(chnkFileNames)) {
    # re-create complete simSetting data
    simSetting <- purrr::map(.x = chnkFileNames, .f = readRDS) %>%
      dplyr::bind_rows() %>%
      addMetaData(timeTag = DATETIME_TAG)

    saveRDS(simSetting, file = rdsName)

    if (file.exists(rdsName) && (! exists('infoRDS') || ! inherits( try(infoRDS(rdsName), silent = TRUE), "try-error"))) {
      message("Removing ", length(chnkFileNames), " intermediate chunked RDS-files!")
      try(file.remove(chnkFileNames))
    } #fi remove RDS-chunk-files

  } else {
    warning("Did not find chunked RDS-output.", call. = FALSE)
  }

} #esle chunking



# example for a visualization of test results
# simSetting %>%
#   filter(dplyr::near(n, 10), dplyr::near(delay_x, 5), dplyr::near(rate_x, .1)) %>%
#   unnest(cols = P) %>%
#   ggplot(mapping = aes(x = P, col = method)) +
#   geom_freqpoly(bins = 12) + xlim(0,1) +
#   coord_trans(y = "log1p") +
#   facet_grid(rows = vars(delay_y), cols = vars(rate_ratio), labeller = label_both) +
#   labs(x = "P-value", title = "Test Results under different group effects", subtitle = "**n = 10**, delay~x~ = 5")



# teardown ----

# output the latest warnings:
cat("\n+++\nThese are warnings from the script:\n+++\n")
warnings()

if (USE_FUTURE && isNamespaceLoaded("future")) {
  future::plan(strategy = future::sequential)
}

cat("It is ***", toString(Sys.time()), "***\n")
cat("\n\n~fine~\n")

