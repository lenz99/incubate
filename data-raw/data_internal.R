#!/usr/bin/env Rscript
# mkuhn, 2023-04-05
# simulate median weights W1, W2 and W3 for the weighed MLE approach (Cousineau, 2009)
####

# init -----

message("Start script at ", toString(Sys.time()))

library("usethis")
library("readr")
library("tibble")
library("tidyr")
library("dplyr")
library("purrr")
library("future")
library("future.callr")
library("furrr")
library("matrixStats", warn.conflicts = FALSE)

suppressPackageStartupMessages(library('R.utils'))

TODAY <- Sys.Date()

# command line arguments -----
cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(resultsDir = getwd(), seed=as.integer(TODAY),
                                                workers=future::availableCores(methods = "system", omit = 4), mcnrep="1001"))


if (any(c('help', 'h') %in% names(cmdArgs))){
  cat('Run Monte-Carlo simulations to estimate the median weights W1, W2 and W3 for weighted maximum likelihood approach (MLEw).\n')
  cat('See as reference Cousineau, 2009.\n')
  cat('  --help\t print this help\n')
  cat('  --resultsDir=\t specify the directory where to put the result files. Defaults to the directory where Rscript is executed.\n')
  cat('  --seed=\t if given, set random seed at the start of the script. Default is date-dependent.\n')
  cat('  --workers=\t number of parallel computations using `future.callr`. The only level of parallelization is for the different numbers of observations (and scale for W3).\n')
  cat('  --mcnrep=\t size of Monte-Carlo study: number of replications which are then aggregated.\n')
  quit(save = 'no')
}

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir),
           # check read & write permission (first octal information)
           (file.mode(myResultsDir) %>% as.character() %>% substr(1,1) %>% as.octmode() & 6) == '6')

mySeed <- cmdArgs[["seed"]]
stopifnot( is.numeric(mySeed), length(mySeed) == 1L, mySeed >= 0L )

myWorkers <- cmdArgs[["workers"]]
stopifnot( is.numeric(myWorkers), length(myWorkers) == 1L, myWorkers >= 1L )

myMCNrep <- readr::parse_number(cmdArgs[["mcnrep"]])
stopifnot( is.numeric(myMCNrep), length(myMCNrep) == 1L, myMCNrep >= 1L )



# set up simulation setting -----

if (mySeed > 0L) set.seed(mySeed)
if (myWorkers > 1L) {
  future::plan(strategy = future.callr::callr, workers = myWorkers)
}


# distribution of W1 is Gamma with shape n and scale 1/n
nObs <- c(1:20, 25, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000, 2500)
shape_W3 <- c(0.01, 0.05, 0.1, 0.25, 0.5, .75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7)

aggFun <- stats::median; isMedian <- TRUE
stopifnot( is.function(aggFun), "na.rm" %in% formalArgs(aggFun) )

# simulation W1 ------------------------

message("Start with W1")

W1_mc <- furrr::future_map_dbl(.x = nObs,
                               .f = ~ aggFun(rgamma(n=myMCNrep, shape = .x, scale = 1/.x)),
                               .options = furrr_options(seed = TRUE)) %>%
  purrr::set_names(nm = nObs)
if ( abs(W1_mc[[1]] - log(2)) > 1e-3 ) {
  warning("For n=1, W1 deviates more than 1e-3 from the true value ln(2)!", call. = FALSE)
}
W1_mc[[1L]] <- log(2)


# simulation of W2 -----------------------------------------

message("Start with W2")
W2_mc <- furrr::future_map_dbl(.x = nObs, .f = ~ aggFun(replicate(n = myMCNrep,
                                                                 expr = {
                                                                   z <- rexp(n=.x)
                                                                   sum(z * log(z))/sum(z) - mean(log(z))
                                                                 })), .options = furrr_options(seed = TRUE)) %>%
  purrr::set_names(nm = nObs)
stopifnot( abs(W2_mc[[1]]) < 1e-5 )
W2_mc[[1L]] <- 0


# simulation of W3 --------------------------------------------------------

message("Start with W3")
#currently, W3 is using simulation on log-transform, aggregates and then backtransform via exp.
#+this works for median but for instance not for mean!
stopifnot(isMedian)

W3_mc_df <- tidyr::expand_grid(nObs = nObs,
                               shape = shape_W3) %>%
  dplyr::mutate(W3 = furrr::future_map2_dbl(.x = nObs, .y = shape,
                                            .f = ~ exp(aggFun(replicate(n = myMCNrep,
                                                                        expr = {
                                                                          z <- rexp(n=.x)
                                                                          res <- NA_real_
                                                                          # W1_mc[as.character(.x)] * sum(z**(-1/.y))/sum(z**((.y-1)/.y)) # on original (=non-log) scale
                                                                          try(
                                                                            expr = res <- log(W1_mc[as.character(.x)]) +
                                                                              matrixStats::logSumExp(lx = -1/.y * log(z)) -
                                                                              matrixStats::logSumExp(lx = (.y-1)/.y * log(z)),
                                                                            silent = TRUE
                                                                          )
                                                                          res
                                                                        }),
                                                              na.rm = TRUE)),
                                            .options = furrr_options(seed = TRUE)))



# save results ------------------------------------------------------------

MLEw_weights <- list(
  W12 = dplyr::inner_join(
    x = tibble::enframe(W1_mc, name = "nObs", value = "W1"),
    y = tibble::enframe(W2_mc, name = "nObs", value = "W2"),
    by = join_by(nObs)
  ),
  W3 = W3_mc_df,
  MCSS_setting = list(seed = mySeed,
                      mcnrep = myMCNrep)
)

#saveRDS(res_mc, file = "MLEweights.rds")
usethis::use_data(MLEw_weights, internal = TRUE, overwrite = FALSE)

# exit --------------------------------------------------------------------

# tear-down
future::plan(future::sequential())

message("~~ Fine ~~")
message("Finished script at ", toString(Sys.time()))
