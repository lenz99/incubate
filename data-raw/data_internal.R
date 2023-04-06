#!/usr/bin/env Rscript
# mkuhn, 2023-04-05
# internal data for the incubate package
#
# simulate median weights W1, W2 and W3 for the weighed MLE approach (Cousineau, 2009)
####

# init -----

message("Start script for internal data at ", toString(Sys.time()))

library("rlang")
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
  cat('  --seed=\t if given, set random seed at the start of the script. Default is date-dependent.\n')
  cat('  --workers=\t number of parallel computations using `future.callr`. The only level of parallelization is for the different numbers of observations (and scale for W3).\n')
  cat('  --mcnrep=\t size of Monte-Carlo study: number of replications which are then aggregated.\n')
  cat('  --resultsDir=\t directory where to save the result files (when not internal) Defaults to the directory where Rscript is executed.\n')
  cat('  --internal/-i\t save as internal package data. Then `resultsDir` is irrelevant.\n')
  cat('  --overwrite/-f\t Set `overwrite=TRUE` when saving data.\n')
  quit(save = 'no')
}


mySeed <- cmdArgs[["seed"]]
stopifnot( is.numeric(mySeed), length(mySeed) == 1L, mySeed >= 0L )

myWorkers <- cmdArgs[["workers"]]
stopifnot( is.numeric(myWorkers), length(myWorkers) == 1L, myWorkers >= 1L )

myMCNrep <- readr::parse_number(cmdArgs[["mcnrep"]])
stopifnot( is.numeric(myMCNrep), length(myMCNrep) == 1L, myMCNrep >= 1L )

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir),
           # check read & write permission (first octal information)
           (file.mode(myResultsDir) %>% as.character() %>% substr(1,1) %>% as.octmode() & 6) == '6')
myInternal <- isTRUE(any(c("internal", "i") %in% tolower(names(cmdArgs))))
myOverwrite <- isTRUE(any(c("overwrite", "ow", "f") %in% tolower(names(cmdArgs))))


# set up simulation setting -----

if (mySeed > 0L) set.seed(mySeed)
if (myWorkers > 1L) {
  future::plan(strategy = future.callr::callr, workers = myWorkers)
}


# distribution of W1 is Gamma with shape n and scale 1/n
nObs <- c(1:25, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000, 2500, 5000, 10000)
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
  warning("For n=1, W1 deviates more than 1e-3 from the true value ln(2)! (We use ln(2) instead, anyhow.)", call. = FALSE)
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

W3_mc_df <- tidyr::expand_grid(nObs = as.integer(nObs),
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


W12 <- dplyr::inner_join(
  x = tibble::enframe(W1_mc, name = "nObs", value = "W1"),
  y = tibble::enframe(W2_mc, name = "nObs", value = "W2"),
  by = join_by(nObs) ) %>%
  dplyr::mutate(nObs = as.integer(nObs))


# approximations ----------------------------------------------------------

if (rlang::is_interactive()) {
  library("ggplot2")
  library("patchwork")

  W12 <- W12 %>%
    dplyr::mutate(lnObs = log(nObs),
                  # approximation for median of gamma(n, 1/n)
                  #+using Wilson-Hilferty transformation (see <https://en.wikipedia.org/wiki/Gamma_distribution>)
                  W1pred = (1 - 1 / (9 * nObs))^3,
                  W1diff = W1 - W1pred)

  ggplot(W12, mapping = aes(x = nObs, y = W1)) +
    geom_point() +
    scale_x_log10() +
    labs(title = "W1 as median") +
    geom_line(mapping = aes(y = W1pred), col = "blue")  |

    ggplot(W12, mapping = aes(x = nObs, y = W1diff)) +
    geom_point() + geom_line() +
    scale_x_log10() +
    labs(title = "Deviation")

  # (0,0) point is taken out
  # SSasymp on log(nObs):
  # W2 = 1 + (R0 - 1) * n**-r
  fm_W2_A <- nls(W2 ~ SSasymp(input = lnObs, Asym = 1, R0, lrc), start = list(R0 = 0, lrc = -.02), data = W12, subset = -1)
  fm_W2_A

  W12 <- W12 %>%
    dplyr::mutate(W2pred = c(NA_real_, predict(fm_W2_A)),
                  W2diff = W2 - W2pred)

  W2_R0 <- coef(fm_W2_A)[["R0"]]
  W2_nr <- -exp(coef(fm_W2_A)[["lrc"]])

  ggplot(W12, mapping = aes(x = nObs, y = W2)) +
    geom_point() +
    scale_x_log10() +
    labs(title = "W2 as median", subtitle = "Point (0,0) does not follow the pattern") +
    ylim(0, NA) +
    geom_line(mapping = aes(y = W2pred), col = "blue") +
    # approximating function (same as predict)
    geom_function(fun = ~ 1 + (W2_R0 - 1) * .x**W2_nr, col = "red") |

    ggplot(W12, mapping = aes(x = nObs, y = W2diff)) +
    geom_point() + geom_line() +
    scale_x_log10()

  W3 <- W3_mc_df %>% dplyr::mutate(lnObs = log(nObs),
                             shapeF = factor(shape))

  W3_smallShape <- W3 %>% filter(shape < 1) %>% droplevels()
  W3_bigShape <- W3 %>% filter(shape >= 1) %>% droplevels()

  ggplot(W3_bigShape, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
    geom_point() + geom_line() +
    scale_x_log10() +
    scale_y_log10() |

    ggplot(W3_smallShape, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
    geom_point() + geom_line() +
    scale_x_log10() +
    scale_y_log10()

  fm_W3_A <- nls(W3 ~ SSasymp(input = lnObs, Asym, R0 = 1, lrc), start = list(Asym = 10, lrc = 0),
                 data = W3_bigShape, subset = shapeF == "1.25")
  summary(fm_W3_A)

  fm_W3_B <- nls(W3 ~ SSasymp(input = lnObs, Asym[shapeF], R0 = 1, lrc[shapeF]),
                 start = list(Asym = c(15, 5, 3, rep.int(2.3815, 12)), lrc = rep.int(-1.32, 15)),
                 data = W3_bigShape)

  fm_W3_C <- lm(log(W3) ~ nObs + (lnObs + sqrt(lnObs) + log(lnObs+1)) * shapeF, data = W3_bigShape)
  summary(fm_W3_C)
  drop1(fm_W3_C, test = "F")
}

# save results ------------------------------------------------------------

MLEw_weights <- list(
  W12 = W12,
  W3 = W3_mc_df,
  MCSS_setting = list(seed = mySeed,
                      aggFun = aggFun,
                      mcnrep = myMCNrep)
)




if (myInternal) {
  if (!inherits(try(expr = usethis::proj_get(), silent = TRUE), what = "try-error")) {
    message("Save weights for MLEw as internal package data.")
    usethis::use_data(MLEw_weights, internal = TRUE, overwrite = myOverwrite)
  } else warning("Unable to save internal package data because there is no active package project! Please run from within a package project..")
} else {
  message("Save weights for MLEw as RDS file.")
  rdsFile <- file.path(myResultsDir, "MLEw_weights.rds")
  if (file.exists(rdsFile) && ! myOverwrite) {
    warning("File ", rdsFile, "already exists! You would need to set overwrite-flag.")
  } else {
    saveRDS(MLEw_weights, file = rdsFile)
  }
}


# exit --------------------------------------------------------------------

# tear-down
future::plan(future::sequential())

message("~~ Fine ~~")
message("Finished script at ", toString(Sys.time()))
