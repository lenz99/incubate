#!/usr/bin/env Rscript
# mkuhn, 2023-04-05
# internal data for the incubate package
#
# simulate median weights W1, W2 and W3 for the weighed MLE approach (Cousineau, 2009)
####

# init -----

message("Script to prepare MLE weights to be stored as internal data of incubate package!")
message("incubate package installed is: ", packageVersion("incubate"))
message("Start at ", toString(Sys.time()))

library("rlang")
library("usethis")
library("readr") #parse_number
library("tibble")
library("tidyr", warn.conflicts = FALSE)
library("dplyr", warn.conflicts = FALSE)
library("purrr", warn.conflicts = FALSE)
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot(font_size = 15))
library("patchwork")

library("gslnls")
library("splines")
library("matrixStats", warn.conflicts = FALSE)

library("future")
library("future.callr")
library("furrr")


suppressPackageStartupMessages(library('R.utils'))

TODAY <- Sys.Date()
DEBUG <- FALSE

# command line arguments -----
cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(resultsDir = getwd(), seed=as.integer(TODAY),
                                                # at most 97 cores
                                                workers = min(97L, future::availableCores(methods = "system", omit = 4)),
                                                # mcnrep as string, so user can rely on parse_number!
                                                mcnrep="1001"))


if (any(c('help', 'h') %in% names(cmdArgs))) {
  cat('Run Monte-Carlo simulations to estimate the median weights W1, W2 and W3 for weighted maximum likelihood approach (MLEw)\n')
  cat('And also find approximating functions for these weights.\n')
  cat('See as reference Cousineau, 2009.\n')
  cat('  --help\t print this help\n')
  cat('  --seed=\t if given, set random seed at the start of the script. Default is date-dependent.\n')
  cat('  --workers=\t number of parallel computations using `future.callr`. The only level of parallelization is for n, the different numbers of observations (and scale for W3).\n')
  cat('  --mcnrep=\t size of Monte-Carlo study: number of replications which are then aggregated. Default value is 1001.\n')
  cat('  --resultsDir=\t directory where to save the result files (when not internal) Defaults to the directory where Rscript is executed.\n')
  cat('  --overwrite/--force\t Overwrite data file when it already exists?\n')
  quit(save = 'no')
}


mySeed <- cmdArgs[["seed"]]
stopifnot(is.numeric(mySeed), length(mySeed) == 1L, mySeed >= 0L)

myWorkers <- cmdArgs[["workers"]]
stopifnot(is.numeric(myWorkers), length(myWorkers) == 1L, is.finite(myWorkers), myWorkers >= 1L)

myMCNrep <- readr::parse_number(cmdArgs[["mcnrep"]])
stopifnot(is.numeric(myMCNrep), length(myMCNrep) == 1L, myMCNrep >= 1L)

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot(is.character(myResultsDir), dir.exists(myResultsDir),
           # check read & write permission (first octal information)
           (file.mode(myResultsDir) %>% as.character() %>% substr(1,1) %>% as.octmode() & 6) == '6')
myOverwrite <- isTRUE(any(c("overwrite", "ow", "force") %in% tolower(names(cmdArgs))))

if (DEBUG) {
  cat(paste(names(cmdArgs), cmdArgs, sep = ": ", collapse = "***"), "\n")
  cat("Overwrite: ", myOverwrite, "\n")
}

# fail early
rdataFile <- file.path(myResultsDir, "MLEw_weights.RData")
if (file.exists(rdataFile) && !myOverwrite) {
  stop("File ", rdataFile, "already exists! You would need to set overwrite-flag.")
}




# set up simulation settings -----

if (mySeed > 0L) set.seed(mySeed)
if (myWorkers > 1L) {
  future::plan(strategy = future.callr::callr, workers = myWorkers)
}


# distribution of W1 is Gamma with shape n and scale 1/n
nObs_vctr <- c(1:50, 55, 60, 65, 70, 75, 80, 90, 100, 125, 150, 200, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 7500, 10000) |>
  unique()
#shape relevant for W3
shape_vctr <- c(0.01, 0.05, 0.1, 0.25, 0.5, .75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7,
                8:15, 20, 25, 30, 40, 50, 75, 100) |>
  unique()

aggFun <- stats::median; isMedian <- TRUE
stopifnot(is.function(aggFun), "na.rm" %in% formalArgs(aggFun))



# simulate W1 ------------------------

message("Start simulation for W1")

# when 3-param Weibull holds then the mean of z values (where z is Exp(1)) are gamma-distributed with parameter shape n and scale 1/n
# hence, the mean of W1 is 1 (independently of n)
W1_mcs <- furrr::future_map_dbl(.x = rlang::set_names(nObs_vctr),
                                .f = ~ aggFun(stats::rgamma(n=myMCNrep, shape = .x, scale = 1/.x)),
                                .options = furrr_options(seed = TRUE))

if (nObs_vctr[[1L]] == 1) {
  if (!dplyr::near(W1_mcs[[1]], log(2), tol = 1e-3)) {
    warning("For n=1, Monte Carlo simulation for W1 deviates more than 1e-3 from the true value ln(2)! (We use ln(2) instead, anyhow.)", call. = FALSE)
  }
  W1_mcs[[1L]] <- log(2)
} else {
  stop("n=1 has not run!!")
}


# simulate W2 -----------------------------------------

message("Start simulation for W2")
W2_mcs <- furrr::future_map_dbl(.x = rlang::set_names(nObs_vctr),
                                .f = ~ aggFun(replicate(n = myMCNrep,
                                                        expr = {
                                                          z <- stats::rexp(n=.x)
                                                          sum(z * log(z))/sum(z) - mean(log(z))
                                                        })),
                                .options = furrr_options(seed = TRUE))

if (nObs_vctr[[1L]] == 1) {
  if (abs(W2_mcs[[1]]) > 1e-5) {
    warning("W2 Monte Carlo simulation for n=1 not close to zer0!. We use 0!")
  }
  W2_mcs[[1L]] <- 0
} else {
  stop("n=1 has not run!!")
}

# gather results for W1 & W2 in a combined dataframe
W12_mcs_df <- dplyr::inner_join(
  x = tibble::enframe(W1_mcs, name = "nObs", value = "W1"),
  y = tibble::enframe(W2_mcs, name = "nObs", value = "W2"),
  by = join_by(nObs)) |>
  dplyr::mutate(nObs = as.integer(nObs))



# simulate W3 --------------------------------------------------------

message("Start simulation for W3")
#currently, W3 is using simulation on log-transform, aggregates and then backtransform via exp.
#+this works for median but for instance not for mean!
stopifnot(isMedian)

W3_mcs_df <- tidyr::expand_grid(nObs = as.integer(nObs_vctr),
                                shape = shape_vctr) |>
  dplyr::mutate(W3 = furrr::future_map2_dbl(.x = nObs, .y = shape,
                                            .f = ~ exp(aggFun(replicate(n = myMCNrep,
                                                                        expr = {
                                                                          z <- stats::rexp(n=.x)
                                                                          res <- NA_real_
                                                                          # on original (=non-log) scale
                                                                          # W1_mcs[as.character(.x)] * sum(z**(-1/.y))/sum(z**((.y-1)/.y))
                                                                          try (
                                                                            expr = res <- log(W1_mcs[as.character(.x)]) +
                                                                              matrixStats::logSumExp(lx = -1/.y * log(z)) -
                                                                              matrixStats::logSumExp(lx = (.y-1)/.y * log(z)),
                                                                            silent = TRUE)
                                                                          res
                                                                        }), na.rm = TRUE)),
                                            .options = furrr_options(seed = TRUE)))


# Monte-Carlo simulation results for median
.MLEw_mcs <- list(
  W12 = W12_mcs_df,
  W3 = W3_mcs_df,
  settings = list(date = TODAY,
                  host = Sys.info()[["nodename"]],
                  R.version = R.version.string,
                  incubate = paste("installed: ", utils::packageVersion("incubate")),
                  seed = mySeed,
                  aggFun = aggFun,
                  mcnrep = myMCNrep)
)

try(expr = rm(W12_mcs_df, W3_mcs_df), silent = FALSE)

# save in between (in case approximation fails/stops for some reason)
saveRDS(.MLEw_mcs, file = file.path(myResultsDir, "MLEw_mcs.rds"))


#XXX approximation could go into data_internal.R (together with checks that the approximation is ok)

# approximation W1 ----------------------------------------------------------
message("Start with building approximations for the weights!")


# MLEw-weights W1
# W1 for given sample sizes (of one group)
# vectorized variant
w1F <- function(nObs) {

  if (missing(nObs) || !is.numeric(nObs) || any(!is.finite(nObs))) {
    stop("Please provide the number of observations within a group!", call. = FALSE)
  }

  nObs <- pmax.int(1L, nObs)

  nObsIdx_direct <- which(nObs <= 17)
  nObsIdx_approx <- which(nObs > 17)

  # first rows of W12 are complete and ordered
  stopifnot(identical(1:17, .MLEw_mcs$W12$nObs[1:17]))

  retV <- numeric(length(nObs))
  retV[nObsIdx_direct] <- .MLEw_mcs[["W12"]]$W1[nObs[nObsIdx_direct]]
  retV[nObsIdx_approx] <- (1 - 1 / (9 * nObs[nObsIdx_approx]))^3

  retV

}#fn w1F


# approximation W2 --------------------------------------------------------

# asymptotic regression model: cf. SSasymp model on log(nObs)
# starting at nObs = 2. nObs = 1 is off.
fm_W2 <- gsl_nls(W2 ~ 1 + (R0 - 1) * nObs**-exp(lr),
                 start = list(R0 = -.25, lr = -.01),
                 data = .MLEw_mcs$W12, subset = nObs > 1)

# check model fit
if (!fm_W2$convInfo$isConv || fm_W2$convInfo$stopCode != 0 || fm_W2$convInfo$nEval[["f"]] > 27 || deviance(fm_W2) > 0.01) {
  stop("Model fit for W2 is bad!")
}

# save infos for MLEw-approximation
.MLEw_approx <- list(
  coef = list(
    W2_R0 = coef(fm_W2)[[1L]], #was -0.44193638
    W2_negRate = -exp(coef(fm_W2)[[2L]]) #was -exp(-0.00624712316)
  )
)


# MLEw weight W2
# For given sample size of one group
# vectorized variant
w2F <- function(nObs) {

  if (missing(nObs) || !is.numeric(nObs) || any(!is.finite(nObs))) {
    stop("Please provide the number of observations within group!", call. = FALSE)
  }

  nObs <- pmax.int(1L, nObs)

  nObsIdx_direct <- which(nObs <= 17)
  nObsIdx_approx <- which(nObs > 17)

  # first rows of W12 are complete and ordered
  stopifnot(identical(1:17, .MLEw_mcs$W12$nObs[1:17]))

  retV <- numeric(length(nObs))
  retV[nObsIdx_direct] <- .MLEw_mcs[["W12"]]$W2[nObs[nObsIdx_direct]]
  # median approximation via asymptotic regression model SSasymp on log(n):
  # We hence model: W2 = 1 + (R0 - 1) * nObs**(-r)
  retV[nObsIdx_approx] <- 1 + (.MLEw_approx[["coef"]][["W2_R0"]]- 1) * nObs[nObsIdx_approx]**.MLEw_approx[["coef"]][["W2_negRate"]]

  retV
}#fn w2F


if (rlang::is_interactive()) {
  plDatW2 <- .MLEw_mcs$W12 |>
    dplyr::filter(nObs > 1) |>
    dplyr::mutate(W2mod = predict(fm_W2),
                  W2approx = w2F(nObs))

  ggplot(data = plDatW2,
         mapping = aes(x = nObs, y = W2)) +
    geom_point(alpha = .7) + #geom_line() +
    geom_point(mapping = aes(y = W2mod), size = .5, col = "darkred", alpha = .4) +
    geom_line(mapping = aes(y = W2mod), col = "darkred", alpha = .3, linetype = "dotted") +
    geom_line(mapping = aes(y = W2approx), col = "darkblue", alpha = .4, linetype = "dotdash") +
    scale_x_log10() +
    labs(title = "W2: median of MC-sim and fitted W2-function",
         subtitle = "Modell: red, Approx-fun: blue") |


  ggplot(data = plDatW2,
         mapping = aes(x = W2, y = W2-W2approx)) +
    geom_point() +
    labs(title = "Agreement of w2-approx and w2-simulation (med)")
} #fi




# approximation W3 --------------------------------------------------------

.MLEw_approx[["coef"]][["W3_richards"]] <- local({
  W3 <- .MLEw_mcs$W3 |>
    dplyr::filter(nObs > 1) |>
    dplyr::mutate(lshape = log(shape))

  ITER_MAX <- 1011

  # for each sample size nObs we fit
  # Richards' generalized logistic function (as function of shape)
  # This way we can estimate W3 for each shape value for those nObs
  #+that we have simulated in .MLEw_mcs$W3
  fm_W3_indiv <- purrr::map(.x = rlang::set_names(unique(W3$nObs)),
                            .f = function(.x) {
                              #W3 ~ A + (K - A) / (1 + Q / shape**B)**(1/nu) same formula, but little worse convergence
                              gsl_nls(W3 ~ A + (K - A) / (1 + Q * exp(-B * lshape))**(1/nu),
                                      data = W3, subset = nObs == {.x}, #jac = TRUE,
                                      start = list(A = 2*.x, K = 1, Q = .15 - .015 * log(.x), B = log(.x+1), nu = log(.x+1)/3),
                                      control = gsl_nls_control(maxiter = ITER_MAX))
                              })

  # check convergence for each model
  stopifnot(all(purrr::map_lgl(fm_W3_indiv, .f = list("convInfo", "isConv"))))
  stopifnot(all(purrr::map_dbl(fm_W3_indiv, .f = list("convInfo", "nEval", "f")) < ITER_MAX))


  if (rlang::is_interactive()) {

    # W3 approaches 1 for high shape values
    #+but W3pred overestimates these W3-values close to 1 (at high shape)
    #+the more the higher the sample sizes
    local({
      # pick some nObs
      myNObs <- intersect(55, unique(W3$nObs))
      stopifnot(is.finite(myNObs))

      plDatW3i <- W3 |>
        dplyr::filter(nObs == {myNObs}) |>
        dplyr::mutate(W3pred = predict(fm_W3_indiv[[as.character(myNObs)]]))

      ggplot(data = plDatW3i,
             mapping = aes(x = shape, y = W3)) +
        geom_point() + #geom_line() +
        geom_point(mapping = aes(y = W3pred), size = .5, col = "darkred") +
        geom_line(mapping = aes(y = W3pred), col = "darkred") +
        scale_x_log10() + scale_y_log10() +
        labs(title = paste("W3 as function of shape | n = ", myNObs)) |

        # Bland-Altman
        ggplot(data = plDatW3i,
               mapping = aes(x = W3, y = W3-W3pred)) +
        geom_hline(yintercept = 0, col = "grey", linetype = "dashed") +
        geom_point() +
        scale_x_log10()

    })

  }# fi

  # gather coefficients of individual generalized logistic functions per nObs
  # residual std. dev increases with nObs
  purrr::map(fm_W3_indiv, .f = coef) |>
    purrr::list_transpose(simplify = TRUE) |>
    as_tibble() |>
    dplyr::mutate(nObs = as.numeric(names(fm_W3_indiv)),
                  resStdDev = purrr::map_dbl(fm_W3_indiv, .f = sigma)) |>
    dplyr::relocate(nObs)
    ## |> dplyr::arrange(desc(resStdDev))
})

if (rlang::is_interactive()) {

  # look at model error for Richards logistic models
  ggplot(data = .MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = resStdDev)) +
    geom_point() +
    labs(title = "Residual std. deviation") +
    scale_x_log10() +
    scale_y_log10()

  local({
    # check how parameters depend on nObs
    # regular pattern desired so that we can approximate this
    nObs_v <- .MLEw_approx[["coef"]][["W3_richards"]]$nObs
    opar <- par(mfrow = c(2,3))
    # walk across parameters
    purrr::iwalk(.MLEw_approx[["coef"]][["W3_richards"]][-1L],
                 .f = ~plot(x = nObs_v, y = .x, main = paste("parameter", .y)))
    par(opar)
  })
} #fi


# # spline interpolation through individual curve coefficients
# dat_ips_A  <- with(.MLEw_approx[["coef"]][["W3_richards"]], interpSpline(A ~ nObs))
# dat_ips_K  <- with(.MLEw_approx[["coef"]][["W3_richards"]], interpSpline(K ~ nObs))
# dat_ips_B  <- with(.MLEw_approx[["coef"]][["W3_richards"]], interpSpline(B ~ nObs))
# dat_ips_Q  <- with(.MLEw_approx[["coef"]][["W3_richards"]], interpSpline(Q ~ nObs))
# dat_ips_nu <- with(.MLEw_approx[["coef"]][["W3_richards"]], interpSpline(nu ~ nObs))

# interpolation spline
.MLEw_approx[["coef"]][["W3_richards_ips"]] <- with(data = .MLEw_approx[["coef"]][["W3_richards"]],
                                                    expr = list(A =  splines::interpSpline(A ~ nObs),
                                                                K =  splines::interpSpline(K ~ nObs),
                                                                Q =  splines::interpSpline(Q ~ nObs),
                                                                B =  splines::interpSpline(B ~ nObs),
                                                                nu = splines::interpSpline(nu ~ nObs)))

if (rlang::is_interactive()) {
  nObs_interp <- sort(unique(c(.MLEw_approx[["coef"]][["W3_richards"]]$nObs[-seq_len(7)],
                               5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000)))

  ggplot(data = .MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = A)) +
    geom_point() +
    geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["A"]], x = nObs_interp)),
              mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
    #  geom_line(data = dat_sp2, mapping = aes(x = x, y = y), col = "blue") +
    #  geom_line(data = dat_sp3, mapping = aes(x = x, y = y), col = "darkgreen") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Parameter A")

  ggplot(data = .MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = K)) +
    geom_point() +
    geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["K"]], x = nObs_interp)),
              mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Parameter K")

  ggplot(.MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = B)) +
    geom_point() +
    geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["B"]], x = nObs_interp)),
              mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
    #  geom_line(data = dat_sp2, mapping = aes(x = x, y = y), col = "blue") +
    #  geom_line(data = dat_sp3, mapping = aes(x = x, y = y), col = "darkgreen") +
    scale_x_log10() +
    #scale_y_log10() +
    labs(title = "Parameter B")

  ggplot(.MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = Q)) +
    geom_point() +
    geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["Q"]], x = nObs_interp)) %>%
                # avoid negative values
                mutate(y = pmax.int(sqrt(.Machine$double.eps), y)),
              mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
    scale_x_log10() +
    #scale_y_log10() +
    labs(title = "Parameter Q")

  ggplot(.MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = nu)) +
    geom_point() +
    geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["nu"]], x = nObs_interp)),
              mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
    scale_x_log10() +
    #scale_y_log10() +
    labs(title = "Parameter nu")
} #fi


#' Factory method to get weight function for W3 for a given sample size
#'
#' Generally, the weight W3 depends on the sample size and the shape parameter.
#' Here, we return a function that returns W3 for given shape parameter.
#' @param nObs sample size for which to return the W3-function
#' @return W3-function for the given sample size. The function returns the W3 weight for the given shape
w3FF <- function(nObs) {

  if (missing(nObs) || length(nObs) != 1L || !is.numeric(nObs) || !is.finite(nObs)) {
    stop("Please provide the number of observations within group!", call. = FALSE)
  }


  # catch all for n = 1
  if (nObs < 2L) return(function(k) 1)

  # get coefficients for a Richards generalized logistic function
  # if we have fit the parameter nObs directly we use the Richards fit
  #+otherwise, we rely on the spline approximation of the fit.
  approx_W3_ind <- which(.MLEw_approx[["coef"]][["W3_richards"]]$nObs == nObs)
  approx_W3_names <- c("A", "K", "Q", "B", "nu")

  # check for match in W3_richards
  approx_W3_coefs <- if (length(approx_W3_ind) == 1L) {
    .MLEw_approx[["coef"]][["W3_richards"]][approx_W3_ind, approx_W3_names]
  } else {
    # nObs was not in MC-sim for W3
    stopifnot(setequal(approx_W3_names, names(.MLEw_approx[["coef"]][["W3_richards_ips"]])))
    purrr::map(.x = .MLEw_approx[["coef"]][["W3_richards_ips"]][approx_W3_names],
               # use predict from splines
               .f = ~ max(0, predict(.x, x = nObs)$y)) %>%
      rlang::set_names(nm = approx_W3_names)
  } #esle

  # return:
  # fn of shape k
  #W3 ~ A + (K - A) / (1 + Q * exp(-B * lshape))**(1/nu)
  function(k) {
    evalq(expr = A + (K - A) / (1 + Q * k**-B)**(1/nu),
          envir = as.list(approx_W3_coefs),
          enclos = rlang::current_env())
  }#fn

}#fn w3FF



# save results ------------------------------------------------------------

# append weight functions
.MLEw_approx <- append(.MLEw_approx,
                       values = list(fun = list(w1F = w1F,
                                                w2F = w2F,
                                                w3FF = w3FF)))

message("Save MLEw-weights info (MCS and approx) as RData-file ", rdataFile)
save(list = c(".MLEw_mcs", ".MLEw_approx"), file = rdataFile)


# exit --------------------------------------------------------------------

# tear-down
future::plan(future::sequential())

# output the latest warnings:
message("\n\n+++\nThese are warnings from the script:\n+++\n")
warnings()


message("~~ Fine ~~")
message("Finished script at ", toString(Sys.time()))
