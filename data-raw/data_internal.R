# mkuhn, 2023-04-11
# adds MLE-weights table as internal data to package
#
# the MLE-weights are established through a Monte-Carlo simulation (mcs), see inst/scripts/simul_MLEweights.R
# Based on this, we build here approximating function.
# Code to explore which are good/best approximations W1, W2 and W3 are in scratch/MLEw_weights2.R.
# The package incubate makes use of these functions in .MLEw_approx[["fun"]]
####


# init --------------------------------------------------------------------

library("usethis")

library("dplyr")
library("patchwork")

library("gslnls")
library("splines")
library("matrixStats", warn.conflicts = FALSE)


FNAME <- "MLEw_mcs.rds"
stopifnot(file.exists(FNAME))
#(load(FNAME))
MLEw_mcs <- readRDS(FNAME)


stopifnot(is.list(MLEw_mcs))
stopifnot(identical(names(MLEw_mcs), c("W12","W3", "settings")),
          is.data.frame(MLEw_mcs$W12), is.data.frame(MLEw_mcs$W3))


# for which Ns do we use direct numbers
N_DIRECT <- 50L
# check that we have all data stored for all consecutive nObs starting from 1!
# we assume this within w1F and w2F below!
stopifnot(NROW(MLEw_mcs$W12) > N_DIRECT, NROW(MLEw_mcs$W3) > N_DIRECT)
stopifnot(identical(MLEw_mcs$W12$nObs[seq_len(N_DIRECT)], seq_len(N_DIRECT)))
stopifnot(identical(MLEw_mcs$W3 |>
                      dplyr::distinct(nObs) |>
                      dplyr::slice_head(n=N_DIRECT) |>
                      dplyr::pull(nObs), seq_len(N_DIRECT)))

# MLEw_mcs is not exported to internal data.
# Instead, we save relevant infos for MLEw-approximation in .MLEw_approx
#+that we export as internal data
.MLEw_approx <- list(
  # store simulation results for W1 and W2
  MCsim = MLEw_mcs$W12 |>
    dplyr::slice_head(n=N_DIRECT) |>
    as.list()
)


# asymptotic regression model: cf. SSasymp model on log(nObs)
# starting at nObs = 2 (as nObs = 1 is off).
fm_W2 <- gsl_nls(W2 ~ 1 + (R0 - 1) * nObs**-exp(lr),
                 start = list(R0 = -.25, lr = -.01),
                 data = MLEw_mcs$W12, subset = nObs > 1)

# check model fit
if (!fm_W2$convInfo$isConv || fm_W2$convInfo$stopCode != 0 || fm_W2$convInfo$nEval[["f"]] > 27 || deviance(fm_W2) > 1e-3) {
  stop("Model fit for W2 is bad!")
}


# add coef for W2-approximation
.MLEw_approx$coef <- list(
  #was c(R0 = -0.44193638, -exp(-0.00624712316)
  W2 = c(R0 = coef(fm_W2)[[1L]],
         negRate = -exp(coef(fm_W2)[[2L]]))
)




# approximation W1 ----------------------------------------------------------

message("Start with building approximations for the weights!")


# MLEw-weights W1
# W1 for given sample sizes (of one group).
# For small `nObs` we use direct results from Monte-Carlo simulation.
# For higher `nObs` we use an approximation (based on Wilson-Hilferty transformation)
# @param nObs numeric. number of observations (vectorized)
# @return numeric. W1-value corrsponding to nObs. Same length as nObs
w1F <- function(nObs) {

  if (missing(nObs) || !is.numeric(nObs) || any(!is.finite(nObs))) {
    stop("Please provide the number of observations within a group. Must be numeric and finite!", call. = FALSE)
  }

  nObs <- pmax.int(1L, nObs)

  # nbr of W1 simulation results we use directly
  W1_MCsim <- .MLEw_approx[["MCsim"]][["W1"]]

  nObsIdx_direct <- which(nObs <= length(W1_MCsim))
  nObsIdx_approx <- which(nObs > length(W1_MCsim))

  retV <- numeric(length(nObs))
  retV[nObsIdx_direct] <- W1_MCsim[nObs[nObsIdx_direct]]
  # approximation for median of gamma(n, 1/n)
  #+using Wilson-Hilferty transformation (see <https://en.wikipedia.org/wiki/Gamma_distribution>)
  retV[nObsIdx_approx] <- (1 - 1 / (9 * nObs[nObsIdx_approx]))^3

  retV

}#fn w1F


# approximation W2 --------------------------------------------------------

# MLEw weight W2
# For given sample size of one group
# @param nObs numeric. Sample size (vectorized)
# @return W2 (same size as nObs)
w2F <- function(nObs) {

  if (missing(nObs) || !is.numeric(nObs) || any(!is.finite(nObs))) {
    stop("Please provide the number of observations within group!", call. = FALSE)
  }

  nObs <- pmax.int(1L, nObs)

  # nbr of W2 simulation results we use directly
  W2_MCsim <- .MLEw_approx[["MCsim"]][["W2"]]
  # median approximation via asymptotic regression model SSasymp on log(n):
  # We hence model: W2 = 1 + (R0 - 1) * nObs**(-r)
  W2_coef <- .MLEw_approx[["coef"]][["W2"]]

  nObsIdx_direct <- which(nObs <= length(W2_MCsim))
  nObsIdx_approx <- which(nObs > length(W2_MCsim))


  retV <- numeric(length(nObs))
  retV[nObsIdx_direct] <- W2_MCsim[nObs[nObsIdx_direct]]
  retV[nObsIdx_approx] <- 1 + (W2_coef[["R0"]] - 1) * nObs[nObsIdx_approx]**W2_coef[["negRate"]]

  retV
}#fn w2F


if (rlang::is_interactive()) {
  plDatW2 <- MLEw_mcs$W12 |>
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
    #scale_y_sqrt() +
    labs(title = "Agreement of W2 (median)", subtitle = "approx vs simulation")
} #fi




# approximation W3 --------------------------------------------------------


# Generalized logistic function
# (due to Richards, 1957)
genLogisticF <- function(theta, xVal) {
  theta <- as.numeric(theta)
  stopifnot(length(theta) == 5)

  A <- theta[1]
  K <- theta[2]
  Q <- theta[3]
  B <- theta[4]
  nu <- theta[5]

  A + (K - A) / (1 + Q * exp(-B * xVal))**(1/nu)
}

genLogisticJ <- function(theta, xVal) {
  #cat("theta:", paste(theta, sep = ","), "\n")
  #cat(class(theta), "\n")
  theta <- as.numeric(theta)

  stopifnot(is.numeric(xVal))
  nObs <- length(xVal)

  stopifnot(is.numeric(theta), length(theta) == 5)
  #theta <- rlang::set_names(theta, nm = c("A", "K", "Q", "B", "nu"))
  A <- theta[1]
  K <- theta[2]
  Q <- theta[3]
  B <- theta[4]
  nu <- theta[5]

  # return
  cbind(
    1 - (1 + Q * exp(-B*xVal))^(-1/nu), #A
    (1 + Q * exp(-B*xVal))^(-1/nu), #K
    -(K-A) * exp(-B * xVal) / (nu * (1+Q * exp(-B * xVal))^(1 + 1/nu)), #Q
    (K-A) * xVal * Q * exp(-B * xVal) / (nu * (1+Q * exp(-B * xVal))^(1 + 1/nu)), #B
    (K-A) * log1p(Q * exp(-B * xVal)) / (nu^2 * (1 + Q * exp(-B * xVal))^(1/nu)) #nu
  )
}

.MLEw_approx[["fun"]] <- list(genLogisticF = genLogisticF,
                              genLogisticJ = genLogisticJ)

W3 <- MLEw_mcs$W3 |>
  dplyr::filter(nObs > 1) |>
  dplyr::mutate(
    lshape = log(shape),
    nlshape = -lshape,
    nlp1shape = -log1p(shape),
    lW3 = log(W3),
    llW3 = log(lW3),
    lp1lW3 = log1p(lW3))


# test jacobian for some sample size nObs
myN <- sample(unique(W3$nObs), size = 1)
W3i <- W3 |> dplyr::filter(nObs == myN)
stopifnot(length(myN) == 1)
stopifnot(exists("W3i"), NROW(W3i) > 1)

# some parameters values
startL <- list(A = 0.01, K = log1p(log(2*myN)), Q = 2, B = 1.5, nu = .75)

# test gradient function
idx <- sort(sample(x = NROW(W3i), size = 5, replace = FALSE))
waldo::compare(
  x = purrr::map(.x = W3i$lshape[idx],
                 .f = ~numDeriv::grad(func = .MLEw_approx$fun$genLogisticF,
                                      x = as.numeric(startL), xVal = .x)) |>
    # convert to single matrix, columns = nbr of parameters
    unlist() |> matrix(ncol = 5, byrow = TRUE),

  y = .MLEw_approx$fun$genLogisticJ(theta = as.numeric(startL),
                                    xVal = W3i$lshape[idx]),
  tolerance = 1e-5)



# check that we do not need too many iterations (sign for difficult/bad fit?!)
ITER_MAX <- 97

.MLEw_approx[["coef"]][["W3_richards"]] <- local({

  #we do not fix A=0 on lp1lW3 scale (even though it might be true)
  #+because it restricts the Richards fit too much?!

  # for each sample size nObs we fit
  # Richards' generalized logistic function (as function of shape)
  # This way we can estimate W3 for each shape value for those nObs
  #+that we have simulated in .MLEw_mcs$W3
  fm_W3_indiv <- purrr::map(.x = rlang::set_names(unique(W3$nObs)),
                            .f = function(n_) {
                              W3i <- W3 |>
                                dplyr::filter(nObs == {n_})
                              stopifnot(NROW(W3i) > 5)
                              gsl_nls(fn = .MLEw_approx$fun$genLogisticF,
                                      jac = .MLEw_approx$fun$genLogisticJ,
                                      y = W3i$lp1lW3, xVal = W3i$nlshape,
                                      start = list(A = 0.01, K = log1p(log(2*{n_})), Q = 2, B = 1.5, nu = .75),
                                      lower = c(A = -.25, K = .25, Q = 1e-3, B = -Inf, nu = 1e-3),
                                      ##weights = sqrt(W3i$shape),
                                      control = gsl_nls_control(maxiter = ITER_MAX))
                            })

  # check convergence for each model
  stopifnot(all(purrr::map_lgl(fm_W3_indiv, .f = list("convInfo", "isConv"))))
  stopifnot(all(purrr::map_dbl(fm_W3_indiv, .f = list("convInfo", "nEval", "f")) < ITER_MAX))


  if (rlang::is_interactive()) {

    # W3 approaches 1 for high shape values
    #+hence, log1p(W3) => 0
    local({
      # reuse previous myN
      stopifnot(exists("myN"), is.finite(myN))

      W3i <- W3i |>
        dplyr::mutate(lp1lW3_pred = predict(fm_W3_indiv[[as.character(myN)]]))

      ggplot(data = W3i,
             mapping = aes(x = shape, y = lp1lW3)) +
        geom_point() + #geom_line() +
        geom_point(mapping = aes(y = lp1lW3_pred), size = .5, col = "darkred") +
        geom_line(mapping = aes(y = lp1lW3_pred), col = "darkred") +
        scale_x_log10() + #scale_y_log10() +
        labs(title = paste("W3 as function of shape | n = ", myN)) |

        # Bland-Altman
        ggplot(data = W3i,
               mapping = aes(x = W3, y = lp1lW3-lp1lW3_pred)) +
        geom_hline(yintercept = 0, col = "grey", linetype = "dashed") +
        geom_point() +
        scale_x_log10()

    })
  }#fi interactive


  # gather coefficients of individual generalized logistic functions per nObs
  # residual std. dev increases with nObs
  purrr::map(fm_W3_indiv, .f = coef) |>
    purrr::list_transpose(simplify = TRUE) |>
    as_tibble() |>
    dplyr::mutate(nObs = as.numeric(names(fm_W3_indiv)),
                  resStdDev = purrr::map_dbl(fm_W3_indiv, .f = stats::sigma)) |>
    dplyr::relocate(nObs)
})

# check
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
                 .f = ~plot(x = nObs_v, y = .x, log = "x",
                            ylab = .y, main = paste("parameter", .y)))
    par(opar)
  })
} #fi interactive


# # interpolation spline
# # I would need predict.XXX from splines-package below
# .MLEw_approx[["coef"]][["W3_richards_ips"]] <- with(data = .MLEw_approx[["coef"]][["W3_richards"]],
#                                                     expr = list(A =  splines::interpSpline(A ~ nObs),
#                                                                 K =  splines::interpSpline(K ~ nObs),
#                                                                 Q =  splines::interpSpline(Q ~ nObs),
#                                                                 B =  splines::interpSpline(B ~ nObs),
#                                                                 nu = splines::interpSpline(nu ~ nObs)))
#
#
# if (rlang::is_interactive()) {
#   nObs_interp <- sort(unique(c(.MLEw_approx[["coef"]][["W3_richards"]]$nObs[-seq_len(17)],
#                                5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000)))
#
#   ggplot(data = .MLEw_approx[["coef"]][["W3_richards"]],
#          mapping = aes(x = nObs, y = A)) +
#     geom_point() +
#     geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["A"]], x = nObs_interp)),
#               mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
#     #  geom_line(data = dat_sp2, mapping = aes(x = x, y = y), col = "blue") +
#     #  geom_line(data = dat_sp3, mapping = aes(x = x, y = y), col = "darkgreen") +
#     scale_x_log10() +
#     #scale_y_log10() +
#     labs(title = "Parameter A")
#
#   ggplot(data = .MLEw_approx[["coef"]][["W3_richards"]],
#          mapping = aes(x = nObs, y = K)) +
#     geom_point() +
#     geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["K"]], x = nObs_interp)),
#               mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
#     scale_x_log10() +
#     scale_y_log10() +
#     labs(title = "Parameter K")
#
#   ggplot(.MLEw_approx[["coef"]][["W3_richards"]],
#          mapping = aes(x = nObs, y = B)) +
#     geom_point() +
#     geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["B"]], x = nObs_interp)),
#               mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
#     #  geom_line(data = dat_sp2, mapping = aes(x = x, y = y), col = "blue") +
#     #  geom_line(data = dat_sp3, mapping = aes(x = x, y = y), col = "darkgreen") +
#     scale_x_log10() +
#     #scale_y_log10() +
#     labs(title = "Parameter B")
#
#   ggplot(.MLEw_approx[["coef"]][["W3_richards"]],
#          mapping = aes(x = nObs, y = Q)) +
#     geom_point() +
#     geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["Q"]], x = nObs_interp)) %>%
#                 # avoid negative values
#                 mutate(y = pmax.int(sqrt(.Machine$double.eps), y)),
#               mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
#     scale_x_log10() +
#     #scale_y_log10() +
#     labs(title = "Parameter Q")
#
#   ggplot(.MLEw_approx[["coef"]][["W3_richards"]],
#          mapping = aes(x = nObs, y = nu)) +
#     geom_point() +
#     geom_line(data = as_tibble(predict(.MLEw_approx[["coef"]][["W3_richards_ips"]][["nu"]], x = nObs_interp)),
#               mapping = aes(x = x, y = y), col = "blue", linetype = "dashed") +
#     scale_x_log10() +
#     #scale_y_log10() +
#     labs(title = "Parameter nu")
# } #fi


#' Factory method to get weight function for W3 for a given sample size
#'
#' Generally, the weight W3 depends on the sample size and the shape parameter.
#' Here, we return a function that returns W3 for given shape parameter.
#' @param nObs sample size for which to return the W3-function
#' @return W3-function for the given sample size. The function returns the W3 weight for the given shape
w3FF <- function(nObs) {

  if (missing(nObs) || length(nObs) != 1L || !is.numeric(nObs) || !is.finite(nObs)) {
    stop("Please provide the single number of observations within group!", call. = FALSE)
  }

  # catch all for n = 1 (or n negative etc)
  if (nObs < 2L) return(function(k) 1)

  # get coefficients for a Richards' generalized logistic function
  # if we have fit the parameter nObs directly we use the Richards fit
  #+otherwise, we rely on the spline approximation of the fit.
  approx_W3_ind <- which(.MLEw_approx[["coef"]][["W3_richards"]]$nObs == nObs)
  approx_W3_names <- c("A", "K", "Q", "B", "nu")

  # check for match in W3_richards
  approx_W3_coefs <- if (length(approx_W3_ind) == 1L) {
    .MLEw_approx[["coef"]][["W3_richards"]][approx_W3_ind, approx_W3_names]
  } else {
    # nObs was not in MC-sim for W3
    list(
      A = stats::spline(x = .MLEw_approx[["coef"]][["W3_richards"]]$nObs,
                        y = .MLEw_approx[["coef"]][["W3_richards"]]$A,
                        method = "natural", xout = {nObs})$y,
      K = stats::spline(x = .MLEw_approx[["coef"]][["W3_richards"]]$nObs,
                        y = .MLEw_approx[["coef"]][["W3_richards"]]$K,
                        method = "natural", xout = {nObs})$y,
      Q = stats::spline(x = .MLEw_approx[["coef"]][["W3_richards"]]$nObs,
                        y = .MLEw_approx[["coef"]][["W3_richards"]]$Q,
                        method = "natural", xout = {nObs})$y,
      B = stats::spline(x = .MLEw_approx[["coef"]][["W3_richards"]]$nObs,
                        y = .MLEw_approx[["coef"]][["W3_richards"]]$B,
                        method = "natural", xout = {nObs})$y,
      nu = stats::spline(x = .MLEw_approx[["coef"]][["W3_richards"]]$nObs,
                         y = .MLEw_approx[["coef"]][["W3_richards"]]$nu,
                         method = "natural", xout = {nObs})$y
    )

  } #esle


  # W3 as fn of shape k
  # @param k shape
  # @return: W3 (same length as k)
  # XXX add gradient to this function as attribute?
  function(k) {

    # undo the transformation:
    #+x (predictor) as neg. log(shape)
    #+y (response) as log1p(lW3)
    exp(expm1(.MLEw_approx$fun$genLogisticF(theta = approx_W3_coefs,
                                            xVal = -log(k))))

    # evalq(expr = A + (K - A) / (1 + Q * k**-B)**(1/nu),
    #       envir = as.list(approx_W3_coefs),
    #       enclos = rlang::current_env())
  }#fn

}#fn w3FF




# append weight functions
.MLEw_approx$fun <- append(.MLEw_approx$fun,
                           values = list(w1F = w1F,
                                         w2F = w2F,
                                         w3FF = w3FF))




# save as internal data ---------------------------------------------------

message("Save MLEw-weights approximation functions as internal data")

usethis::use_data(.MLEw_approx, internal = TRUE, overwrite = TRUE)


cat("\n~~ Fine ~~\n")

#q(save = "no")
