# mkuhn, 2023-04-11
# adds MLE-weights table as internal data to package as list 'MLEw_approx'
#
# the MLE-weights are established through a Monte-Carlo simulation (MCS),
# see inst/scripts/simul_MLEweights.R which produces 'MLEw_mcs.rds'
# Here, we build approximating function based on this MCS.
# Code to explore which are good/best approximations W1, W2 and W3 are in scratch/MLEw_weights2.R.
# The package incubate makes use of these functions in MLEw_approx[["fun"]]
####


# init --------------------------------------------------------------------

library("usethis")

library("dplyr")
library("patchwork")

library("gslnls")
library("splines")
#library("matrixStats", warn.conflicts = FALSE)


# start from directory "data-raw/"
FNAME <- "MLEw_mcs.rds"
stopifnot(file.exists(FNAME))
#(load(FNAME))
MLEw_mcs <- readRDS(FNAME)


stopifnot(is.list(MLEw_mcs))
stopifnot(identical(names(MLEw_mcs), c("W12","W3", "settings")),
          is.data.frame(MLEw_mcs$W12), is.data.frame(MLEw_mcs$W3))


W3 <- MLEw_mcs$W3 |>
  dplyr::filter(nObs > 1) |>
  dplyr::mutate(
    lshape = log(shape),
    nlshape = -lshape,
    nlp1shape = -log1p(shape),
    lW3 = log(W3),
    llW3 = log(lW3),
    lp1lW3 = log1p(lW3))



# for which Ns do we use direct numbers
N_DIRECT <- 49L
# check that we have all data stored for all consecutive nObs starting from 1!
# we assume this within w1F and w2F below!
stopifnot(NROW(MLEw_mcs$W12) > N_DIRECT, NROW(MLEw_mcs$W3) > N_DIRECT)
stopifnot(identical(MLEw_mcs$W12$nObs[seq_len(N_DIRECT)], seq_len(N_DIRECT)))
stopifnot(identical(MLEw_mcs$W3 |>
                      dplyr::distinct(nObs) |>
                      dplyr::slice_head(n=N_DIRECT) |>
                      dplyr::pull(nObs),
                    seq_len(N_DIRECT)))


# read in weights from publication of Cousineau (2009):
# Cousineau did a relatively small MC-simulation study
W_cousineau2009 <- local({
  # from Cousineau Table 2
  W1_mcss_str <- "1 1.000 0.561 0.693 2 1.000 0.763 0.839 3 1.000 0.839 0.891 4 1.000 0.878 0.918 5 1.000 0.902 0.934 6 1.000 0.918 0.945 7 1.000 0.929 0.953 8 1.000 0.938 0.959 9 1.000 0.945 0.963 10 1.000 0.950 0.967 11 1.000 0.955 0.970 12 1.000 0.959 0.972 13 1.000 0.962 0.974 14 1.000 0.965 0.976 15 1.000 0.967 0.978 16 1.000 0.969 0.979"
  # from Cousineau Table 3
  W2_mcss_str <- "1 0.000 0.000 0.000 2 0.500 0.163 0.275 3 0.667 0.409 0.517 4 0.750 0.553 0.638 5 0.800 0.642 0.711 6 0.833 0.702 0.759 7 0.857 0.742 0.791 8 0.875 0.775 0.817 9 0.889 0.800 0.838 10 0.900 0.820 0.853 11 0.909 0.835 0.867 12 0.917 0.849 0.877 13 0.923 0.860 0.886 14 0.929 0.871 0.895 15 0.933 0.879 0.902 16 0.938 0.887 0.908"
  # from Cousineau Table 4
  W3_mcss_str <- "1 12.429 20.157 11.371 12.483 24.796 2 19.452 11.567 4.183 3.147 2.771 3 33.320 14.372 3.701 2.596 2.225 4 73.132 36.570 3.480 2.411 2.043 5 66.530 14.812 3.431 2.309 1.948 6 87.540 11.751 3.297 2.235 1.888 7 124.230 21.331 3.270 2.198 1.852 8 99.608 13.230 3.192 2.170 1.831 9 97.148 14.622 3.178 2.154 1.808 10 447.600 19.335 3.244 2.143 1.794 11 105.660 13.879 3.195 2.120 1.779 12 164.510 13.765 3.154 2.113 1.769 13 136.390 12.762 3.109 2.109 1.759 14 342.220 14.270 3.111 2.099 1.755 15 168.430 14.737 3.110 2.091 1.746 16 198.130 13.641 3.101 2.093 1.742 1 1.004 1.001 1.001 1.006 1.002 2 2.395 1.851 1.524 1.360 1.268 3 3.682 2.375 1.775 1.520 1.383 4 4.854 2.753 1.934 1.603 1.438 5 6.005 3.046 2.042 1.665 1.479 6 7.097 3.295 2.119 1.704 1.506 7 8.103 3.497 2.184 1.740 1.528 8 9.145 3.665 2.229 1.766 1.543 9 10.133 3.842 2.278 1.778 1.552 10 11.104 3.966 2.312 1.800 1.564 11 12.190 4.083 2.354 1.817 1.572 12 13.019 4.192 2.378 1.829 1.583 13 13.898 4.280 2.403 1.839 1.582 14 14.857 4.367 2.423 1.842 1.591 15 15.819 4.493 2.440 1.854 1.595 16 16.604 4.561 2.464 1.862 1.598 1 1.001 0.999 0.999 0.995 0.998 2 2.096 1.668 1.456 1.339 1.262 3 3.081 2.082 1.680 1.479 1.367 4 3.950 2.381 1.822 1.567 1.428 5 4.806 2.631 1.920 1.625 1.464 6 5.631 2.808 2.004 1.669 1.492 7 6.433 2.982 2.056 1.698 1.509 8 7.150 3.114 2.105 1.722 1.525 9 7.931 3.252 2.151 1.739 1.537 10 8.643 3.365 2.180 1.758 1.552 11 9.319 3.462 2.207 1.774 1.555 12 10.051 3.560 2.239 1.782 1.565 13 10.746 3.642 2.262 1.793 1.570 14 11.379 3.713 2.285 1.804 1.578 15 12.069 3.780 2.301 1.813 1.581 16 12.743 3.854 2.324 1.820 1.586"

  W1_mcss <- strsplit(W1_mcss_str, split = " ", fixed = TRUE)[[1]] |>
    as.numeric() |>
    matrix(ncol = 4, byrow = TRUE, dimnames = list(1:16, c("n", "E", "G", "J"))) |>
    as.data.frame() |>
    tidyr::pivot_longer(cols = c(E, G, J), names_to = "location", values_to = "value") |>
    dplyr::mutate(shape = NA_real_, .before = value) |>
    dplyr::relocate(location)


  W2_mcss <- strsplit(W2_mcss_str, split = " ", fixed = TRUE)[[1]] |>
    as.numeric() |>
    matrix(ncol = 4, byrow = TRUE, dimnames = list(1:16, c("n", "E", "G", "J"))) |>
    as.data.frame() |>
    tidyr::pivot_longer(cols = c(E, G, J), names_to = "location", values_to = "value") |>
    dplyr::mutate(shape = NA_real_, .before = value) |>
    dplyr::relocate(location)

  W3_mcss <- strsplit(W3_mcss_str, split = " ", fixed = TRUE)[[1]] |>
    as.numeric() |>
    matrix(ncol = 6, byrow = TRUE, dimnames = list(paste0(rep(c("E", "G", "J"), each = 16), rep.int(1:16, times = 3)),
                                                   c("n", paste("shape", c(0.5, 1, 1.5, 2, 2.5), sep = "_")))) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "rown") |>
    dplyr::mutate(location = stringr::str_sub(rown, end = 1), .before = n,
                  rown = NULL) |>
    tidyr::pivot_longer(cols = starts_with("shape_"),
                        names_to = "shape", names_prefix = "^shape_", names_transform = as.numeric,
                        values_to = "value")


  dplyr::bind_rows(list(W1 = W1_mcss, W2 = W2_mcss, W3 = W3_mcss), .id = "type")
})


# build internal data -----------------------------------------------------

message("Start with building approximations for the weights!")

#MLEw_approx is exported as internal data in the end!

# MLEw_mcs is not exported in its entirety to internal data.
# Instead, we save only relevant bits for MLEw-approximation
MLEw_approx <- list(
  # store simulation results for W1 and W2
  MCsim = MLEw_mcs$W12 |>
    dplyr::slice_head(n=N_DIRECT) |>
    dplyr::mutate(W1gamma = stats::qgamma(p = 0.5, shape = nObs, rate = nObs),
                  .after = W1) |>
    as.list(),
  MCsim_cousineau2009 = W_cousineau2009
)




# approximation W2 --------------------------------------------------------

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
MLEw_approx$coef <- list(
  #was c(R0 = -0.44193638, -exp(-0.00624712316)
  W2 = c(R0 = coef(fm_W2)[[1L]],
         negRate = -exp(coef(fm_W2)[[2L]]))
)








if (rlang::is_interactive()) {
  plDatW2 <- MLEw_mcs$W12 |>
    dplyr::filter(nObs > 1) |>
    dplyr::mutate(W2mod = predict(fm_W2))

  ggplot(data = plDatW2,
         mapping = aes(x = nObs, y = W2)) +
    geom_point(alpha = .57, col = "lightgrey", size = 2) + #geom_line() +
    geom_point(mapping = aes(y = W2mod), size = .5, col = "darkred", alpha = .2) +
    geom_line(mapping = aes(y = W2mod), col = "darkred", alpha = .1, linetype = "dotted") +
    geom_vline(xintercept = N_DIRECT, col = "darkgrey") +
    scale_x_log10() +
    labs(title = "W2: median of MC-sim and fitted W2-function",
         subtitle = "Modell: darkred") |


    ggplot(data = plDatW2 |> dplyr::filter(nObs > N_DIRECT),
           mapping = aes(x = W2, y = W2mod-W2)) +
    geom_point() +
    geom_hline(yintercept = 0, col = "darkgrey") +
    #scale_y_sqrt() +
    labs(title = "Bias of modelled W2 vs MC-simulation (median)",
         subtitle = paste("beyond n =", N_DIRECT))
} #fi




# approximation W3 --------------------------------------------------------


# Generalized logistic function
# (due to Richards, 1957)
# in the unified formulation (Tjorve, 2010)
genLogisticF <- function(theta, xVal) {
  theta <- as.numeric(theta)
  stopifnot(length(theta) == 5)

  L <- theta[1]  #lower asymp
  A <- theta[2]  #upper asymp
  d <- theta[3]  #inflection value is (L + (A-L) * d^(1/(1-d)))
  K <- theta[4]  #actual relative growth rate (slope at infl point is (A-L)*K)
  Xi <- theta[5] #inflection point (x-value)

  L + (A - L) * (1 + (d-1) * exp(-K * (xVal - Xi)/d^(d/(1-d))))^(1/(1-d))
}

# Derivative of generalized logistic Richards function
# as function of x-value (for given parameters theta)
genLogisticD <- function(xVal, theta) {
  theta <- as.numeric(theta)
  stopifnot(length(theta) == 5)

  L <- theta[1]  #lower asymp
  A <- theta[2]  #upper asymp
  d <- theta[3]  #inflection value is (L + (A-L) * d^(1/(1-d)))
  K <- theta[4]  #actual relative growth rate (slope at infl point is (A-L)*K)
  Xi <- theta[5] #inflection point (x-value)

  dExpV <- d^(d/(1-d))
  expV <- exp(-K * (xVal - Xi)/dExpV)

  (A-L) * (1 + (d-1) * expV)^(d/(1-d)) * expV * K / dExpV
}#fn

# Gradient of generalized logistic Richards function
# as function of parameters (for given x value)
# using unified parametrization (Tjorve)
# Is used for fitting Richards function to data
genLogisticJ <- function(theta, xVal) {
  theta <- as.numeric(theta)
  stopifnot(is.numeric(theta), length(theta) == 5)

  L <- theta[1]  #lower asymp
  A <- theta[2]  #upper asymp
  d <- theta[3]  #inflection value is (Ad^(1/(1-d)))
  K <- theta[4]  #actual relative growth rate (slope at infl point is AK)
  Xi <- theta[5] #inflection point (x-value)

  stopifnot(is.numeric(xVal))
  nObs <- length(xVal)

  dExpV <- d^(d/(1-d))
  expV <- exp(-K * (xVal - Xi)/dExpV)
  bracV <- 1 + (d-1) * expV

  # return gradient: nObs x 5 matrix
  cbind(
    L = 1 - bracV^(1/(1-d)), #L
    A = bracV^(1/(1-d)),   #A
    d = (A-L) *  bracV^(1/(1-d)) * 1/(1-d) * (1/(1-d) * log(bracV) + expV * (K * (xVal - Xi)/dExpV * (- 1/(1-d) * log(d) - 1) + 1)/bracV), #d
    K = (A-L) *  bracV^(d/(1-d)) * expV *(xVal - Xi) / dExpV, #K
    Xi = -(A-L) * bracV^(d/(1-d)) * expV * K / dExpV #Xi
  )
}#fn grad

MLEw_approx$fun <- list(genLogisticF = genLogisticF,
                        genLogisticD = genLogisticD,
                        genLogisticJ = genLogisticJ)



# check that we do not need too many iterations (sign for difficult/bad fit?!)
ITER_MAX <- 97

MLEw_approx[["coef"]][["W3_richards"]] <- local({

  #we do not fix A=0 on lp1lW3 scale (even though it might be true)
  #+because we want best fit and not too many restrictions for the Richards fit

  # for each sample size nObs we fit
  # Richards' generalized logistic function (as function of shape)
  # This way we can estimate W3 for each shape value for those nObs
  #+that we have simulated in .MLEw_mcs$W3
  fm_W3_indiv <- purrr::map(.x = rlang::set_names(unique(W3$nObs)),
                            .f = function(n_) {
                              W3i <- W3 |>
                                dplyr::filter(nObs == {n_})
                              stopifnot(NROW(W3i) > 5)

                              gsl_nls(fn = MLEw_approx$fun$genLogisticF,
                                      jac = MLEw_approx$fun$genLogisticJ,
                                      y = W3i$lp1lW3, xVal = W3i$nlshape,
                                      start = list(L = 0, A = log1p(log(2*{n_})),
                                                   d = max(2,log({n_})), K = .5, Xi = 0),
                                      lower = c(L = -.1, A = .2, d = 2, K = .2, Xi = -10),
                                      ##weights = sqrt(W3i$shape),
                                      control = gsl_nls_control(maxiter = ITER_MAX))
                            })

  # check convergence for each model
  stopifnot(all(purrr::map_lgl(fm_W3_indiv, .f = list("convInfo", "isConv"))))
  stopifnot(all(purrr::map_dbl(fm_W3_indiv, .f = list("convInfo", "nEval", "f")) < ITER_MAX))


  # gather coefficients of individual generalized logistic functions per nObs
  # residual std. dev increases with nObs
  purrr::map(fm_W3_indiv, .f = coef) |>
    purrr::list_transpose(simplify = TRUE) |>
    as_tibble() |>
    dplyr::mutate(nObs = as.numeric(names(fm_W3_indiv)),
                  resStdDev = purrr::map_dbl(fm_W3_indiv, .f = stats::sigma)) |>
    dplyr::relocate(nObs)
})



# checking ----------------------------------------------------------------

# test jacobian for some sample size nObs
myN <- sample(unique(W3$nObs), size = 1)
stopifnot(length(myN) == 1)
myShapes <- stats::rlnorm(n=7, meanlog = .51, sdlog = 2) |>
  sort()


# some parameters values
startL <- list(L = .01,
               A = log1p(log(2*myN)),
               d = 3.5+stats::rnorm(n=1, mean = .2, sd = .1),
               K = .01 + abs(stats::rnorm(n=1, mean = .5, sd = .1)),
               Xi = stats::rnorm(n=1, mean = -1, sd = .1))

# test gradient function (partial derivatives of parameters)
all.equal(purrr::map(.x = myShapes,
                     .f = ~numDeriv::grad(func = MLEw_approx$fun$genLogisticF,
                                          x = as.numeric(startL), xVal = .x)) |>
            # convert to single matrix, columns = nbr of parameters
            unlist() |> matrix(ncol = 5, byrow = TRUE,
                               dimnames = list(NULL, c("L", "A", "d", "K", "Xi"))),

          #current=
          MLEw_approx$fun$genLogisticJ(theta = as.numeric(startL),
                                       xVal = myShapes),
          tolerance = 1e-7) |>
  isTRUE() |>
  stopifnot()

all.equal(
  numDeriv::grad(func = MLEw_approx$fun$genLogisticF,
                 x = myShapes, theta = startL),
  MLEw_approx$fun$genLogisticD(xVal = myShapes, theta = startL),
  tolerance = 1e-7
) |> isTRUE() |>
  stopifnot()

if (rlang::is_interactive()) {

  # check how parameters depend on nObs
  # regular pattern desired so that we can approximate this
  local({
    nObs_v <- MLEw_approx[["coef"]][["W3_richards"]]$nObs
    opar <- par(mfrow = c(2,3))
    # walk across parameters
    purrr::iwalk(MLEw_approx[["coef"]][["W3_richards"]][-1L],
                 .f = ~plot(x = nObs_v, y = .x, log = "x",
                            type = "p", pch = 16, cex = 0.5,
                            ylab = .y, main = paste("parameter", .y)))
    par(opar)
  })

  # focus on std. dev.
  # look at model error for Richards logistic models
  ggplot(data = MLEw_approx[["coef"]][["W3_richards"]],
         mapping = aes(x = nObs, y = resStdDev)) +
    geom_point() +
    labs(title = "Residual std. deviation") +
    scale_x_log10() +
    scale_y_log10()


  # W3 approaches 1 for high shape values
  #+hence, log1p(W3) => 0
  local({
    # reuse previous myN
    stopifnot(exists("myN"), is.finite(myN))

    myTheta <- MLEw_approx$coef$W3_richards |>
      dplyr::filter(nObs == myN) |>
      dplyr::select(!c(nObs, resStdDev)) |>
      unlist()

    W3i <- W3 |>
      dplyr::filter(nObs == {{myN}}) |>
      #dplyr::mutate(lp1lW3_pred = predict(fm_W3_indiv[[as.character(myN)]]))
      dplyr::mutate(lp1lW3_pred = MLEw_approx$fun$genLogisticF(theta = myTheta, xVal = nlshape))

    ggplot(data = W3i,
           mapping = aes(x = shape, y = lp1lW3)) +
      geom_point() + #geom_line() +
      geom_point(mapping = aes(y = lp1lW3_pred), size = .5, col = "darkred") +
      geom_line(mapping = aes(y = lp1lW3_pred), col = "darkred") +
      scale_x_log10() + #scale_y_log10() +
      labs(title = paste("W3 as function of shape || n = ", myN)) |

      # Bland-Altman
      ggplot(data = W3i,
             mapping = aes(x = W3, y = lp1lW3-lp1lW3_pred)) +
      geom_hline(yintercept = 0, col = "grey", linetype = "dashed") +
      geom_point() +
      scale_x_log10()

  })
} #fi interactive





# save as internal data ---------------------------------------------------

message("Save MLEw-weights approximation functions as internal data")

usethis::use_data(MLEw_approx, internal = TRUE, overwrite = TRUE)


message("~~ Fine ~~")


#q(save = "no")
