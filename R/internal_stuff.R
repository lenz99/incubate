# internal MLEw weight functions


#' internal MLEw-weights W1 function
#' W1 for given sample sizes (of one group).
#' For small `nObs` we use direct results from Monte-Carlo simulation.
#' For higher `nObs` we use an approximation (based on Wilson-Hilferty transformation)
#' @param nObs numeric. number of observations (vectorized)
#' @returns numeric. W1-value corrsponding to nObs. Same length as nObs
w1Fint <- function(nObs) {

  if (missing(nObs) || !is.numeric(nObs) || any(!is.finite(nObs))) {
    stop("Please provide the number of observations within a group. Must be numeric and finite!", call. = FALSE)
  }

  nObs <- pmax.int(1L, nObs)

  # nbr of W1 simulation results we use directly
  W1_MCsim <- MLEw_approx[["MCsim"]][["W1"]]

  nObsIdx_direct <- which(nObs <= length(W1_MCsim))
  nObsIdx_approx <- which(nObs > length(W1_MCsim))

  retV <- numeric(length(nObs))
  retV[nObsIdx_direct] <- W1_MCsim[nObs[nObsIdx_direct]]
  # approximation for median of gamma(n, 1/n)
  #+using Wilson-Hilferty transformation (see <https://en.wikipedia.org/wiki/Gamma_distribution>)
  retV[nObsIdx_approx] <- (1 - 1 / (9 * nObs[nObsIdx_approx]))^3

  retV

}#fn w1Fint


#' internal MLEw weight W2
#' For given sample size of one group
#' @param nObs numeric. Sample size (vectorized)
#' @returns W2 (same size as nObs)
w2Fint <- function(nObs) {

  if (missing(nObs) || !is.numeric(nObs) || any(!is.finite(nObs))) {
    stop("Please provide the number of observations within group!", call. = FALSE)
  }

  nObs <- pmax.int(1L, nObs)

  # nbr of W2 simulation results we use directly
  W2_MCsim <- MLEw_approx[["MCsim"]][["W2"]]
  # median approximation via asymptotic regression model SSasymp on log(n):
  # We hence model: W2 = 1 + (R0 - 1) * nObs**(-r)
  W2_coef <- MLEw_approx[["coef"]][["W2"]]

  nObsIdx_direct <- which(nObs <= length(W2_MCsim))
  nObsIdx_approx <- which(nObs > length(W2_MCsim))


  retV <- numeric(length(nObs))
  retV[nObsIdx_direct] <- W2_MCsim[nObs[nObsIdx_direct]]
  retV[nObsIdx_approx] <- 1 + (W2_coef[["R0"]] - 1) * nObs[nObsIdx_approx]**W2_coef[["negRate"]]

  retV
}#fn w2Fint


#' Internal factory method to get weight function for W3 for a given sample size
#'
#' Generally, the weight W3 depends on the sample size and the shape parameter.
#' The sample size of a group is fixed. Hence, we return a function that returns
#' W3 for provided shape parameter as argument. If the sample size occured
#' during the Monte-Carlo simulation study the coefficients of generalized
#' logistic curve are directly returned. Otherwise a natural cubic spline is fit
#' on the fly. This guarantees that the spline function of the R-version of the
#' current user is used. Drawback is that performance is maybe not optimal
#' (`w3FFint` is not precompiled but run by the [objFunFactory()] once per
#' group)
#' @param nObs sample size for which to return the W3-function
#' @returns W3-function for the given sample size. The function returns the W3 weight for the given shape
w3FFint <- function(nObs) {

  if (missing(nObs) || length(nObs) != 1L || !is.numeric(nObs) || !is.finite(nObs)) {
    stop("Please provide the single number of observations within group!", call. = FALSE)
  }

  # catch all for n = 1 (or even n negative)
  if (nObs < 2L) return(function(k) 1)

  # get coefficients for a Richards' generalized logistic function
  # if we have fit the parameter nObs directly we use this Richards fit
  #+otherwise, we rely on the spline approximation for each parameter intrapolating the given nObs
  W3richCoef <- MLEw_approx[["coef"]][["W3_richards"]]
  approx_W3_names <- c("A", "K", "Q", "B", "nu")
  stopifnot(is.data.frame(W3richCoef), all(c("nObs", approx_W3_names) %in% names(W3richCoef)))
  approx_W3_ind <- which(W3richCoef$nObs == {nObs})

  # check for match in W3_richards
  approx_W3_coefs <- if (length(approx_W3_ind) == 1L) {
    W3richCoef[approx_W3_ind, approx_W3_names]
  } else {
    # nObs was not in MC-sim for W3, use interpolation per Richards coefficient
    list(
      A = stats::spline(x = W3richCoef$nObs,
                        y = W3richCoef$A,
                        method = "natural", xout = {nObs})$y,
      K = stats::spline(x = W3richCoef$nObs,
                        y = W3richCoef$K,
                        method = "natural", xout = {nObs})$y,
      Q = stats::spline(x = W3richCoef$nObs,
                        y = W3richCoef$Q,
                        method = "natural", xout = {nObs})$y,
      B = stats::spline(x = W3richCoef$nObs,
                        y = W3richCoef$B,
                        method = "natural", xout = {nObs})$y,
      nu = stats::spline(x = W3richCoef$nObs,
                         y = W3richCoef$nu,
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
    exp(expm1(MLEw_approx$fun$genLogisticF(theta = approx_W3_coefs,
                                           xVal = -log(k))))

    # evalq(expr = A + (K - A) / (1 + Q * k**-B)**(1/nu),
    #       envir = as.list(approx_W3_coefs),
    #       enclos = rlang::current_env())
  }#fn

}#fn w3FFint
