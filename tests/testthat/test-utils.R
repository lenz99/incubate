# mkuhn, 2021-10-11
# test utility functions of this package

test_that("Internal package data (related to weights functions for MLEw)", {
  # load internal package data
  FNAME_SYSD <- xfun::magic_path("sysdata.rda")
  stopifnot(file.exists(FNAME_SYSD))
  load(FNAME_SYSD)
  expect_true(exists(".MLEw_approx"))
  expect_named(.MLEw_approx, expected = c("MCsim", "coef", "fun"))
  expect_named(.MLEw_approx[["MCsim"]], expected = c("nObs", "W1", "W2"))
  expect_named(.MLEw_approx[["coef"]], expected = c("W2", "W3_richards"))
  expect_named(.MLEw_approx[["fun"]], expected = c("genLogisticF", "genLogisticJ", "w1F", "w2F", "w3FF"))
  # all functions in .MLEw_approx[["fun"]]
  purrr::walk(.x = names(.MLEw_approx[["fun"]]),
              .f = \(nam) expect_type(.MLEw_approx[["fun"]][[nam]], "closure"))

  w1F <- .MLEw_approx[["fun"]][["w1F"]]
  w2F <- .MLEw_approx[["fun"]][["w2F"]]
  w3FF <- .MLEw_approx[["fun"]][["w3FF"]]

  expect_length(w1F(seq_len(3)), n = 3)
  expect_type(w1F(5), type = "double")
  expect_equal(w1F(1), expected = log(2))
  expect_equal(w1F(-1), expected = log(2))
  expect_equal(w1F(1.5), expected = log(2))
  # W1 is (mostly) monotonically increasing
  expect_gte(min(diff(w1F(seq_len(40)))), expected = 0)

  # Cousineau: Nearly unbiased.. (Table 2)
  expect_equal(w1F(6), expected = 0.945, tolerance = 1e-3)
  expect_equal(w1F(14), expected = 0.976, tolerance = 1e-3)

  expect_length(w2F(seq_len(3)), n = 3)
  expect_type(w2F(5), type = "double")
  expect_equal(w2F(1), expected = 0)
  expect_equal(w2F(-1), expected = 0)
  expect_equal(w2F(1.5), expected = 0)
  # W2 is monotonically increasing
  expect_gte(min(diff(w2F(seq_len(100)))), expected = 0)

  # Cousineau: Nearly unbiased.. (Table 3)
  expect_equal(w2F(8), expected = 0.817, tolerance = 1e-3)
  expect_equal(w2F(15), expected = 0.902, tolerance = 1e-3)

  # Cousineau: Nearly unbaised.. (Table 4)
  # agreement here is not as high: W3 is more challenging, in particular for small shape
  #+our median is based on higher sample size in our MC-sim
  shapes <- seq.int(0.5, 2.5, by=.5)
  expect_equal(w3FF(6)(shapes),  expected = c( 5.631, 2.808, 2.004, 1.669, 1.492), tolerance = 5e-2)
  expect_equal(w3FF(11)(shapes), expected = c( 9.319, 3.462, 2.207, 1.774, 1.555), tolerance = 4e-2)
  expect_equal(w3FF(12)(shapes), expected = c(10.051, 3.560, 2.239, 1.782, 1.565), tolerance = 3e-2)
  expect_equal(w3FF(16)(shapes), expected = c(12.743, 3.854, 2.324, 1.820, 1.586), tolerance = 2e-2)

  # W3 for neighbouring nObs and some shapes
  # we use also fractional nObs to test if the spline interpolation of coefficients works properly
  w3Ex_mat <- c(15.999, 16, 16.1, 16.11, 16.25, 16.3, 16.5, 16.6, 16.7, 16.9, 17, 17.1, 17.2, 17.3, 17.5, 18, 19.1,
                50, 50.01, 50.1, 50.2, 51, 52, 53, 54, 55, 55.1, 55.11, 55.13, 59, 60, 61, 61.01, 61.02, 61.5, 61.6, 61.9, 61.99, 62) |>
    sort.int() |>
    purrr::map(.f = \(n_) w3FF(n_)(shapes)) |>
    unlist() |> matrix(ncol = length(shapes), byrow = TRUE,
                       dimnames = list(list(), shape = paste0("k=", shapes)))
  # W3 decreases as function of shape (for given n)
  expect_lte(w3Ex_mat |>
               apply(MARGIN = 1, FUN = diff) |>
               # diffs between shapes are put in columns: hence, next apply with MARGIN=2
               apply(MARGIN = 2, FUN = max, simplify = TRUE) |>
               max(), expected = 0)

  # W3 increases as function of n (for given shape)
  expect_gte(w3Ex_mat |>
    apply(MARGIN = 2, FUN = diff) |>
    apply(MARGIN = 2, FUN = min, simplify = TRUE) |>
    min(), expected = 0)

})

test_that("Estimate rounding error from sample", {
  set.seed(1234)
  # some random data (around 0)
  obsList <- list(obs1 = rnorm(31L),
                  obs2 = rt(37, df = 3),
                  obs3 = sqrt(rpois(51, lambda = 11)))


  # cross observation vectors and rounding digits 0:3
  purrr::walk2(.x = rep(obsList, 4),
               # rounding digits
               .y = rep(0:3, each = length(obsList)),
               .f = ~ expect_identical(estimRoundingError(round(.x, .y)), expected = 10**-.y))

  # rounds everything to 0
  expect_identical(estimRoundingError(round(obsList$obs1, -2)), expected = 1)
  expect_identical(estimRoundingError(round(obsList$obs2, -1)), expected = 1)
  expect_identical(estimRoundingError(round(obsList$obs3, -1)), expected = 1)


  # zeros at the end
  expect_identical(estimRoundingError(round(obsList$obs1, 3) * 10000), expected = 10)
  expect_identical(estimRoundingError(round(obsList$obs2, 2) * 1000), expected = 10)
  expect_identical(estimRoundingError(round(obsList$obs3, 1) * 100), expected = 10)
  expect_identical(estimRoundingError(round(obsList$obs1, 0) * 10), expected = 10)

  expect_identical(estimRoundingError(round(obsList$obs1, 1) * 1000), expected = 100)
  expect_identical(estimRoundingError(round(obsList$obs2, 1) * 1000), expected = 100)
  expect_identical(estimRoundingError(round(obsList$obs3, 1) * 1000), expected = 100)

  # exceeding the specified precision (e.g. here max(roundDigits) is 5) we expect to have one more
  expect_identical(estimRoundingError(obsList$obs1, roundDigits = -2:5), expected = 1e-6)
  expect_identical(estimRoundingError(obsList$obs2, roundDigits = -13:5), expected = 1e-6)
  expect_identical(estimRoundingError(obsList$obs1, roundDigits = 4:8), expected = 1e-9)
  expect_identical(estimRoundingError(obsList$obs3, roundDigits = -5:6), expected = 1e-7)
  expect_identical(estimRoundingError(obsList$obs3, roundDigits = -3:7), expected = 1e-8)
  expect_identical(estimRoundingError(obsList$obs3, roundDigits = -1:8), expected = 1e-9)

  # exceeding the specified precision on the negative side
  expect_identical(estimRoundingError(round(obsList$obs3,0)*1000, roundDigits = -2:5), expected = 10**3)
})


test_that("Ties in data", {

  set.seed(2023-04-17)
  # draw random data
  x <- sqrt(5 + stats::rpois(n = 17L, lambda = 9))
  # add ties
  x <- sort(sample(x = x, size = length(x)+1, replace = TRUE))
  xs <- sort(survival::Surv(time = x, event = sample(x = c(0, 1, 1), size = length(x), replace = TRUE), type = "right"))

  expect_error(objFunFactory(x = x,
                             distO = buildDist("exponential"),
                             control = buildControl(verbose = 0, profiled = FALSE, pen_shape = FALSE, ties = "error")),
               regexp = "ties")

  # objFunEqui1 <- objFunFactory(x = x, control = buildControl(ties = "equi"))
  # x_pp <- rlang::env_get(rlang::fn_env(objFunEqui1), nm = "x")
  # # there are no duplicates any more!
  # expect_false(any(duplicated(x_pp)))
  # #waldo::compare(x, y=rlang::env_get(rlang::fn_env(objFunEqui1), nm = "x"))
  # # deviations through tie-break are below and above tie and sum to zer0
  # expect_identical(sum(x_pp - x), expected = 0)

  # objFunEqui2 <- objFunFactory(x = xs, control = buildControl(ties = "equi"))
  # xs_pp <- rlang::env_get(rlang::fn_env(objFunEqui2), nm = "x")
  # # there are no duplicates any more!
  # #+ for array/matrix: it means no duplicate rows, here: no time + status duplicates!
  # #+ this means, we allow for same times that are once as observed time and once as a (right-) censoring time
  # expect_false(any(duplicated(xs_pp)))
})
