# mkuhn, 2021-10-11
# test utility functions of this package

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

  set.seed(20230417)
  # draw random data
  x <- sqrt(5 + stats::rpois(n = 17L, lambda = 9))
  # add ties
  x <- sort(sample(x = x, size = length(x)+1, replace = TRUE))

  xs <- sort(survival::Surv(time = x, event = sample(x = c(0, 1, 1), size = length(x), replace = TRUE), type = "right"))

  expect_error(objFunFactory(x = x, ties = "error"))

  objFunEqui1 <- objFunFactory(x = x, ties = "equi")
  x_pp <- rlang::env_get(rlang::fn_env(objFunEqui1), nm = "x")
  # there are no duplicates any more!
  expect_false(any(duplicated(x_pp)))
  #waldo::compare(x, y=rlang::env_get(rlang::fn_env(objFunEqui1), nm = "x"))
  # deviations through tie-break are below and above tie and sum to zer0
  expect_identical(sum(x_pp - x), expected = 0)

  objFunEqui2 <- objFunFactory(x = xs, ties = "equi")
  xs_pp <- rlang::env_get(rlang::fn_env(objFunEqui2), nm = "x")
  # there are no duplicates any more!
  #+ for array/matrix: it means no duplicate rows, here: no time + status duplicates!
  #+ this means, we allow for same times that are once as observed time and once as a (right-) censoring time
  expect_false(any(duplicated(xs_pp)))
})
