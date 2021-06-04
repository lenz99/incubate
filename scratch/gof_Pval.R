# benchmark different statistical tests in the delay-test
#
# test delay parameter

library(incubate)

library(dplyr)
library(purrr)
library(tidyr)

library(future.callr)
library(future.apply)

future::plan(strategy = future.callr::callr, workers = 7L)

set.seed(12)

# same delay
te_H0 <- future.apply::future_replicate(n = 60, expr = {
  x <- 10 + rexp(23, rate = .2)
  y <- 10 + rexp(25, rate = .1)

  te_diff <- NULL
  try(expr = {te_diff <- test_delay_diff(x = x, y = y, distribution = "exp", param = "delay", R = 200)},
      silent = TRUE)

  te_diff
}, simplify = FALSE, future.seed = TRUE)

# drop NULLs
te_H0 <- purrr::compact(te_H0)

extract_Pvals <- function(testList){
  # drop NULLs
  testList <- purrr::compact(testList)

  tibble(
    boot = purrr::map_dbl(testList, ~ chuck(., "P", "boot")),
    gof_pearson = purrr::map_dbl(testList, ~ chuck(., "P", "gof_pearson")),
    gof_ad = purrr::map_dbl(testList, ~ chuck(., "P", "gof_ad")),
    lr = purrr::map_dbl(testList, ~ chuck(., "P", "lr")),
    lr_pp = purrr::map_dbl(testList, ~ chuck(., "P", "lr_pp"))
  ) %>% pivot_longer(cols = everything(), names_to = "method", values_to = "P")
}

ggplot(extract_Pvals(te_H0), mapping = aes(x = P, col = method)) +
  geom_freqpoly(bins = 21) +
  xlim(c(0,1)) +
  labs(title = "P-values under H0: no difference in delay", x = "P-value")


# delay diff 2 (10 vs 12)
te_H1 <- future.apply::future_replicate(n = 300, expr = {
  x <- 10 + rexp(20, rate = .1)
  y <- 12 + rexp(25, rate = .1)

  te_diff <- NULL
  try(expr = {te_diff <- test_delay_diff(x = x, y = y, distribution = "exp", param = "delay", R = 400)},
      silent = TRUE)
  te_diff
}, simplify = FALSE, future.seed = TRUE)

te_H1 <- purrr::compact(te_H1)

ggplot(extract_Pvals(te_H1), mapping = aes(x = P, col = method)) +
  geom_freqpoly(bins = 21) +
  xlim(c(0,1)) +
  labs(title = "P-values under a difference of 2 time units in delay, with rate λ~1~ = .2 and λ~2~ = .1", x = "P-value")


# teardown
future::plan(strategy = future::sequential)
