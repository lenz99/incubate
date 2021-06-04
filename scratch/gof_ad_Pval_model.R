# critical values for goodness of fit tests
# model derived from given critical values and significance levels
# for Anderson Darling test with exponential and weibull distributions

library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)



# exponential ----

# cf «GOF-techniques» book
# P-value for A2_mod based on upper tail percentage points (table 4.14, p.138)
A2_mod_crit <- tribble(
  ~sig_level, ~crit_val,
  .25, .736,
  .15, .916,
  .10, 1.062,
  .05, 1.321,
  .025, 1.591,
  .01, 1.959)

ggplot(A2_mod_crit, mapping = aes(x = crit_val, y=sig_level)) +
  scale_y_log10() +
  geom_point() +
  geom_smooth(method = lm, formula = y ~x, se = F, size = .2, linetype = "dotted", fullrange = TRUE) +
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = F, size = .2, col = "red", fullrange = TRUE) +
  expand_limits(x=c(.5, 2.5))


ggplot(A2_mod_crit, mapping = aes(x = sqrt(crit_val), y = log(sig_level/(1-sig_level)))) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~x, se = F, size = .2, linetype = "dotted", fullrange = TRUE) +
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = F, size = .2, col = "red", linetype = "dashed", fullrange = TRUE) +
  #geom_smooth(method = lm, formula = y ~ poly(x,3), se = F, size = .2, col = "blue", fullrange = TRUE) +
  expand_limits(x=c(.01, 2.5))

summary(lm(log(sig_level) ~ crit_val, data = A2_mod_crit))
summary(lm(log(sig_level / (1-sig_level)) ~ sqrt(crit_val), data = A2_mod_crit))
summary(lm(log(sig_level / (1-sig_level)) ~ crit_val + sqrt(crit_val), data = A2_mod_crit))


# weibull -----

# read in critical values for Weibull GOF
# cf Lockhart, 1994
crit_weib <- read.table(file = "gof_ad_weibull.txt", sep = " ", header = T) %>%
  pivot_longer(cols = starts_with("crit"),
               names_to = "sig_level",
               names_prefix = "crit[_]", names_transform = list(sig_level = as.numeric),
               values_to = "crit_val") %>%
  mutate(sig_level = 1L - sig_level)

ggplot(crit_weib, mapping = aes(x = sqrt(crit_val), y = log(sig_level/(1-sig_level)))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F, size = .2, linetype = "dotted", fullrange = TRUE) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F, size = .2, col = "red", fullrange = TRUE) +
  facet_wrap(vars(c), labeller = label_both) +
  expand_limits(x=c(.4, 1.2))

fm_cat <- lm(log(sig_level/(1-sig_level)) ~ sqrt(crit_val) * factor(c), data = crit_weib)
summary(fm_cat)

fm_c <- lm(log(sig_level/(1-sig_level)) ~ sqrt(crit_val) * c, data = crit_weib)
summary(fm_c)
fm_c2 <- lm(log(sig_level/(1-sig_level)) ~ crit_val + sqrt(crit_val) * c, data = crit_weib)
summary(fm_c2)
fm_c3 <- lm(log(sig_level/(1-sig_level)) ~ (crit_val + sqrt(crit_val)) * c, data = crit_weib)
summary(fm_c3)
# fm_cat and fm_c both have R^2 of 0.999
# fm_c allows for easier extrapolation of c
# although fm_cat comes out as winner in F-test
anova(fm_c, fm_cat)

ggplot(crit_weib, mapping = aes(x = log(crit_val), y = sig_level)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x,2), se = F, size = .2, linetype = "dotted") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3), se = F, size = .2, col = "red") +
  facet_wrap(vars(c), labeller = label_both)



