# mkuhn, 2023-04-11
# add MLE-weights table as internal data to package
# the MLE-weights are established through simulation, see inst/scripts/simul_MLEweights.R
###

library("usethis")

.MLEw_weights <- readRDS("inst/scripts/MLEw_weights.rds")

stopifnot( is.list(.MLEw_weights), is.data.frame(.MLEw_weights$W3) )



library("ggplot2")
library("patchwork")
library("nlsr")
library("gslnls")


# approximate W12 ----------------------------------------------------------

# currently still exploring approximation fits
stopifnot(rlang::is_interactive())


W12 <- .MLEw_weights$W12 %>%
  dplyr::mutate(lnObs = log(nObs),
                # approximation for median of gamma(n, 1/n)
                #+using Wilson-Hilferty transformation (see <https://en.wikipedia.org/wiki/Gamma_distribution>)
                W1pred = (1 - 1 / (9 * nObs))^3,
                W1resid = W1 - W1pred)


ggplot(W12, mapping = aes(x = nObs, y = W1)) +
  geom_point() +
  scale_x_log10() +
  labs(title = "W1 as median") +
  geom_line(mapping = aes(y = W1pred), col = "blue")  |

  ggplot(W12, mapping = aes(x = nObs, y = W1resid)) +
  geom_point() + geom_line() +
  scale_x_log10() +
  labs(title = "Deviation")

# SSasymp on log(nObs):
# W2 = 1 + (R0 - 1) * n**-r
fm_W2_A <- nls(W2 ~ SSasymp(input = lnObs, Asym = 1, R0, lrc), start = list(R0 = 0, lrc = -.02),
               # (0,0) point is taken out (as it does not follow the overall functional pattern)
               data = W12, subset = -1)
summary(fm_W2_A)

# equivalent code with nlsr::nlxb
fm_W2_B <- nlxb(W2 ~ 1 + (R0-1) * nObs**-exp(lrc), start = list(R0=0, lrc = -.02),
                data = W12[-1,])
print(fm_W2_B)

fm_W2_B <- nlxb(W2 ~ 1 + (R0-1) * nObs**-exp(lrc),
                start = list(R0=0, lrc = -.02),
                data = W12[-1,])
print(fm_W2_B)

# rate directly on original scale (and hence with lower bound = 0). same result
fm_W2_C <- nlxb(W2 ~ 1 + (R0-1) * nObs**-rc,
                start = list(R0=0, rc = exp(-.02)), lower = c(-Inf, 0),
                data = W12[-1,])
print(fm_W2_C)

fm_W2_D <- gsl_nls(W2 ~ 1 + (R0-1) * nObs**-exp(lrc),
                   start = list(R0=0, lrc = -.02),
                   data = W12[-1,])
summary(fm_W2_D)

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



# approximation W3 --------------------------------------------------------


W3 <- .MLEw_weights$W3 %>%
  dplyr::mutate(lnObs = log(nObs),
                shapeF = factor(shape))

W3XS <- W3 %>% filter(shape <= .25) %>% droplevels()
W3S <- W3 %>% filter(between(shape, .251, .75)) %>% droplevels()
W3M <- W3 %>% filter(between(shape, .751, 1))
W3L <- W3 %>% filter(shape > 1) %>% droplevels()

# facet per n, ##W3 %>% dplyr::filter(between(nObs, 2, 13))
ggplot(W3 %>% filter(nObs > 1), mapping = aes(x = shape, y = W3, col = ordered(nObs))) +
  geom_point() + geom_line() +
  #facet_wrap(facets = cut_number(W3$nObs[W3$nObs>1], n = 9), scales = "free_y") +
  scale_x_log10()
#scale_y_log10()

ggplot(W3 %>% filter(nObs >= 1000), mapping = aes(x = shape, y = W3, col = ordered(nObs))) +
  geom_point() + geom_line() +
  scale_x_log10() +
  scale_y_log10()

# facet per shape
ggplot(W3L, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
  geom_point() + geom_line() |
  #scale_x_log10() |
  #scale_y_log10() |

  ggplot(W3M, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
  geom_point() + geom_line() |
  #scale_x_log10() |
  #scale_y_log10() |

  ggplot(W3S, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
  geom_point() + geom_line() |
  #scale_x_log10()
  #scale_y_log10()

  ggplot(W3XS, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
  geom_point() + geom_line()
#scale_x_log10()
#scale_y_log10()

# logistic fit as function of shape k
W3logi2 <- W3 %>%
  dplyr::filter(nObs > 1) %>%
  dplyr::mutate(lshape = log(shape))

# problematic: small nObs OR shape=1
W3logi10 <- W3logi2 %>%
  dplyr::filter(nObs >= 10)


# 5-parametric logistic model
fm_W3_A <- gsl_nls(W3 ~ D0 + D1 * nObs + D2 * lnObs + (A0 + A1 * nObs - D0 - D1 * nObs - D2 * lnObs) / (1 + exp((B0 + B1 * nObs + B2 * lnObs - lshape) / (C0 + C1 * lnObs)))**(E0 + E1 * nObs + E2 * lnObs),
                   start = list(A0 = 0, A1 = 1, B0 = 0.5, B1 = 0, B2 = 0, C0 = -1, C1 = 0, D0 = 0, D1 = 0, D2 = 0, E0 = 1, E1 = 0, E2 = 0),
                   data = W3logi2) #, subset = nObs >= 10)
summary(fm_W3_A)

fm_W3_B <- gsl_nls(W3 ~ D0 + D1 * nObs + D2 * lnObs + (A0 + A1 * nObs - D0 - D1 * nObs - D2 * lnObs) / (1 + exp((B0 + B1 * nObs + B2 * lnObs - lshape) / (C0 + C1 * lnObs)))**(E0 + E1 * nObs + E2 * lnObs),
                   start = list(A0 = 0, A1 = 1, B0 = 0.5, B1 = 0, B2 = 0, C0 = -1, C1 = 0, D0 = 0, D1 = 0, D2 = 0, E0 = 1, E1 = 0, E2 = 0),
                   data = W3logi10)
summary(fm_W3_B)

# Gompertz model
fm_W3_C <- gsl_nls(W3 ~ D0 + D1 * nObs + D2 * lnObs + (A0 + A1 * nObs - D0 - D1 * nObs - D2 * lnObs) * exp(-b2 * exp(-b3*shape)), data = W3logi2,
                   start = list(A0 = 0, A1 = 1, D0 = 0, D1 = 0, D2 = 0, b2 = .1, b3 = .1),
                   control = gsl_nls_control(maxiter = 101))
summary(fm_W3_C)

# generalized logistic function (Richard's curve)
fm_W3_D <- gsl_nls(W3 ~ A0 + A1 * nObs + (K0 + K1 * lnObs - A0 - A1 * nObs) / (1 + (Q0 + Q1 * lnObs) * exp(-(B1 * lnObs) * lshape))**(1/(nu0 + nu1 * lnObs)),
                   data = W3logi2, #jac = TRUE,
                   start = list(A0 = 1, A1 = 2, K0 = 1, K1 = 0,
                                Q0 = .001, Q1 = 0.001, B1 = 1, nu0 = .1, nu1 = .1), #2, nu1 = 0),
                   control = gsl_nls_control(maxiter = 1010))
summary(fm_W3_D)

W3logi2 <- W3logi2 %>%
  mutate(W3pred = predict(fm_W3_D),
         W3resid = W3 - W3pred)
W3logi10 <- W3logi10 %>%
  mutate(W3pred = predict(fm_W3_D),
         W3resid = W3 - W3pred)

ggplot(data = W3logi2, mapping = aes(x = shape, y = W3, col = ordered(nObs))) +
  geom_point() + geom_line() +
  geom_line(mapping = aes(y = W3pred, group = nObs), col = "darkred") +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 15)) +
  facet_wrap(facets = cut_number(W3logi2$nObs, n = 4), scales = "free_y") +
  guides(col = "none")

# relative deviation
ggplot(data = W3logi2, mapping = aes(x = shape, y = W3resid / W3, col = ordered(nObs))) +
  geom_point() + geom_hline(yintercept = 0, col = "darkgrey", linetype = "dashed") +
  scale_y_continuous(labels = scales::label_percent())

W3logi2 %>%
  dplyr::arrange(desc(abs(W3resid/W3)))


fm_W3_E1 <- gsl_nls(W3 ~ A + (K - A) / (1 + Q * exp(-B * lshape))**(1/nu),
                    data = W3logi2, subset = nObs == 1500, #jac = TRUE,
                    start = list(A = 2500, K = 1, Q = .005, B = 8, nu = 1.5),
                    control = gsl_nls_control(maxiter = 1010))
summary(fm_W3_E1)

nObs_vctr <- W3logi2 %>% distinct(nObs) %>% pull(nObs) #c(2:20, 25, 50, 75)
nObs_vctr <- nObs_vctr[1:30]
fm_W3_indiv <- purrr::map(.x = nObs_vctr,
                          .f = ~gsl_nls(W3 ~ A + (K - A) / (1 + Q * exp(-B * lshape))**(1/nu),
                                        data = W3logi2, subset = nObs == .x, #jac = TRUE,
                                        start = list(A = 2*.x, K = 1, Q = .05, B = log(.x+1), nu = log(.x+1)/3),
                                        control = gsl_nls_control(maxiter = 1010)))

coef_lowN <- purrr::map(fm_W3_indiv, .f = coef) %>%
  purrr::list_transpose(simplify = TRUE) %>%
  append(values = list(nObs=nObs_vctr), after = 0) %>%
  as.data.frame()

opar <- par(mfrow = c(2, 3))
plot(x = nObs_vctr, y = map_dbl(fm_W3_indiv, .f = sigma), main = "Residual std. deviation")
iwalk(coef_lowN[-1], .f = ~plot(x = nObs_vctr, y = .x, main = paste("parameter", .y)))
par(opar)

# modeling parameter dependence on nObs

# A -- upper asymptote
fm_coef_A <- lm(A ~ nObs, data = coef_lowN)
summary(fm_coef_A)
ggplot(broom::augment(fm_coef_A), mapping = aes(x = nObs, y = A)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred") +
  geom_abline(data = tibble::enframe(purrr::set_names(coef(fm_coef_A), nm = c("intercept", "slope"))) %>%
                tidyr::pivot_wider(),
              mapping = aes(intercept = intercept, slope = slope), col = "grey", linetype = "dashed")

# K -- lower asymptote
fm_coef_K1 <- lm(K ~ log(nObs), data = coef_lowN)
summary(fm_coef_K1)
plot(fm_coef_K1)
fm_coef_K2 <- lm(K ~ poly(log(nObs),2), data = coef_lowN)
summary(fm_coef_K2)
plot(fm_coef_K2)

ggplot(broom::augment(fm_coef_K1), mapping = aes(x = `log(nObs)`, y = K)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred") +
  geom_abline(data = tibble::enframe(purrr::set_names(coef(fm_coef_K1), nm = c("intercept", "slope"))) %>%
                tidyr::pivot_wider(),
              mapping = aes(intercept = intercept, slope = slope), col = "grey", linetype = "dashed")

ggplot(broom::augment(fm_coef_K2) %>% bind_cols(nObs = nObs_vctr),
       mapping = aes(x = nObs, y = K)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred", size = .33) + geom_line(mapping = aes(y = .fitted), col = "darkred")

# Q -- related to Y(0)
fm_coef_Q1 <- lm(Q ~ log(nObs), data = coef_lowN)
summary(fm_coef_Q1)

fm_coef_Q2 <- lm(Q ~ poly(log(nObs),2), data = coef_lowN)
summary(fm_coef_Q2)

ggplot(broom::augment(fm_coef_Q1), mapping = aes(x = `log(nObs)`, y = Q)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred", size = .5) +
  geom_abline(data = tibble::enframe(purrr::set_names(coef(fm_coef_Q1), nm = c("intercept", "slope"))) %>%
                tidyr::pivot_wider(),
              mapping = aes(intercept = intercept, slope = slope), col = "grey", linetype = "dashed")

# nObs = 2 and 3 are problematic (high residuals)
ggplot(broom::augment(fm_coef_Q2) %>% bind_cols(nObs = nObs_vctr),
       mapping = aes(x = nObs, y = Q)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred", size = .33) + geom_line(mapping = aes(y = .fitted), col = "darkred")


# B -- growth rate
fm_coef_B1 <- lm(B ~ log(nObs), data = coef_lowN)
summary(fm_coef_B1)
plot(fm_coef_B1)

fm_coef_B2 <- lm(B ~ poly(log(nObs),2), data = coef_lowN)
summary(fm_coef_B2)
plot(fm_coef_B2)

ggplot(broom::augment(fm_coef_B1), mapping = aes(x = `log(nObs)`, y = B)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred") +
  geom_abline(data = tibble::enframe(purrr::set_names(coef(fm_coef_B1), nm = c("intercept", "slope"))) %>%
                tidyr::pivot_wider(),
              mapping = aes(intercept = intercept, slope = slope), col = "grey", linetype = "dashed")

# nu -- maximal growth near which asymptote

fm_coef_nu1 <- lm(nu ~ log(nObs), data = coef_lowN)
summary(fm_coef_nu1)
plot(fm_coef_nu1)

ggplot(broom::augment(fm_coef_nu1), mapping = aes(x = `log(nObs)`, y = nu)) +
  geom_point() +
  geom_point(mapping = aes(y = .fitted), col = "darkred") +
  geom_abline(data = tibble::enframe(purrr::set_names(coef(fm_coef_nu1), nm = c("intercept", "slope"))) %>%
                tidyr::pivot_wider(),
              mapping = aes(intercept = intercept, slope = slope), col = "grey", linetype = "dashed")



W3logi2 %>% filter(nObs == 1500) %>%
  mutate(W3predi = predict(fm_W3_E1)) %>%
  ggplot(mapping = aes(x = shape, y = W3)) +
  geom_point() + geom_line() +
  geom_line(mapping = aes(y = W3pred, group = nObs), col = "darkred") +
  geom_line(mapping = aes(y = W3predi, group = nObs), col = "darkorange") +
  scale_x_log10()



# global fit with nls
fm_W3L_A <- nls(W3 ~ SSasymp(input = lnObs, Asym, R0 = 1, lrc), start = list(Asym = 10, lrc = 0),
                data = W3L, subset = shapeF == "1.25")
summary(fm_W3L_A)

# per group with nls
fm_W3L_B <- nls(W3 ~ SSasymp(input = lnObs, Asym[shapeF], R0 = 1, lrc[shapeF]),
                start = list(Asym = c(15, 5, 3, rep.int(2.3815, 12)), lrc = rep.int(-1.32, 15)),
                data = W3L)

fm_W3L_C <- nlxb(W3 ~ A0 + A1 * shape + (1-A0-A1 * shape) * nObs**-exp(lrc),
                 start = list(A0 = 10, A1 = -.2, lrc = -1),
                 data = W3L)
print(fm_W3L_C)

fm_W3L_D <- nlxb(W3 ~ shape / (shape-1) + (1-shape / (shape-1)) * nObs**-exp(lrc),
                 start = list(lrc = -1),
                 data = W3L)
print(fm_W3L_D)

fm_W3L_E <- nlxb(W3 ~ 1 + exp(A0 + A1 * shape + A2 * shape^2) - exp(A0 + A1 * shape + A2 * shape^2) * nObs**-exp(lrc),
                 start = list(A0 = 1, A1 = -1, A2 = 0, lrc = -1),
                 data = W3L)
print(fm_W3L_E)

fm_W3L_F <- nlxb(W3 ~ A0 + A1 / (shape+1) + (1-A0 - A1/ (shape+1)) * nObs**-exp(lrc),
                 start = list(A0 = 1, A1 = 2, lrc = -1),
                 data = W3L)
print(fm_W3L_F)

W3L <- W3L %>% mutate(W3pred = predict(fm_W3L_F, newdata = W3L))
ggplot(W3L, mapping = aes(x = nObs, y = W3, col = ordered(shape))) +
  geom_point() + geom_line() +
  geom_line(mapping = aes(y = W3pred, group = shape), col = "darkred")


W3L1.25 <- W3L %>% filter(near(shape, 1.25))
fm_W3L_C <- nls(W3 ~ SSasymp(input = lnObs, Asym, R0 = 1, lrc),
                start = list(Asym = 7, lrc = -1),
                data = W3L1.25)
summary(fm_W3L_C)
W3L1.25 <- W3L1.25 %>% mutate(W3predC = predict(fm_W3L_C))

ggplot(W3L1.25, mapping = aes(x = nObs, y = W3 - W3predC)) +
  geom_point()

W3L2 <- W3L %>% filter(near(shape, 2))
fm_W3L_D <- nls(W3 ~ SSasymp(input = lnObs, Asym, R0 = 1, lrc),
                start = list(Asym = 7, lrc = -1),
                data = W3L2)
summary(fm_W3L_D)
W3L2 <- W3L2 %>% mutate(W3predD = predict(fm_W3L_D))

ggplot(W3L2, mapping = aes(x = nObs)) +
  geom_point(mapping = aes(y = W3)) +
  geom_point(mapping = aes(y = W3predD), col = "red") +
  geom_line(mapping = aes(y = W3predD), col = "red") |
  ggplot(W3L2, mapping = aes(x = nObs, y = W3 - W3predD)) +
  geom_point()

# W3L linear models
fm_W3L_C <- lm(log(W3) ~ nObs + (lnObs + sqrt(lnObs) + log(lnObs+1)) * shapeF, data = W3L)
summary(fm_W3L_C)
drop1(fm_W3L_C, test = "F")

# XS
fm_W3XS_A <- lm(W3 ~ poly(nObs,2) + nObs * poly(shape,2) , data = W3XS)
summary(fm_W3XS_A)

W3XS <- W3XS %>%
  dplyr::mutate(W3predA = predict(fm_W3XS_A))

ggplot(W3XS, mapping = aes(x = W3, y = W3 - W3predA, col = ordered(shape))) +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed")
#ylim(-1, 1)

fm_W3S_A <- lm(W3 ~ poly(nObs,2) * shape, data = W3S)
summary(fm_W3S_A)

W3S <- W3S %>%
  dplyr::mutate(W3predA = predict(fm_W3S_A))

ggplot(W3S, mapping = aes(x = W3, y = W3 - W3predA, col = ordered(shape))) +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed")



# XXX add approximations to MLEw_weights list


# save as internal data ---------------------------------------------------

usethis::use_data(.MLEw_weights, internal = TRUE, overwrite = myOverwrite)


message("~~Fine~~")

