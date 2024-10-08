# mkuhn, 2021-05-25
# problematic data fitting,
# comes from 3rd iteration in:
# set.seed(12)
#
# te_H0 <- future.apply::future_replicate(n = 10, expr = {
#   x <- 10 + rexp(23, rate = .2)
#   y <- 10 + rexp(25, rate = .1)
#
#   te_diff <- test_delay_diff(x = x, y = y, distribution = "exp", param = "delay", R = 600)
# }, simplify = FALSE, future.seed = TRUE)

library(incubate)

distribution = "expo"
param <- "delay"

x <- c(11.9626910821535, 20.9178791996845, 10.29089452466, 18.1200836687804,
       14.442985006398, 10.5530852288939, 14.4648662160844, 18.4570179409803,
       12.7243826235645, 11.5971722337417, 20.7013441484499, 19.0127636208655,
       16.7845788968489, 15.4034043427319, 14.5589126148034, 16.2521717810875,
       13.7549493957888, 15.364592699124, 11.6033756686375, 12.4792049615644,
       21.500390683602, 10.3356814872802, 10.5235869324863)

y <- c(27.5840741573813, 10.2312779319737, 17.7072776498797, 17.3926557596297,
       12.5368687696755, 17.3279593787612, 23.9633602011729, 12.7061508976514,
       10.7451286529783, 11.7158383736387, 16.0329165263101, 13.2301631057635,
       37.3267119178984, 11.581134502776, 14.0266492217779, 29.6148683308604,
       11.281704692845, 18.5797028384218, 10.6765562470765, 15.9468249417841,
       28.6504221512916, 15.2674137521535, 48.1991241239864, 17.5364672760935,
       28.4304038613083)


fit0 <- delay_model(x = x, y = y, distribution = distribution, bind = param)
fit0 <- delay_model(x = x - 10, y = y - 10, distribution = distribution, bind = param)
fit1 <- delay_model(x = x, y = y, distribution = distribution)
# this currently works (incubate, 2021-06-04)
te_diff <- NULL
te_diff <- test_delay_diff(x = x, y = y, distribution = "exp", param = "delay", R = 100)
# problematic are small delays or large delays
te_diff <- test_delay_diff(x = x - 10, y = y - 10, distribution = "exp", param = "delay", R = 100)
te_diff <- test_delay_diff(x = x + 100, y = y + 100, distribution = "exp", param = "delay", R = 100)
#XXX find out why?!

plot(fit0)
plot(fit1)

plot(te_diff)
