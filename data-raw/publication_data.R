# Bearing fatigue data originally reported by
# McCool, J.I., 1974. Inferential techniques for Weibull populations. Technical Report TR 74-0180, Wright Patterson Air Force Base, Ohio.
# data taken from Nagatsuka et al., "A consistent method of estimation for the three-parameter Weibull distribution", 2012

fatigue <- c(152.7, 172.0, 172.5, 173.3, 193.0,
             204.7, 216.5, 234.9, 262.6, 422.6)




# Apparently artificial example from
# Howard Rockette, "Maximum Likelihood Estimation with the Weibull Model"  (1974)
rockette <- c(3.1, 4.6, 5.6, 6.8)


usethis::use_data(fatigue, rockette, overwrite = TRUE)
