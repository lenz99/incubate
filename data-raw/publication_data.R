# Apparently artificial example from
# Howard Rockette, "Maximum Likelihood Estimation with the Weibull Model"  (1974)
publication_examples <- c(3.1, 4.6, 5.6, 6.8)


# maximum flood level (in millions of cubic feet per second)
# for the Susquehanna River of Harrisburg (Pennsylvania, USA) over 20 4-year periods
# from Dumonceaux and Antle (1973)
# as cited by Cheng (1982)
susquehanna <- c(
  0.654, 0.613, 0.315, 0.449, 0.297,
  0.402, 0.379, 0.423, 0.379, 0.3235,
  0.269, 0.740, 0.418, 0.412, 0.494,
  0.416, 0.338, 0.392, 0.484, 0.265
)

# Bearing fatigue data originally reported by
# McCool, J.I., 1974. Inferential techniques for Weibull populations. Technical Report TR 74-0180, Wright Patterson Air Force Base, Ohio.
# data taken from Nagatsuka et al., "A consistent method of estimation for the three-parameter Weibull distribution", 2012
fatigue <- c(152.7, 172.0, 172.5, 173.3, 193.0,
             204.7, 216.5, 234.9, 262.6, 422.6)


# beach pollution
# from Steen & Stickler (1976)
# as cited by Cheng (1982)
pollution <- c(
  1364, 2154, 2236, 2518, 2527,
  2600, 3009, 3045, 4109, 5500,
  5800, 7200, 8400, 8400, 8900,
  11500, 12700, 15300, 18300, 20400
)


usethis::use_data(publication_examples, fatigue, susquehanna, pollution, overwrite = TRUE)
