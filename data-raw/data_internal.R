# mkuhn, 2023-04-11
# add MLE-weights table as internal data to package
# the MLE-weights are established through simulation, see inst/scripts/simul_MLEweights.R
# There, the approximating functions are build, code to explore which are good/best approximations W1, W2 and W3
# are in scratch/MLEw_weights2.R. The package incubate makes use of these functions in .MLEw_approx[["fun"]] etc.
###

# init --------------------------------------------------------------------

library("usethis")
library("dplyr")
library("ggplot2")
library("patchwork")
library("nlsr")
library("gslnls")


#.MLEw_weights <- readRDS("MLEw_weights.rds")
load("MLEw_weights.RData")

stopifnot( is.list(.MLEw_mcs), is.list(.MLEw_approx),
           identical(names(.MLEw_mcs), c("W12","W3", "setting")), is.data.frame(.MLEw_mcs$W3),
           identical(names(.MLEw_approx), c("coef", "fun")) )

# save as internal data ---------------------------------------------------

usethis::use_data(.MLEw_mcs, .MLEw_approx, internal = TRUE, overwrite = TRUE)


message("~~Fine~~")

q(save = "no")
