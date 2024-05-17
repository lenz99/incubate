# mkuhn, 2023-04-11
# adds MLE-weights table as internal data to package
#
# the MLE-weights are established through simulation, see inst/scripts/simul_MLEweights.R
# There, the approximating functions are built.
# Code to explore which are good/best approximations W1, W2 and W3 are in scratch/MLEw_weights2.R.
# The package incubate makes use of these functions in .MLEw_approx[["fun"]] etc.
###

# init --------------------------------------------------------------------

library("usethis")

FNAME <- "MLEw_weights.RData"
stopifnot(file.exists(FNAME))
(load(FNAME))

stopifnot(is.list(.MLEw_mcs), is.list(.MLEw_approx))
stopifnot(identical(names(.MLEw_mcs), c("W12","W3", "settings")),
          is.data.frame(.MLEw_mcs$W12), is.data.frame(.MLEw_mcs$W3))
stopifnot(identical(names(.MLEw_approx), c("coef", "fun")))

# save as internal data ---------------------------------------------------

usethis::use_data(.MLEw_mcs, .MLEw_approx, internal = TRUE, overwrite = TRUE)


message("~~Fine~~")

#q(save = "no")
