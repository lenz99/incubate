# mkuhn, 2025-04-24:
# read in data from Long 2017 publication:
# "Adjuvant Dabrafenib plus Trametinib in Stage III BRAF-Mutated Melanoma"
# relapse free survival in melanoma patients
# taken from Kaplan-Meier Survival Plot (Fig1A)
# sent by Sean Devlin (sedevlin@gmail.com, email 2025-04-18)
# digitizing technique as described by Guyot et al., 2012; Satagopan et al., 2017

library("readr")
library("dplyr")
library("survival")

try(setwd('data-raw/'), silent = FALSE)

FNAME <- "digitized_data_Long_2017.csv"
stopifnot(file.exists(FNAME))
long2017 <- readr::read_csv(file = FNAME,
                            col_names = TRUE, col_types = cols(
                              ID = col_character(),
                              time = col_double(),
                              event = col_double(),
                              trtmt_nmbr = col_skip(),
                              trtmt = col_character()
                            )) |>
  dplyr::mutate(trtmt = dplyr::case_match(trtmt,
                                          "dabrafenib" ~ "Dabrafenib+Trametinib",
                                          .default = trtmt, .ptype = "X")) |>
  dplyr::rename(status = event)

# check data
kmfit_long2017 <- survfit(Surv(time, status) ~ trtmt, data = long2017)
# KM-survival plot
plot(kmfit_long2017, mark.time = TRUE, main = "Relapse Free Survival",
     sub = "Long (2017):\nAdjuvant Dabrafenib plus Trametinib\nin Stage III BRAF-Mutated Melanoma")

# Long (2017) reports:
# The estimated rates of relapse-free survival were
# combination-therapy group:
# 88% at 1 year,
# 67% at 2 years, and
# 58% at 3 years
# in the placebo group:
# rates of 56%, 44%, and 39%, respectively.
# in fig1:
# Hazard ratio for relapse or death, HR=0.47 (95% CI, 0.39–0.58) P<0.001 by stratified log- rank test
# As of the data cutoff at a median of 2.8 years of follow‐up,
# disease recurrence or death had been reported in 166 of 438 patients (38%) in the combination‐therapy group
# and in 248 of 432 patients (57%) in the placebo group.
summary(kmfit_long2017, times = c(1, 2, 2.8,3)*12, scale = 12)


# save data in package
usethis::use_data(long2017, overwrite = TRUE)
