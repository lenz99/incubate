# mkuhn, 2021-09-10:
# read in data from Stankovic 2018 publication,
# taken from Figure 6J and Figure 6K

library(readr)
library(dplyr)

try(setwd('data-raw/'), silent = FALSE)

#' Prepare survival data for analysis with MSE-method.
#' The data was digitized from the survival plots in the publication of Stankovic et al (2018)
#' @param fileName file name of CSV-file within folder 'stankovic_2018'
prepSurvData <- function(fileName) {
  fnPath <- file.path(fileName)
  stopifnot( file.exists(fnPath) )

  dat <- readr::read_csv(fnPath, comment = "#", show_col_types = FALSE) %>%
    mutate(Time = round(Time), r = round(n * Survival)) %>%
    group_by(Group) %>%
    # d: number of events
    mutate(d=lag(r, default = first(n))-r) %>%
    ungroup

  with(dat,
       bind_cols(Time=rep.int(Time, times = d),
                 Status = 1L, # all observations are events
                 Group = rep.int(Group, times = d),
                 Colour = rep.int(Colour, times = d)))
}

fig6J <- prepSurvData("stankovic_fig6J_U87.csv")
fig6K <- prepSurvData('stankovic_fig6K_U87.csv')

library(survival)
# visual test
plot(survfit(Surv(Time, Status) ~ Group, data = fig6J), col = 1:4, main = 'Figure 6J')
plot(survfit(Surv(Time, Status) ~ Group, data = fig6K), col = 1:3, main = 'Figure 6K')

# check P-values from figure
# fig1J
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6J, subset = Group %in% c('CTRL', 'Anti-E7'))$chisq, df = 1, lower.tail = F)
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6J, subset = Group %in% c('CTRL', 'Anti-VEGF'))$chisq, df = 1, lower.tail = F)
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6J, subset = Group %in% c('CTRL', 'Combo'))$chisq, df = 1, lower.tail = F)
# p=0.293. MS Fig1J writes p=0.2!
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6J, subset = Group %in% c('Anti-VEGF', 'Combo'))$chisq, df = 1, lower.tail = F)

# fig1K: P-values are roughly OK, but
#+rounding is incorrect in MS, P-values in MS are generally rounded downwards (to 0).
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6K, subset = Group != 'Combo')$chisq, df = 1, lower.tail = F)
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6K, subset = Group != 'Anti-VEGF')$chisq, df = 1, lower.tail = F)
pchisq(survdiff(Surv(Time, Status) ~ Group, data = fig6K, subset = Group != 'TMD')$chisq, df = 1, lower.tail = F)


stankovic <- dplyr::bind_rows(fig6J=fig6J, fig6K=fig6K, .id = 'figure')
# save as one dataframe, including column 'figure' to describe source
#saveRDS(stankovic, file = 'stankovic_fig6JK_U87.rds')

#usethis::use_data(stankovic, overwrite = FALSE)
