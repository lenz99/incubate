# mkuhn, 2021-10-13
# debug moran's GOF test

library(dplyr)
library(ggplot2)
library(incubate)
packageVersion('incubate')


# quantiles (cheng p3) ---------------------------------------------------------------

n <- 20
EUL_MAS <- -digamma(1L)
mo_m <- (n+1L) * (log(n+1) + EUL_MAS) - .5 - (12L*(n + 1L))**-1L
mo_v <- (n+1L) * (pi**2L / 6L - 1L) - .5 - (6L*(n + 1L))**-1L

C1 <- mo_m - sqrt(.5 * n * mo_v)
C2 <- sqrt(mo_v / (2L*n))

(C1 + C2 * qchisq(p = c(.9, .95, .99), df = n, lower.tail = TRUE)) / n




# delayed exponential -----------------------------------------------------

# from a call to 9 + rexp(13, rate = 0.5)
xx <- c(9.37584220431745, 9.43687826953828, 9.44079428166151, 9.63324003852904,
        9.76594032067806, 9.80526794679463, 9.90732720028609, 10.3573373407125,
        10.596041733315, 10.6229753642434, 11.1074129025543, 11.5750403608287,
        16.3800762999327)

fd_exp <- delay_model(xx, distribution = "expon")
coef_exp <- coef(fd_exp)

n <- length(xx)
k <- length(coef_exp)


cumDiffs <- diff(c(0, pexp_delayed(q = xx, delay = coef_exp[['delay']], rate = coef_exp[['rate']]), 1))
# no ties
stopifnot( length(which(dplyr::near(cumDiffs, 0))) == 0L )

# respect the machine's numerical lower limit
cumDiffs[which(cumDiffs < .Machine$double.eps)] <- .Machine$double.eps

mseCrit <- -mean(log(cumDiffs))
# matches the model best fit
stopifnot( dplyr::near(fd_exp$val, mseCrit) )


EUL_MAS <- -digamma(1L)
mo_m <- (n+1L) * (log(n) + EUL_MAS) - .5 - (12L*(n + 1L))**-1L
mo_v <- (n+1L) * (pi**2L / 6L - 1L) - .5 - (6L*(n + 1L))**-1L

C1 <- mo_m - sqrt(.5 * n * mo_v)
C2 <- sqrt(mo_v / (2L*n))

T_moran <- (mseCrit * (n+1L) + .5 * k - C1) / C2
pchisq(q = T_moran, df = n, lower.tail = FALSE)



# Normal example from Cheng ------------------------------------------------

h590 <- sort.int(c(27.55, 31.82, 33.74, 34.15, 35.32, 36.78,
29.89, 32.23, 33.74, 34.44, 35.44, 37.07,
30.07, 32.28, 33.86, 34.62, 35.61, 37.36,
30.65, 32.69, 33.86, 34.74, 35.61, 37.36,
31.23, 32.98, 33.86, 34.74, 35.73, 37.36,
31.53, 33.28, 34.15, 35.03, 35.90, 40.28,
31.53, 33.28, 34.15, 35.03, 36.20))

mean(h590); var(h590)
qplot(x=h590, binwidth = .3)

n.h590 <- length(h590)
k.h590 <- 2L

# break ties from incubate package
preprocess <- function(obs, verbose = 0, ties = c('equidist', 'random', 'density')) {
  if (is.null(obs)) return(NULL)
  ties <- match.arg(ties)
  if (ties == 'density' ) return(obs) ##|| ties == 'groupedML') # groupedML not implemented yet


  diffobs <- diff(obs)
  stopifnot( all(diffobs >= 0L) ) # i.e. sorted obs

  tiesDiffInd <- which(diffobs == 0L) # < .Machine$double.xmin

  if (length(tiesDiffInd)){
    #rl <- rle(diff(tiesDiffInd))
    if (verbose > 0L){
      cat(length(tiesDiffInd) + sum(diff(tiesDiffInd)>1L) + 1L, 'tied observations in',
          sum(diff(tiesDiffInd)>1L) + 1L, #length(which(rl$values > 1L))+1L,
          'group(s) within data vector.\n')
    }

    roundOffPrecision <- incubate:::estimRoundingError(obs, maxObs = 1000L)
    if (verbose > 0L){
      cat("Round-off error has magnitude", roundOffPrecision, "\n")
    }

    # rounding radius can't be wider than smallest observed diff.
    # plogis to mitigate the effect of sample size: the larger the sample the more we can «trust» the observed minimal diff
    # obs[1L] = min(obs) = diff of minimal obs with 0
    rr <- .5 * min(stats::plogis(q = length(obs), scale = 11) * diffobs[which(diffobs > 0L)],
                   # rounding precision here
                   roundOffPrecision, obs[1L], na.rm = TRUE)

    ## modify tied observations per group of ties
    # sort() ensures that data after tie-break is still sorted from small to large
    startInd <- endInd <- 1L
    repeat {
      #proceed to end of tie-group
      while (endInd < length(tiesDiffInd) && tiesDiffInd[endInd+1L] == tiesDiffInd[endInd] + 1L) {endInd <- endInd+1L}
      #include adjacent index to complete tie-group
      obsInd <- c(tiesDiffInd[startInd:endInd], tiesDiffInd[endInd]+1L)
      stopifnot( sd(obs[obsInd]) == 0L ) #check: tie-group
      obs[obsInd] <- obs[obsInd] + if (ties == 'random') sort(runif(n = length(obsInd), min = -rr, max = +rr)) else
        #QQQ Cheng (1989) on Moran test statistic proposes to have equidist on CDF-transformed values.
        #+They first use the ties = 'dens' approach for estimation of parameters for Moran's statistic
        seq.int(from = -rr, to = +rr, length.out = length(obsInd)) #evenly spread
      startInd <- endInd <- endInd+1L
      if ( startInd > length(tiesDiffInd) ) break
    }

    if (verbose > 1L){ cat("New data: ", paste(obs, collapse = ", "), "\n")}
  } #fi tiesdiff

  # we have broken all ties
  stopifnot( !any(diff(obs)==0L) )

  obs
} #fn preprocess

# fix ties
h590c <- preprocess(obs = h590, verbose = 1, ties = 'equi')
stopifnot( length(h590) == length(h590c), cor(h590, h590c) > .999 )


cumDiffs.h590c <- diff(c(0, pnorm(q = h590c, mean = 34.072, sd = sqrt(6.874)), 1))
# no ties
stopifnot( length(which(dplyr::near(cumDiffs.h590c, 0))) == 0L )
stopifnot( length(which(cumDiffs.h590c < .Machine$double.eps)) == 0L )

mseCrit.h590c <- -mean(log(cumDiffs.h590c))


EUL_MAS <- -digamma(1L)
mo_m.h590 <- (n.h590+1L) * (log(n.h590+1) + EUL_MAS) - .5 - (12L*(n.h590 + 1L))**-1L
mo_v.h590 <- (n.h590+1L) * (pi**2L / 6L - 1L) - .5 - (6L*(n.h590 + 1L))**-1L

C1.h590 <- mo_m.h590 - sqrt(.5 * n.h590 * mo_v.h590)
C2.h590 <- sqrt(mo_v.h590 / (2L*n.h590))

T_moran.h590 <- (mseCrit.h590c * (n.h590+1L) + .5 * k.h590 - C1.h590) / C2.h590
pchisq(q = T_moran.h590, df = n.h590, lower.tail = FALSE)

# upper point of chisq (cf Cheng, p3)
qchisq(p = .05, df = n.h590, lower.tail = FALSE)
