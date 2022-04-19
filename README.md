
<!-- README.md is generated from README.Rmd. Please edit that file -->

# incubate

<!-- badges: start -->
<!-- badges: end -->

Parametric time-to-event analysis where groups show an incubation period
with different hazard.

## Example

With `incubate`, you can statistically compare the survival experience
of two groups with respect to delay and other parameters.

``` r
library(incubate)

# simulate data from exponential distribution with delay
x <- rexp_delayed(n = 13, delay = 1.1)
y <- rexp_delayed(n = 11, delay = 1.6)

# fit delay model 
fm <- delay_model(x, y)

plot(fm)
```

<img src="man/figures/README-example-1.png" width="60%" style="display: block; margin: auto;" />

``` r
# confidence interval for delay-parameters
confint(fm, parm = c('delay.x', 'delay.y'))
#>            2.5%  97.5%
#> delay.x 0.98069 1.2925
#> delay.y 1.45697 2.1287

# test on difference in delay, using only R=75 bootstrap draws
delay_test <- test_diff(x, y, R = 75)
plot(delay_test)
```

<img src="man/figures/README-example-2.png" width="60%" style="display: block; margin: auto;" />

## Parallel computation

To switch on parallel computation, e.g. for bootstrap parameter tests or
power simulations, simply set up a suitable computation plan via the
Future-API. For instance, do the following to activate four R-sessions
in the background of your local computer for computer-intensive tasks in
`incubate`:

``` r
library(future)
plan(multisession, workers = 4)
```

That’s it. You do *not* have to change any function calls. `incubate` is
`future`-aware. Consult the [`future`-package on
CRAN](https://CRAN.R-project.org/package=future) for more information
about futures and about supported computation plans.

When you are done with the heavy computing, it is best practice to
release the parallel connections via `plan(sequential)`.

## Installation

The `incubate` package is hosted publicly at
[Gitlab](https://gitlab.com/imb-dev/incubate). To install its latest
**released version** use from within an R-session:

``` r
remotes::install_gitlab("imb-dev/incubate")
```

To install a specific version, add the version tag after the name,
separated by a `@`, e.g. to install `incubate` in version `v1.0` use

``` r
remotes::install_gitlab("imb-dev/incubate@v1.0")
```

The suffix `@develop` points to the latest **development version** from
Gitlab.
