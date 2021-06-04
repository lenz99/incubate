
<!-- README.md is generated from README.Rmd. Please edit that file -->

# incubate

<!-- badges: start -->
<!-- badges: end -->

Parametric time-to-event analysis where groups show an incubation period
where hazard is different.

## Installation

You can install the released version of incubate from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("incubate")
```

## Example

With incubate, you can statistically compare the survival experience of
two groups with respect to delay and other parameters.

``` r
library(incubate)

# simulate data from exponential distribution with delay
x <- rexp_delayed(n = 13, delay = 1.1)
y <- rexp_delayed(n = 11, delay = 1.6)

# fit delay model 
fm <- delay_model(x, y)

plot(fm)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
delay_test <- test_delay_diff(x, y, R = 100)
plot(delay_test)
```

<img src="man/figures/README-example-2.png" width="100%" />

## Installation

You can install the `sscn` package from
[Gitlab](https://gitlab.com/imb-dev/incubate) with:

``` r
remotes::install_gitlab("imb-dev/incubate")
```

To install a specific version, add the version tag after the name,
separated by a `@`, e.g. to install `incubate` in version `v0.0.1` use

``` r
remotes::install_gitlab("imb-dev/incubate@v0.0.1")
```
