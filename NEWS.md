# incubate 0.6.1
* make handling of ties in calculation of spacings within the objective function more robust

# incubate 0.6
* implement AD-GOF test within `test_GOF`

# incubate 0.5
* add separate function `test_GOF` for goodness of fit (GOF) tests based on a fitted model
* add Moran's GOF-test

# incubate 0.4
* bug fix in MSE-criterion
* bug fix in tie-handling
* function `estimRoundingError` to estimate rounding error (for tie handling)

# incubate 0.3
* simplified power function (dropped stuff and conventions from `sscn`-package)

# incubate 0.2
* `ties='equidist'` is default method to handle ties
  
# incubate 0.1
* more choice in handling of ties: `equidist` and `random`
* rename function `test_delay_diff` to `test_diff` 

# incubate 0.0.5
* add confidence intervals (based on basic bootstrap) for model parameters

# incubate 0.0.4
* Weibull: `bind=` is implemented for all possible parameter combinations
* internal clean up

# incubate 0.0.3
* more warnings when MSE-optimization fails
* bug fix for power simulation
* bug fix in testing

# incubate 0.0.2
* initial start as a move out from package `sscn` (git history is not moved along, though)
* 'Try' to make code more robust: add `try` code at lower level, directly at where `optim` is called.

