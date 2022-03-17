# incubate 0.9
* confint: separate bootstrap data generation and inference step. This allows for more efficient confint-simulations.

# incubate 0.8
* enhance test-coverage for package
* reorganize tests to better match the R-script

# incubate 0.7.6
* `test_diff`: do GOF-tests also for unrestricted model (e.g. `gof_pearson1`) besides for the restricted null model (renamed to `gof_pearson0`)

# incubate 0.7.5
* remove dependency on `dplyr`

# incubate 0.7.4
* add `bs_infer = 't0'` 
* add +3 to degrees of freedom for t-quantiles in `bs_infer='t'` and `='t0'`
* version bump belated: these features were temporarily also released as 0.7.3

# incubate 0.7.3
* fixed mapping of bootstrap inference names for `boot`-package

# incubate 0.7.2
* allow to use `boot` package to calculate confidence interval

# incubate 0.7.1
* allow to choose which tests to perform in `test_diff` (helpful when calculation of AD GOF-tests fails on older R-installations)

# incubate 0.7
* change of S3-classes:
    * for delay model fitting the new class is `incubate_fit`
    * for tests the new class is `incubate_test`
* confidence intervals supports now simple data generation to draw from the data with replacement (besides parametric bootstrap)
* preliminary support for MLE-fitting: this allows to compare confidence intervals based on MSE vs MLE

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

