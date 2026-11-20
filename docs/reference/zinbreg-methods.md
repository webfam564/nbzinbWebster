# Methods for `zinbreg` objects

Standard extractor, summary, and prediction methods for zero-inflated
negative binomial fits returned by
[`zinbreg()`](https://webfam564.github.io/nbzinbWebster/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'zinbreg'
predict(
  object,
  newdata,
  type = c("response", "count", "zero", "prob"),
  at = NULL,
  ...
)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently ignored).

- object, x:

  An object of class `"zinbreg"` (or `"summary.zinbreg"` for the print
  method of a summary).

- signif.stars:

  Logical; show significance stars?

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"`, `"count"`, or
  `"zero"`.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` or `"pearson"`. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`,
  `"count"`, `"zero"`, or `"prob"`.

- newdata:

  Optional data frame in which to look for the covariates (and any
  count-offset variables). If omitted, the fitted values are used.

- at:

  For `predict(type = "prob")`: the counts at which to evaluate the
  probabilities. Default `0:max(y)`.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return numeric
vectors; [`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`"logLik"`; [`nobs()`](https://rdrr.io/r/stats/nobs.html) an integer;
[`predict()`](https://rdrr.io/r/stats/predict.html) a vector, or a
matrix for `type = "prob"`;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`"summary.zinbreg"`. The print methods return their argument invisibly.

## Details

- [`coef()`](https://rdrr.io/r/stats/coef.html): \\\hat\beta\\ and
  \\\hat\gamma\\, either together (`model = "full"`, named `count_*` and
  `zero_*` as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)) or
  one part at a time.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html): the matching block of
  the inverse observed information (without \\\log\theta\\, whose SE is
  reported separately).

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): the maximized
  observed-data log-likelihood with attribute `df` = \\p + q + 1\\, so
  [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work.

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): the fitted
  mean \\(1-\hat\pi_i)\hat\mu_i\\.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html): response
  residuals \\y_i - (1-\hat\pi_i)\hat\mu_i\\, or Pearson residuals that
  divide them by \\\sqrt{\mathrm{Var}(Y_i)}\\, \\\mathrm{Var}(Y_i) =
  (1-\pi_i)\mu_i(1 + \pi_i\mu_i + \mu_i/\theta)\\.

- [`predict()`](https://rdrr.io/r/stats/predict.html):
  `type = "response"` gives \\(1-\hat\pi)\hat\mu\\, `"count"` gives
  \\\hat\mu\\, `"zero"` gives \\\hat\pi\\, and `"prob"` gives the matrix
  of \\P(Y_i = k)\\ for `k` in `at`.

- [`summary()`](https://rdrr.io/r/base/summary.html): coefficient tables
  for both parts with Wald z tests; the count table also shows
  \\\log\theta\\.

## Examples

``` r
set.seed(2)
n <- 500
x <- rnorm(n)
z <- rnorm(n)
y <- ifelse(rbinom(n, 1, plogis(-1 + z)) == 1, 0,
            rnbinom(n, mu = exp(1 + 0.5 * x), size = 2))
fit <- zinbreg(y ~ x | z)
coef(fit)
#> count_(Intercept)           count_x  zero_(Intercept)            zero_z 
#>         0.9657780         0.4698595        -0.9633195         1.0990002 
coef(fit, model = "zero")
#> (Intercept)           z 
#>  -0.9633195   1.0990002 
sqrt(diag(vcov(fit)))
#> count_(Intercept)           count_x  zero_(Intercept)            zero_z 
#>        0.06541699        0.05525722        0.22574745        0.20567462 
logLik(fit)
#> 'log Lik.' -871.4509 (df=5)
AIC(fit)
#> [1] 1752.902
nd <- data.frame(x = c(-1, 0, 1), z = c(0, 0, 2))
predict(fit, newdata = nd, type = "response")
#>         1         2         3 
#> 1.1884607 1.9012631 0.9470491 
predict(fit, newdata = nd, type = "zero")
#>         1         2         3 
#> 0.2762141 0.2762141 0.7746368 
# expected number of zeros
sum(predict(fit, type = "prob", at = 0))
#> [1] 224.6451
```
