# Methods for `nbreg` objects

Standard extractor, summary, and prediction methods for negative
binomial fits returned by
[`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("response", "pearson", "deviance"), ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("link", "response"), se.fit = FALSE, ...)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently ignored).

- object, x:

  An object of class `"nbreg"` (or `"summary.nbreg"` for the print
  method of a summary).

- signif.stars:

  Logical; show significance stars?

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"`, `"pearson"`, or `"deviance"`. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"link"` or
  `"response"`.

- newdata:

  Optional data frame in which to look for the covariates (and any
  offset variables). If omitted, the fitted values are used.

- se.fit:

  Logical; if `TRUE`,
  [`predict()`](https://rdrr.io/r/stats/predict.html) also returns
  standard errors (delta method on the response scale).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return numeric
vectors; [`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`"logLik"`; [`nobs()`](https://rdrr.io/r/stats/nobs.html) an integer;
[`predict()`](https://rdrr.io/r/stats/predict.html) a vector, or with
`se.fit = TRUE` a list with components `fit` and `se.fit`;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`"summary.nbreg"`. The print methods return their argument invisibly.

## Details

- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): estimated \\\beta\\ and
  its covariance matrix \\(X^\top\hat W X)^{-1}\\ (with \\\theta\\
  treated as known).

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): the maximized
  log-likelihood with attribute `df` = \\p + 1\\ (the coefficients plus
  \\\theta\\), so [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work.

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): fitted means
  \\\hat\mu_i\\.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html): response
  residuals \\y_i - \hat\mu_i\\, Pearson residuals \\\sqrt{a_i}(y_i -
  \hat\mu_i)/\sqrt{V(\hat\mu_i)}\\ with \\V(\mu) = \mu + \mu^2/\theta\\,
  or deviance residuals.

- [`predict()`](https://rdrr.io/r/stats/predict.html): the linear
  predictor \\\hat\eta\\ (`type = "link"`) or the mean \\\hat\mu =
  e^{\hat\eta}\\ (`type = "response"`), for the original data or for
  `newdata`, optionally with standard errors.

- [`summary()`](https://rdrr.io/r/base/summary.html): coefficient table
  with Wald z tests, plus \\\hat\theta\\, its SE, the log-likelihood,
  AIC, and deviance.

## Examples

``` r
set.seed(1)
n <- 300
x <- rnorm(n)
y <- rnbinom(n, mu = exp(0.5 + 0.8 * x), size = 2)
fit <- nbreg(y ~ x)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x)
#> 
#> Coefficients (negative binomial with log link):
#> (Intercept)            x  
#>      0.4508       0.7669  
#> 
#> Theta: 1.95   Log-likelihood: -526.892 on 3 df
#> 
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1939 -0.7791 -0.2593  0.4389  5.2098 
#> 
#> Coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.45081    0.06613   6.817 9.32e-12 ***
#> x            0.76687    0.06699  11.448  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 1.95  (SE 0.3699)    Log(theta): 0.6678  (SE 0.1897) 
#> Log-likelihood: -526.892 on 3 df   AIC: 1059.78 
#> Residual deviance: 325.54 on 298 degrees of freedom
#> Outer iterations (IRLS + Newton-Raphson): 3 (converged) 
#> 
coef(fit)
#> (Intercept)           x 
#>   0.4508058   0.7668655 
vcov(fit)
#>              (Intercept)            x
#> (Intercept)  0.004373642 -0.001688731
#> x           -0.001688731  0.004487089
logLik(fit)
#> 'log Lik.' -526.8916 (df=3)
AIC(fit)
#> [1] 1059.783
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -0.8050751  2.2472901 -0.7620322  1.2691926 -0.5032198  0.1494196 
predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "response")
#>        1        2        3 
#> 0.729016 1.569576 3.379309 
```
