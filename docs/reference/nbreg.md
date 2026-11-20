# Negative binomial regression

Fits a negative binomial (NB2) regression model with log link,
\\\log\mu_i = x_i^\top\beta + o_i\\ and \\\mathrm{Var}(Y_i) = \mu_i +
\mu_i^2/\theta\\. The coefficients are estimated by iteratively
reweighted least squares (IRLS) with \\\theta\\ held fixed, alternated
with Newton-Raphson updates of \\\log\theta\\ with \\\beta\\ held fixed,
until the log-likelihood stops changing. No general-purpose optimizer is
used.

## Usage

``` r
nbreg(
  formula,
  data,
  offset,
  control = nbreg_control(),
  weights,
  subset,
  na.action
)
```

## Arguments

- formula:

  A model formula, `response ~ terms`. Offsets may be given with
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the
  formula.

- data:

  An optional data frame containing the variables in the model.

- offset:

  An optional vector of offsets (e.g. log exposure), added to any
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the
  formula.

- control:

  A list of control parameters; see
  [`nbreg_control()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg_control.md).

- weights:

  Optional non-negative prior weights \\a_i\\.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values (default: the `na.action` option, usually
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html)).

## Value

An object of class `"nbreg"`, a list with components

- coefficients:

  Named vector of estimated \\\beta\\.

- vcov:

  Estimated covariance matrix of \\\hat\beta\\.

- theta, SE.theta:

  \\\hat\theta\\ and its standard error.

- log_theta, SE.log_theta:

  \\\log\hat\theta\\ and its standard error (delta method).

- loglik:

  Maximized log-likelihood.

- df:

  Number of estimated parameters (\\p + 1\\).

- deviance:

  Residual deviance at \\\hat\theta\\.

- iter:

  Number of outer iterations.

- converged:

  Logical convergence flag.

- fitted.values, linear.predictors, residuals:

  \\\hat\mu_i\\, \\\hat\eta_i\\, and response residuals \\y_i -
  \hat\mu_i\\.

- weights, prior.weights:

  Final IRLS working weights and the prior weights.

- y, offset, nobs, df.residual:

  Response, offset, number of observations with positive weight, and
  \\n - p\\.

- call, formula, terms, model, xlevels, contrasts, na.action, control:

  Information about the fit.

## Algorithm

1.  **Starting values.** \\\beta\\ from a Poisson fit computed by IRLS
    (working weights \\a_i\mu_i\\); \\\theta\\ from the moment estimate
    \\\sum a_i\mu_i^2 / \sum a_i\\(y_i-\mu_i)^2-\mu_i\\\\.

2.  **IRLS for \\\beta\\** (\\\theta\\ fixed): with \\w_i =
    a_i\mu_i\theta/(\theta+\mu_i)\\ and \\z_i = (\eta_i - o_i) +
    (y_i-\mu_i)/\mu_i\\, update \\\beta \leftarrow (X^\top W
    X)^{-1}X^\top W z\\.

3.  **Newton-Raphson for \\\phi=\log\theta\\** (\\\beta\\ fixed): \\\phi
    \leftarrow \phi - U\_\phi/H\_\phi\\ with \\U\_\phi = \theta
    U(\theta)\\ and \\H\_\phi = \theta^2 H(\theta) + \theta U(\theta)\\.

4.  Repeat 2-3 until the log-likelihood changes by less than
    `control$tol`, then do a final IRLS pass at the final \\\theta\\.

Standard errors: \\\widehat{\mathrm{Var}}(\hat\beta) = (X^\top\hat W
X)^{-1}\\ (treating \\\theta\\ as known) and \\\mathrm{SE}(\hat\theta) =
\\-H(\hat\theta)\\^{-1/2}\\.

## References

Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics
with S* (4th ed.). Springer.

## See also

[`nbreg_control()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg_control.md)

## Examples

``` r
set.seed(1)
n <- 300
x <- rnorm(n)
y <- rnbinom(n, mu = exp(0.5 + 0.8 * x), size = 2)
fit <- nbreg(y ~ x)
fit$coefficients
#> (Intercept)           x 
#>   0.4508058   0.7668655 
fit$theta
#> [1] 1.950002
fit$loglik
#> [1] -526.8916

# with an exposure offset
t <- runif(n, 1, 5)
y2 <- rnbinom(n, mu = t * exp(0.5 + 0.8 * x), size = 2)
nbreg(y2 ~ x, offset = log(t))$coefficients
#> (Intercept)           x 
#>   0.5102294   0.6970507 
```
