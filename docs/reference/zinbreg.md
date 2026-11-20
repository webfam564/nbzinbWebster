# Zero-inflated negative binomial regression

Fits a zero-inflated negative binomial (ZINB) regression model. With
probability \\\pi_i\\ observation \\i\\ is a structural zero; with
probability \\1-\pi_i\\ it comes from an NB(\\\mu_i, \theta\\)
distribution, where \\\log\mu_i = x_i^\top\beta + o_i\\ and
\\\mathrm{logit}\\\pi_i = z_i^\top\gamma\\. The model is fitted by an EM
algorithm whose M-steps are IRLS fits. No general-purpose optimizer is
used.

## Usage

``` r
zinbreg(
  formula,
  data,
  offset,
  control = zinbreg_control(),
  weights,
  subset,
  na.action
)
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before the `|`, zero-model covariates after it. Without `|`, the same
  covariates are used in both parts. Offsets (count part only) may be
  given with [`offset()`](https://rdrr.io/r/stats/offset.html) terms
  before the `|`.

- data:

  An optional data frame containing the variables in the model.

- offset:

  An optional offset for the count model (e.g. log exposure), added to
  any [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the
  formula.

- control:

  A list of control parameters; see
  [`zinbreg_control()`](https://webfam564.github.io/nbzinbWebster/reference/zinbreg_control.md).

- weights:

  Optional non-negative prior weights \\a_i\\.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values (default: the `na.action` option, usually
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html)).

## Value

An object of class `"zinbreg"`, a list with components

- coefficients:

  List with elements `count` (\\\beta\\) and `zero` (\\\gamma\\).

- vcov:

  Covariance matrix of \\(\hat\beta, \hat\gamma, \log\hat\theta)\\;
  [`vcov.zinbreg()`](https://webfam564.github.io/nbzinbWebster/reference/zinbreg-methods.md)
  returns the block for the coefficients.

- theta, SE.theta, log_theta, SE.log_theta:

  \\\hat\theta\\, \\\log\hat\theta\\, and their standard errors.

- loglik, loglik_path:

  Maximized observed-data log-likelihood, and its value at the start and
  after each EM iteration.

- df:

  Number of estimated parameters (\\p + q + 1\\).

- iter, converged:

  Number of EM iterations and convergence flag.

- fitted.values:

  \\(1-\hat\pi_i)\hat\mu_i\\, the fitted mean.

- mu, pi, tau:

  \\\hat\mu_i\\, \\\hat\pi_i\\, and the final posterior probabilities
  \\\tau_i\\.

- residuals:

  Response residuals \\y_i - (1-\hat\pi_i)\hat\mu_i\\.

- y, offset, prior.weights, nobs, df.residual:

  Response, count offset, prior weights, number of observations with
  positive weight, and \\n - p - q - 1\\.

- call, formula, terms, model, xlevels, contrasts, na.action, control:

  Information about the fit; `terms`, `xlevels` and `contrasts` are
  lists with `count` and `zero` elements.

## Algorithm

Starting values: \\(\beta, \theta)\\ from an NB fit (as in
[`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md));
\\\gamma\\ from a logistic regression of \\I(y_i = 0)\\ on \\z_i\\,
fitted by IRLS. Then repeat:

1.  **E-step.** Posterior probability of a structural zero, \\\tau_i =
    \pi_i / \\\pi_i + (1-\pi_i) f\_{NB}(0;\mu_i,\theta)\\\\ if \\y_i =
    0\\, and \\\tau_i = 0\\ if \\y_i \> 0\\.

2.  **M-step, zero part.** Maximize \\\sum_i a_i\\\tau_i\log\pi_i +
    (1-\tau_i)\log(1-\pi_i)\\\\: a logistic regression with fractional
    response \\\tau_i\\, fitted by IRLS with \\w_i = a_i\pi_i(1-\pi_i)\\
    and working response \\z_i^\top\gamma +
    (\tau_i-\pi_i)/\\\pi_i(1-\pi_i)\\\\.

3.  **M-step, count part.** Maximize \\\sum_i a_i(1-\tau_i)\log
    f\_{NB}(y_i;\mu_i,\theta)\\: the NB fit of
    [`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md)
    with prior weights \\a_i(1-\tau_i)\\.

until the observed-data log-likelihood changes by less than
`control$tol`. EM never decreases the log-likelihood; the path is stored
in `loglik_path`, and a warning is given if it ever decreases.

Standard errors come from the inverse of the observed information,
\\-\partial^2\ell/\partial\psi\partial\psi^\top\\, with \\\psi = (\beta,
\gamma, \log\theta)\\, computed from the analytic second derivatives of
the observed-data log-likelihood at the MLE. \\\mathrm{SE}(\hat\theta) =
\hat\theta\\\mathrm{SE}(\log\hat\theta)\\ (delta method).

## References

Lambert, D. (1992). Zero-inflated Poisson regression, with an
application to defects in manufacturing. *Technometrics*, 34(1), 1-14.

Zeileis, A., Kleiber, C. and Jackman, S. (2008). Regression models for
count data in R. *Journal of Statistical Software*, 27(8), 1-25.

## See also

[`zinbreg_control()`](https://webfam564.github.io/nbzinbWebster/reference/zinbreg_control.md),
[`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md)

## Examples

``` r
set.seed(2)
n <- 500
x <- rnorm(n)
z <- rnorm(n)
structural_zero <- rbinom(n, 1, plogis(-1 + 1 * z))
y <- ifelse(structural_zero == 1, 0,
            rnbinom(n, mu = exp(1 + 0.5 * x), size = 2))
fit <- zinbreg(y ~ x | z)
fit
#> 
#> Call:
#> zinbreg(formula = y ~ x | z)
#> 
#> Count model coefficients (negbin with log link):
#> (Intercept)            x  
#>      0.9658       0.4699  
#> Theta = 2.291 
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#> (Intercept)            z  
#>     -0.9633       1.0990  
#> 
#> Log-likelihood: -871.451 on 5 df
#> 
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | z)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1592 -0.6860 -0.3809  0.4042  4.9631 
#> 
#> Count model coefficients (negbin with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.96578    0.06542  14.763  < 2e-16 ***
#> x            0.46986    0.05526   8.503  < 2e-16 ***
#> Log(theta)   0.82895    0.22871   3.624  0.00029 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.9633     0.2257  -4.267 1.98e-05 ***
#> z             1.0990     0.2057   5.343 9.12e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta = 2.291  (SE 0.524) 
#> Log-likelihood: -871.451 on 5 df   AIC: 1752.9 
#> EM iterations: 40 (converged) 
#> 
```
