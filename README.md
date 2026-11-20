# nbzinbWebster

<!-- badges: start -->
<!-- CI badge placeholder: the R-CMD-check badge goes here. -->
<!-- badges: end -->

**nbzinbWebster** fits negative binomial (NB2) and zero-inflated negative
binomial (ZINB) regression models with algorithms written from scratch:

* `nbreg()`: NB regression by iteratively reweighted least squares (IRLS)
  for the coefficients, alternated with Newton–Raphson for log(θ).
* `zinbreg()`: ZINB regression by the EM algorithm, with IRLS M-steps (a
  weighted logistic regression for the zero part and a weighted NB fit for
  the count part). Standard errors come from the analytic Hessian of the
  observed-data log-likelihood.

No general-purpose optimizer (`optim()`, `nlm()`, `glm()`, ...) is used in the
fitting code. The results reproduce `MASS::glm.nb()` and
`pscl::zeroinfl(dist = "negbin")`.

This package is the STA 651 (Fall 2026) final project of Simon Webster.

* Website: <https://webfam564.github.io/nbzinbWebster/>
* Report: [`report/report.pdf`](report/report.pdf)

## Installation

```r
# install.packages("devtools")
devtools::install_github("webfam564/nbzinbWebster")
```

## Usage

The `bioChemists` data (from the **pscl** package) record the number of
articles published by 915 biochemistry PhD students.

```r
library(nbzinbWebster)
data("bioChemists", package = "pscl")

# Negative binomial regression
fit_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
fit_nb
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#>
#> Coefficients (negative binomial with log link):
#> (Intercept)     femWomen   marMarried         kid5          phd         ment
#>     0.25614     -0.21642      0.15049     -0.17642      0.01527      0.02908
#>
#> Theta: 2.264   Log-likelihood: -1560.96 on 7 df

# Zero-inflated negative binomial regression:
# count covariates before the |, zero-inflation covariates after it
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                    data = bioChemists)
summary(fit_zinb)
#> Count model coefficients (negbin with log link):
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.404014   0.141715   2.851  0.00436 **
#> femWomen    -0.211906   0.071922  -2.946  0.00322 **
#> marMarried   0.139462   0.081192   1.718  0.08586 .
#> kid5        -0.167624   0.052457  -3.195  0.00140 **
#> phd          0.001963   0.035586   0.055  0.95602
#> ment         0.024393   0.003518   6.934 4.08e-12 ***
#> Log(theta)   1.002910   0.142830   7.022 2.19e-12 ***
#>
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  -0.8067     0.3532  -2.284   0.0224 *
#> ment         -0.6094     0.2458  -2.480   0.0131 *
#>
#> Theta = 2.726  (SE 0.3894)
#> Log-likelihood: -1553.27 on 9 df   AIC: 3124.54
#> EM iterations: 112 (converged)

# The usual S3 methods work
AIC(fit_nb, fit_zinb)
#>          df      AIC
#> fit_nb    7 3135.917
#> fit_zinb  9 3124.542
head(predict(fit_zinb, type = "zero"))   # P(structural zero)
```

Both fitted objects support `print`, `summary`, `coef`, `vcov`, `logLik`,
`AIC`/`BIC`, `fitted`, `residuals` (response and Pearson), `nobs` and
`predict` (with `newdata`). For `zinbreg`, `predict()` takes
`type = "response"`, `"count"`, `"zero"` or `"prob"`. Tolerances and
iteration limits are set with `nbreg_control()` and `zinbreg_control()`.

## Benchmarks and simulation

`inst/benchmark/benchmark.R` compares the package with `MASS::glm.nb()` and
`pscl::zeroinfl()` on `pscl::bioChemists` and `AER::NMES1988`.
`inst/benchmark/simulation.R` runs a simulation study of `zinbreg()`. The
results are discussed in the report.
