# Control parameters for `nbreg()`

Tolerances and iteration limits for the three loops used by
[`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md):
the outer loop that alternates the two updates, the inner IRLS loop for
the regression coefficients, and the inner Newton-Raphson loop for
\\\log\theta\\.

## Usage

``` r
nbreg_control(
  tol = 1e-08,
  maxit = 100,
  irls_tol = 1e-10,
  irls_maxit = 100,
  nr_tol = 1e-10,
  nr_maxit = 100,
  theta_max = 1e+08,
  theta_warn = 10000,
  trace = FALSE
)
```

## Arguments

- tol:

  Outer-loop tolerance on the absolute change in the log-likelihood.
  Default `1e-8`.

- maxit:

  Maximum number of outer iterations. Default `100`.

- irls_tol:

  IRLS tolerance on the relative change in the log-likelihood. Default
  `1e-10`.

- irls_maxit:

  Maximum number of IRLS iterations per outer iteration. Default `100`.

- nr_tol:

  Newton-Raphson tolerance on the absolute step in \\\log\theta\\.
  Default `1e-10`.

- nr_maxit:

  Maximum number of Newton-Raphson iterations per outer iteration.
  Default `100`.

- theta_max:

  Upper bound for \\\theta\\. Data with no overdispersion push
  \\\theta\\ towards infinity (the Poisson limit); the bound keeps the
  computations numerically stable. Default `1e8`.

- theta_warn:

  A warning is issued when the final \\\hat\theta\\ exceeds this value,
  because the model is then essentially Poisson. Default `1e4`.

- trace:

  Logical. If `TRUE`, print the log-likelihood and \\\theta\\ after each
  outer iteration.

## Value

A list of control parameters.

## Details

The outer loop stops when the absolute change in the log-likelihood is
below `tol`. The IRLS loop uses the relative change \\\|\ell\_{new} -
\ell\_{old}\| / (\|\ell\_{new}\| + 0.1)\\ \< `irls_tol`. The
Newton-Raphson loop stops when the step in \\\log\theta\\ is smaller
than `nr_tol`.

## See also

[`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md)

## Examples

``` r
nbreg_control(tol = 1e-10, trace = TRUE)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 100
#> 
#> $irls_tol
#> [1] 1e-10
#> 
#> $irls_maxit
#> [1] 100
#> 
#> $nr_tol
#> [1] 1e-10
#> 
#> $nr_maxit
#> [1] 100
#> 
#> $theta_max
#> [1] 1e+08
#> 
#> $theta_warn
#> [1] 10000
#> 
#> $trace
#> [1] TRUE
#> 
```
