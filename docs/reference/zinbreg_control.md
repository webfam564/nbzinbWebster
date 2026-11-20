# Control parameters for `zinbreg()`

Tolerances and iteration limits for the EM algorithm used by
[`zinbreg()`](https://webfam564.github.io/nbzinbWebster/reference/zinbreg.md)
and for the IRLS fits inside its M-steps.

## Usage

``` r
zinbreg_control(
  tol = 1e-08,
  maxit = 1000,
  mstep_tol = 1e-10,
  mstep_maxit = 100,
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

  EM tolerance on the absolute change in the observed-data
  log-likelihood. Default `1e-8`.

- maxit:

  Maximum number of EM iterations. Default `1000`.

- mstep_tol:

  Outer tolerance (absolute change in the weighted NB log-likelihood)
  for the count-part M-step. Default `1e-10`.

- mstep_maxit:

  Maximum number of outer iterations in the count-part M-step. Default
  `100`.

- irls_tol, irls_maxit:

  Tolerance (relative change in the log-likelihood) and iteration limit
  for the IRLS loops (both parts).

- nr_tol, nr_maxit:

  Tolerance (step in \\\log\theta\\) and iteration limit for the
  Newton-Raphson update of \\\theta\\.

- theta_max, theta_warn:

  Upper bound for \\\theta\\, and the value above which a warning about
  the Poisson limit is given; see
  [`nbreg_control()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg_control.md).

- trace:

  Logical. If `TRUE`, print the log-likelihood after each EM iteration.

## Value

A list of control parameters.

## Details

The EM algorithm stops when the observed-data log-likelihood changes by
less than `tol` (absolute change). Each M-step is fitted to convergence:
the zero-part logistic regression by IRLS (relative log-likelihood
change below `irls_tol`), and the count part by the NB algorithm of
[`nbreg()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg.md)
with outer tolerance `mstep_tol`.

## See also

[`zinbreg()`](https://webfam564.github.io/nbzinbWebster/reference/zinbreg.md),
[`nbreg_control()`](https://webfam564.github.io/nbzinbWebster/reference/nbreg_control.md)

## Examples

``` r
zinbreg_control(tol = 1e-10, maxit = 2000)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 2000
#> 
#> $mstep_tol
#> [1] 1e-10
#> 
#> $mstep_maxit
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
#> [1] FALSE
#> 
```
