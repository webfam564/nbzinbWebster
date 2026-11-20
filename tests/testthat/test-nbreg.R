# Tests for nbreg().  Reference comparisons use MASS::glm.nb() with tight
# tolerances; required agreement (project Table 4, NB column):
#   coefficients 1e-5, standard errors 1e-5, theta 1e-4, log-likelihood 1e-4.

ref_control <- function() stats::glm.control(epsilon = 1e-10, maxit = 100)

max_abs_diff <- function(a, b) max(abs(unname(a) - unname(b)))

load_bio <- function() {
  e <- new.env()
  utils::data("bioChemists", package = "pscl", envir = e)
  e$bioChemists
}

# ---------------------------------------------------------------------------
# Agreement with MASS::glm.nb() on pscl::bioChemists
# ---------------------------------------------------------------------------

test_that("nbreg matches glm.nb on bioChemists (full model)", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- load_bio()
  f <- art ~ fem + mar + kid5 + phd + ment

  ref <- MASS::glm.nb(f, data = d, control = ref_control())
  fit <- nbreg(f, data = d)

  expect_true(fit$converged)
  expect_identical(names(fit$coefficients), names(coef(ref)))

  # Coefficients
  expect_lt(max_abs_diff(fit$coefficients, coef(ref)), 1e-5)
  # Standard errors of beta
  se_ref <- summary(ref)$coefficients[, "Std. Error"]
  expect_lt(max_abs_diff(sqrt(diag(fit$vcov)), se_ref), 1e-5)
  # theta and its SE
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$SE.theta - ref$SE.theta), 1e-4)
  # Log-likelihood
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
  # Deviance and fitted values
  expect_lt(abs(fit$deviance - ref$deviance), 1e-4)
  expect_lt(max_abs_diff(fit$fitted.values, fitted(ref)), 1e-5)
})

test_that("nbreg matches glm.nb with a single covariate", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- load_bio()

  ref <- MASS::glm.nb(art ~ ment, data = d, control = ref_control())
  fit <- nbreg(art ~ ment, data = d)

  expect_lt(max_abs_diff(fit$coefficients, coef(ref)), 1e-5)
  expect_lt(max_abs_diff(sqrt(diag(fit$vcov)),
                         summary(ref)$coefficients[, 2]), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
})

test_that("intercept-only model gives log(mean(y)) and matches glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- load_bio()

  ref <- MASS::glm.nb(art ~ 1, data = d, control = ref_control())
  fit <- nbreg(art ~ 1, data = d)

  # With only an intercept the NB score equation for beta gives
  # sum(y - mu) = 0, so exp(intercept) = mean(y).
  expect_equal(unname(fit$coefficients), log(mean(d$art)), tolerance = 1e-8)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
})

test_that("offsets work, as an argument or in the formula", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- load_bio()
  d$expo <- log(d$phd)  # phd > 0 for all rows

  ref <- MASS::glm.nb(art ~ fem + kid5 + ment + offset(expo), data = d,
                      control = ref_control())
  fit_formula <- nbreg(art ~ fem + kid5 + ment + offset(expo), data = d)
  fit_arg <- nbreg(art ~ fem + kid5 + ment, data = d, offset = expo)

  expect_lt(max_abs_diff(fit_formula$coefficients, coef(ref)), 1e-5)
  expect_lt(abs(fit_formula$theta - ref$theta), 1e-4)
  expect_lt(abs(fit_formula$loglik - as.numeric(logLik(ref))), 1e-4)
  expect_equal(fit_arg$coefficients, fit_formula$coefficients,
               tolerance = 1e-10)
  expect_equal(fit_arg$offset, d$expo)
})

# ---------------------------------------------------------------------------
# Internal consistency
# ---------------------------------------------------------------------------

test_that("log-likelihood equals the sum of dnbinom log-densities", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = d)
  ll <- sum(stats::dnbinom(d$art, size = fit$theta, mu = fit$fitted.values,
                           log = TRUE))
  expect_equal(fit$loglik, ll, tolerance = 1e-10)
  expect_identical(fit$df, 7L)
  expect_identical(fit$nobs, nrow(d))
})

test_that("score equations are zero at the MLE", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = d)
  X <- stats::model.matrix(fit$terms, fit$model)
  mu <- fit$fitted.values
  y <- d$art
  th <- fit$theta
  # Score for beta: X' [(y - mu) theta / (theta + mu)].  The raw score
  # depends on the scale of the covariates, so check the Newton step it
  # implies, (X'WX)^{-1} U, which is on the scale of the coefficients.
  score_beta <- crossprod(X, (y - mu) * th / (th + mu))
  expect_lt(max(abs(fit$vcov %*% score_beta)), 1e-6)
  # Score for theta: the implied Newton step U / (-H) on the theta scale
  a <- rep(1, length(y))
  U <- nbzinbWebster:::nb_score_theta(y, mu, th, a)
  H <- nbzinbWebster:::nb_hess_theta(y, mu, th, a)
  expect_lt(abs(U / H), 1e-6)
})

test_that("integer prior weights are equivalent to replicating rows", {
  set.seed(42)
  n <- 150
  x <- rnorm(n)
  y <- rnbinom(n, mu = exp(0.3 + 0.6 * x), size = 1.5)
  k <- sample(1:3, n, replace = TRUE)
  dw <- data.frame(y = y, x = x, k = k)
  dr <- dw[rep(seq_len(n), k), ]

  fit_w <- nbreg(y ~ x, data = dw, weights = k)
  fit_r <- nbreg(y ~ x, data = dr)

  expect_equal(fit_w$coefficients, fit_r$coefficients, tolerance = 1e-7)
  expect_equal(fit_w$theta, fit_r$theta, tolerance = 1e-6)
  expect_equal(fit_w$loglik, fit_r$loglik, tolerance = 1e-8)
})

test_that("weighted fit matches glm.nb with weights", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- load_bio()
  set.seed(7)
  d$w <- runif(nrow(d), 0.5, 2)

  ref <- MASS::glm.nb(art ~ fem + kid5 + ment, data = d, weights = w,
                      control = ref_control())
  fit <- nbreg(art ~ fem + kid5 + ment, data = d, weights = w)

  expect_lt(max_abs_diff(fit$coefficients, coef(ref)), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
})

test_that("subset and missing values are handled", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  d_na <- d
  d_na$ment[1:5] <- NA

  fit_na <- nbreg(art ~ ment, data = d_na)
  fit_sub <- nbreg(art ~ ment, data = d, subset = -(1:5))

  expect_identical(fit_na$nobs, nrow(d) - 5L)
  expect_equal(fit_na$coefficients, fit_sub$coefficients, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Edge cases, warnings and errors
# ---------------------------------------------------------------------------

test_that("under-dispersed data trigger the Poisson-limit warning", {
  set.seed(1)
  x <- rnorm(200)
  y <- rbinom(200, 3, 0.5)  # variance < mean: no overdispersion
  expect_warning(fit <- nbreg(y ~ x), "theta is very large")
  expect_gt(fit$theta, 1e4)
})

test_that("hitting maxit gives a non-convergence warning", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  expect_warning(fit <- nbreg(art ~ ment, data = d,
                              control = nbreg_control(maxit = 1)),
                 "did not converge")
  expect_false(fit$converged)
})

test_that("trace prints iteration output", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  expect_output(nbreg(art ~ ment, data = d,
                      control = nbreg_control(trace = TRUE)),
                "logLik")
})

test_that("bad input gives informative errors", {
  d <- data.frame(y = c(0, 1, 2, 3, 5, 1), x = 1:6)

  expect_error(nbreg(y ~ x, data = transform(d, y = y - 1)),
               "non-negative integer")
  expect_error(nbreg(y ~ x, data = transform(d, y = y + 0.5)),
               "non-negative integer")
  expect_error(nbreg(y ~ x + I(2 * x), data = d), "rank deficient")
  expect_error(nbreg(y ~ x, data = transform(d, y = 0)), "All responses")
  expect_error(nbreg(y ~ x, data = d, weights = c(-1, 1, 1, 1, 1, 1)),
               "non-negative")
  expect_error(nbreg(y ~ x, data = d, offset = 1:3))
  expect_error(nbreg(y ~ x, data = d, control = "fast"), "must be a list")
  expect_error(nbreg_control(tol = -1), "tol")
  expect_error(nbreg_control(maxit = 0), "maxit")
  expect_error(nbreg_control(trace = "yes"), "trace")
})

# ---------------------------------------------------------------------------
# Project rule: no prohibited functions in the package code
# ---------------------------------------------------------------------------

test_that("package functions do not call prohibited fitting functions", {
  prohibited <- c("optim", "optimize", "optimise", "nlm", "nlminb",
                  "uniroot", "glm", "glm.fit", "lm", "lm.fit", "lm.wfit",
                  "glm.nb", "theta.ml", "zeroinfl")
  ns <- asNamespace("nbzinbWebster")
  fns <- Filter(is.function, mget(ls(ns, all.names = TRUE), envir = ns))
  used <- unique(unlist(lapply(fns, function(f) all.names(body(f)))))
  expect_length(intersect(used, prohibited), 0)
})
