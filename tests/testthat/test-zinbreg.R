# Tests for zinbreg() and its methods.  Reference comparisons use
# pscl::zeroinfl(dist = "negbin") with a tight tolerance; required agreement
# (project Table 4, ZINB column): coefficients, SEs, log(theta) and
# log-likelihood all within 1e-3.

tol_zinb <- 1e-3
max_abs_diff <- function(a, b) max(abs(unname(a) - unname(b)))

load_bio <- function() {
  e <- new.env()
  utils::data("bioChemists", package = "pscl", envir = e)
  e$bioChemists
}

ref_fit <- function(formula, data, ...) {
  pscl::zeroinfl(formula, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12), ...)
}

# Fit the full bioChemists model once and reuse it in several tests.
full_cache <- new.env()
full_pair <- function() {
  if (is.null(full_cache$fit)) {
    d <- load_bio()
    f <- art ~ fem + mar + kid5 + phd + ment
    full_cache$d <- d
    full_cache$fit <- zinbreg(f, data = d)
    full_cache$ref <- ref_fit(f, d)
  }
  as.list(full_cache)
}

# ---------------------------------------------------------------------------
# Agreement with pscl::zeroinfl() on pscl::bioChemists
# ---------------------------------------------------------------------------

test_that("zinbreg matches zeroinfl on bioChemists (full model)", {
  skip_if_not_installed("pscl")
  p <- full_pair()

  expect_true(p$fit$converged)
  expect_identical(names(coef(p$fit)), names(coef(p$ref)))
  # Coefficients (count and zero parts)
  expect_lt(max_abs_diff(coef(p$fit), coef(p$ref)), tol_zinb)
  # Standard errors of the coefficients
  expect_lt(max_abs_diff(sqrt(diag(vcov(p$fit))), sqrt(diag(vcov(p$ref)))),
            tol_zinb)
  # log(theta) and its SE
  expect_lt(abs(p$fit$log_theta - log(p$ref$theta)), tol_zinb)
  expect_lt(abs(p$fit$SE.log_theta - p$ref$SE.logtheta), tol_zinb)
  # Log-likelihood
  expect_lt(abs(p$fit$loglik - as.numeric(logLik(p$ref))), tol_zinb)
})

test_that("zinbreg matches zeroinfl with different covariates per part", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  f <- art ~ fem + mar + kid5 + ment | ment
  fit <- zinbreg(f, data = d)
  ref <- ref_fit(f, d)

  expect_identical(names(coef(fit)), names(coef(ref)))
  expect_lt(max_abs_diff(coef(fit), coef(ref)), tol_zinb)
  expect_lt(max_abs_diff(sqrt(diag(vcov(fit))), sqrt(diag(vcov(ref)))),
            tol_zinb)
  expect_lt(abs(fit$log_theta - log(ref$theta)), tol_zinb)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), tol_zinb)
})

test_that("count offsets match zeroinfl, as argument or in the formula", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  d$expo <- log(d$phd)  # phd > 0 for all rows

  ref <- ref_fit(art ~ fem + kid5 + ment + offset(expo) | ment, d)
  fit_f <- zinbreg(art ~ fem + kid5 + ment + offset(expo) | ment, data = d)
  fit_a <- zinbreg(art ~ fem + kid5 + ment | ment, data = d, offset = expo)

  expect_lt(max_abs_diff(coef(fit_f), coef(ref)), tol_zinb)
  expect_lt(abs(fit_f$loglik - as.numeric(logLik(ref))), tol_zinb)
  expect_equal(coef(fit_a), coef(fit_f), tolerance = 1e-8)
  expect_equal(fit_a$offset, d$expo)
})

test_that("a formula without | uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  fit1 <- zinbreg(art ~ kid5 + ment, data = d)
  fit2 <- zinbreg(art ~ kid5 + ment | kid5 + ment, data = d)
  expect_identical(names(coef(fit1, "zero")), c("(Intercept)", "kid5", "ment"))
  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# EM behavior and standard errors
# ---------------------------------------------------------------------------

test_that("EM never decreases the log-likelihood", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  path <- p$fit$loglik_path
  expect_identical(length(path), p$fit$iter + 1L)
  expect_true(all(diff(path) > -1e-8))
  expect_equal(p$fit$loglik, path[length(path)])
})

test_that("log-likelihood agrees with a direct computation from dnbinom", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  y <- p$d$art
  pi0 <- p$fit$pi
  f <- stats::dnbinom(y, size = p$fit$theta, mu = p$fit$mu)
  ll <- sum(log(ifelse(y == 0, pi0 + (1 - pi0) * f, (1 - pi0) * f)))
  expect_equal(p$fit$loglik, ll, tolerance = 1e-10)
})

test_that("E-step gives tau = 0 for positive counts and Bayes' rule at 0", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  y <- p$d$art
  expect_true(all(p$fit$tau[y > 0] == 0))
  pi0 <- p$fit$pi[y == 0]
  f0 <- stats::dnbinom(0, size = p$fit$theta, mu = p$fit$mu[y == 0])
  expect_equal(unname(p$fit$tau[y == 0]), unname(pi0 / (pi0 + (1 - pi0) * f0)),
               tolerance = 1e-10)
})

test_that("analytic Hessian matches a numerical Hessian", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("numDeriv")
  p <- full_pair()
  fit <- p$fit
  X <- stats::model.matrix(fit$terms$count, fit$model)
  Z <- stats::model.matrix(fit$terms$zero, fit$model)
  y <- fit$y
  a <- fit$prior.weights
  np <- ncol(X)
  nq <- ncol(Z)
  ll <- function(psi) {
    nbzinbWebster:::zinb_loglik(
      y, exp(as.vector(X %*% psi[1:np])), exp(psi[np + nq + 1]),
      as.vector(Z %*% psi[np + 1:nq]), a)
  }
  psi <- c(coef(fit), fit$log_theta)
  H_num <- numDeriv::hessian(ll, psi)
  H_ana <- nbzinbWebster:::zinb_hessian(X, Z, y, fit$mu, fit$theta,
                                        stats::qlogis(fit$pi), a)
  expect_lt(max(abs(H_ana - H_num)) / max(abs(H_num)), 1e-6)
  # SE(theta) by the delta method
  expect_equal(fit$SE.theta, fit$theta * fit$SE.log_theta)
})

# ---------------------------------------------------------------------------
# S3 methods, compared with zeroinfl where possible
# ---------------------------------------------------------------------------

test_that("coef, vcov, logLik, AIC, BIC, nobs agree with zeroinfl", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  expect_lt(max_abs_diff(coef(p$fit, "count"), coef(p$ref, "count")),
            tol_zinb)
  expect_lt(max_abs_diff(coef(p$fit, "zero"), coef(p$ref, "zero")),
            tol_zinb)
  expect_identical(names(coef(p$fit, "count")), names(coef(p$ref, "count")))
  expect_identical(dim(vcov(p$fit, "zero")), c(6L, 6L))
  expect_identical(dimnames(vcov(p$fit)), dimnames(vcov(p$ref)))

  ll <- logLik(p$fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"), attr(logLik(p$ref), "df"))
  expect_lt(abs(AIC(p$fit) - AIC(p$ref)), 1e-2)
  expect_lt(abs(BIC(p$fit) - BIC(p$ref)), 1e-2)
  # pscl has no nobs() method for zeroinfl; it stores n instead.
  expect_identical(nobs(p$fit), p$ref$n)
})

test_that("fitted and residuals agree with zeroinfl", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  expect_lt(max_abs_diff(fitted(p$fit), fitted(p$ref)), tol_zinb)
  expect_lt(max_abs_diff(residuals(p$fit), residuals(p$ref, "response")),
            tol_zinb)
  expect_lt(max_abs_diff(residuals(p$fit, "pearson"),
                         residuals(p$ref, "pearson")), tol_zinb)
  expect_error(residuals(p$fit, type = "deviance"))
})

test_that("predict agrees with zeroinfl for every type, with newdata", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  nd <- p$d[c(1, 100, 500, 915), ]

  for (tp in c("response", "count", "zero")) {
    expect_lt(max_abs_diff(predict(p$fit, type = tp),
                           predict(p$ref, type = tp)), tol_zinb)
    expect_lt(max_abs_diff(predict(p$fit, newdata = nd, type = tp),
                           predict(p$ref, newdata = nd, type = tp)), tol_zinb)
  }
  expect_equal(predict(p$fit, type = "response"), fitted(p$fit))

  pr <- predict(p$fit, newdata = nd, type = "prob", at = 0:5)
  expect_identical(dim(pr), c(4L, 6L))
  expect_lt(max_abs_diff(pr, predict(p$ref, newdata = nd, type = "prob",
                                     at = 0:5)), tol_zinb)
  # Probabilities over 0..max(y) for the data sum to (almost) 1.
  pr_all <- predict(p$fit, type = "prob", at = 0:200)
  expect_true(all(abs(rowSums(pr_all) - 1) < 1e-6))
})

test_that("predict uses the count offset with newdata", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  d$expo <- log(d$phd)
  nd <- d[1:5, ]
  fit_f <- zinbreg(art ~ ment + offset(expo) | ment, data = d)
  fit_a <- zinbreg(art ~ ment | ment, data = d, offset = expo)
  manual <- exp(as.vector(cbind(1, nd$ment) %*% coef(fit_f, "count") +
                            nd$expo))
  expect_equal(unname(predict(fit_f, newdata = nd, type = "count")), manual,
               tolerance = 1e-10)
  expect_equal(unname(predict(fit_a, newdata = nd, type = "count")), manual,
               tolerance = 1e-6)
})

test_that("summary tables agree with zeroinfl and print works", {
  skip_if_not_installed("pscl")
  p <- full_pair()
  s <- summary(p$fit)
  s_ref <- summary(p$ref)$coefficients
  expect_s3_class(s, "summary.zinbreg")
  expect_identical(rownames(s$coefficients$count), rownames(s_ref$count))
  expect_lt(max_abs_diff(s$coefficients$count[, 1:2], s_ref$count[, 1:2]),
            tol_zinb)
  expect_lt(max_abs_diff(s$coefficients$zero[, 1:2], s_ref$zero[, 1:2]),
            tol_zinb)

  expect_output(out <- print(p$fit), "Zero-inflation model")
  expect_identical(out, p$fit)
  expect_output(out <- print(s), "Log\\(theta\\)")
  expect_identical(out, s)
})

# ---------------------------------------------------------------------------
# Edge cases, warnings and errors
# ---------------------------------------------------------------------------

test_that("hitting maxit gives a non-convergence warning", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  expect_warning(fit <- zinbreg(art ~ ment | ment, data = d,
                                control = zinbreg_control(maxit = 2)),
                 "did not converge")
  expect_false(fit$converged)
})

test_that("trace prints EM iterations", {
  skip_if_not_installed("pscl")
  d <- load_bio()
  expect_output(zinbreg(art ~ ment | 1, data = d,
                        control = zinbreg_control(trace = TRUE)),
                "EM iter")
})

test_that("bad input gives informative errors", {
  d <- data.frame(y = c(0, 1, 2, 3, 5, 1, 0, 4), x = 1:8, z = c(1, 0))

  expect_error(zinbreg(y ~ x | z, data = transform(d, y = y + 1)),
               "no zeros")
  expect_error(zinbreg(y ~ x | z, data = transform(d, y = y - 1)),
               "non-negative integer")
  expect_error(zinbreg(y ~ x | z + I(2 * z), data = d), "rank deficient")
  expect_error(zinbreg(y ~ x | z + offset(x), data = d),
               "zero-inflation part")
  expect_error(zinbreg(~ x, data = d), "two-sided")
  expect_error(zinbreg(y ~ x | z, data = d, control = 1), "must be a list")
  expect_error(zinbreg_control(tol = 0), "tol")
  expect_error(zinbreg_control(trace = NA), "trace")
})
