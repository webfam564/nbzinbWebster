# Tests for the S3 methods of class "nbreg".  Where MASS::glm.nb has the
# same method, the results are compared with it.

fit_pair <- function() {
  e <- new.env()
  utils::data("bioChemists", package = "pscl", envir = e)
  d <- e$bioChemists
  f <- art ~ fem + mar + kid5 + phd + ment
  list(
    d = d,
    fit = nbreg(f, data = d),
    ref = MASS::glm.nb(f, data = d,
                       control = stats::glm.control(epsilon = 1e-10,
                                                    maxit = 100))
  )
}

test_that("coef, vcov, fitted, nobs agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  p <- fit_pair()
  expect_equal(coef(p$fit), coef(p$ref), tolerance = 1e-5)
  expect_equal(vcov(p$fit), vcov(p$ref), tolerance = 1e-5)
  expect_equal(fitted(p$fit), fitted(p$ref), tolerance = 1e-5)
  expect_identical(nobs(p$fit), nobs(p$ref))
})

test_that("logLik has df, so AIC and BIC agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  p <- fit_pair()
  ll <- logLik(p$fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"), 7L)
  expect_equal(as.numeric(ll), as.numeric(logLik(p$ref)), tolerance = 1e-8)
  expect_equal(AIC(p$fit), AIC(p$ref), tolerance = 1e-8)
  expect_equal(BIC(p$fit), BIC(p$ref), tolerance = 1e-8)
})

test_that("residuals of each type agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  p <- fit_pair()
  expect_equal(residuals(p$fit), residuals(p$ref, type = "response"),
               tolerance = 1e-5)
  expect_equal(residuals(p$fit, type = "pearson"),
               residuals(p$ref, type = "pearson"), tolerance = 1e-5)
  expect_equal(residuals(p$fit, type = "deviance"),
               residuals(p$ref, type = "deviance"), tolerance = 1e-5)
  # Squared deviance residuals sum to the deviance.
  expect_equal(sum(residuals(p$fit, type = "deviance")^2), p$fit$deviance,
               tolerance = 1e-10)
  expect_error(residuals(p$fit, type = "working"))
})

test_that("predict agrees with glm.nb, with and without newdata", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  p <- fit_pair()
  nd <- p$d[c(1, 50, 300, 915), ]

  # No newdata: linear predictor and fitted values.
  expect_equal(predict(p$fit), predict(p$ref), tolerance = 1e-5)
  expect_equal(predict(p$fit, type = "response"), fitted(p$fit))

  # newdata, both scales, with standard errors.
  for (tp in c("link", "response")) {
    mine <- predict(p$fit, newdata = nd, type = tp, se.fit = TRUE)
    theirs <- predict(p$ref, newdata = nd, type = tp, se.fit = TRUE)
    expect_equal(mine$fit, theirs$fit, tolerance = 1e-5)
    expect_equal(mine$se.fit, theirs$se.fit, tolerance = 1e-5)
  }

  # newdata with factor levels given as character strings.
  nd2 <- data.frame(fem = "Women", mar = "Married", kid5 = 1, phd = 3,
                    ment = 10)
  expect_equal(predict(p$fit, newdata = nd2, type = "response"),
               predict(p$ref, newdata = nd2, type = "response"),
               tolerance = 1e-5)
})

test_that("predict uses offsets from the formula and the offset argument", {
  skip_if_not_installed("pscl")
  e <- new.env()
  utils::data("bioChemists", package = "pscl", envir = e)
  d <- e$bioChemists
  d$expo <- log(d$phd)
  nd <- d[1:5, ]

  fit_f <- nbreg(art ~ ment + offset(expo), data = d)
  fit_a <- nbreg(art ~ ment, data = d, offset = expo)

  manual <- as.vector(cbind(1, nd$ment) %*% coef(fit_f) + nd$expo)
  expect_equal(unname(predict(fit_f, newdata = nd)), manual,
               tolerance = 1e-10)
  expect_equal(unname(predict(fit_a, newdata = nd)), manual,
               tolerance = 1e-8)
})

test_that("summary gives the glm.nb coefficient table", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  p <- fit_pair()
  s <- summary(p$fit)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(s$coefficients, summary(p$ref)$coefficients, tolerance = 1e-5)
  expect_identical(colnames(s$coefficients),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(s$theta, p$ref$theta, tolerance = 1e-5)
  expect_equal(s$SE.theta, p$ref$SE.theta, tolerance = 1e-4)
  expect_equal(s$aic, AIC(p$ref), tolerance = 1e-8)
})

test_that("print methods produce output and return invisibly", {
  skip_if_not_installed("pscl")
  e <- new.env()
  utils::data("bioChemists", package = "pscl", envir = e)
  p <- list(fit = nbreg(art ~ ment, data = e$bioChemists))
  expect_output(out <- print(p$fit), "Theta")
  expect_identical(out, p$fit)
  s <- summary(p$fit)
  expect_output(out <- print(s), "Pr\\(>\\|z\\|\\)")
  expect_identical(out, s)
  expect_output(print(s), "Log\\(theta\\)")
})
