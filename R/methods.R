# S3 methods for objects of class "nbreg".

#' Methods for `nbreg` objects
#'
#' Standard extractor, summary, and prediction methods for negative binomial
#' fits returned by [nbreg()].
#'
#' * `coef()`, `vcov()`: estimated \eqn{\beta} and its covariance matrix
#'   \eqn{(X^\top\hat W X)^{-1}} (with \eqn{\theta} treated as known).
#' * `logLik()`: the maximized log-likelihood with attribute `df` =
#'   \eqn{p + 1} (the coefficients plus \eqn{\theta}), so [stats::AIC()] and
#'   [stats::BIC()] work.
#' * `fitted()`: fitted means \eqn{\hat\mu_i}.
#' * `residuals()`: response residuals \eqn{y_i - \hat\mu_i}, Pearson
#'   residuals \eqn{\sqrt{a_i}(y_i - \hat\mu_i)/\sqrt{V(\hat\mu_i)}} with
#'   \eqn{V(\mu) = \mu + \mu^2/\theta}, or deviance residuals.
#' * `predict()`: the linear predictor \eqn{\hat\eta} (`type = "link"`) or
#'   the mean \eqn{\hat\mu = e^{\hat\eta}} (`type = "response"`), for the
#'   original data or for `newdata`, optionally with standard errors.
#' * `summary()`: coefficient table with Wald z tests, plus
#'   \eqn{\hat\theta}, its SE, the log-likelihood, AIC, and deviance.
#'
#' @param object,x An object of class `"nbreg"` (or `"summary.nbreg"` for
#'   the print method of a summary).
#' @param type For `residuals()`: `"response"`, `"pearson"`, or
#'   `"deviance"`. For `predict()`: `"link"` or `"response"`.
#' @param newdata Optional data frame in which to look for the covariates
#'   (and any offset variables). If omitted, the fitted values are used.
#' @param se.fit Logical; if `TRUE`, `predict()` also returns standard
#'   errors (delta method on the response scale).
#' @param digits Number of significant digits to print.
#' @param signif.stars Logical; show significance stars?
#' @param ... Further arguments (currently ignored).
#'
#' @return `coef()`, `fitted()`, `residuals()` return numeric vectors;
#'   `vcov()` a matrix; `logLik()` an object of class `"logLik"`; `nobs()`
#'   an integer; `predict()` a vector, or with `se.fit = TRUE` a list with
#'   components `fit` and `se.fit`; `summary()` an object of class
#'   `"summary.nbreg"`. The print methods return their argument invisibly.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- rnorm(n)
#' y <- rnbinom(n, mu = exp(0.5 + 0.8 * x), size = 2)
#' fit <- nbreg(y ~ x)
#' fit
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#' AIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "response")
#' @name nbreg-methods
#' @importFrom stats coef vcov logLik nobs fitted residuals predict
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients (negative binomial with log link):\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Log-likelihood:", format(x$loglik, digits = max(digits, 6L)),
      "on", x$df, "df\n")
  if (!x$converged) cat("Warning: the algorithm did not converge.\n")
  cat("\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  # Wald z statistic and two-sided normal p-value.
  z <- est / se
  p <- 2 * stats::pnorm(-abs(z))
  coef_table <- cbind(Estimate = est, `Std. Error` = se,
                      `z value` = z, `Pr(>|z|)` = p)

  out <- list(
    call = object$call,
    coefficients = coef_table,
    pearson = stats::residuals(object, type = "pearson"),
    theta = object$theta,
    SE.theta = object$SE.theta,
    log_theta = object$log_theta,
    SE.log_theta = object$SE.log_theta,
    loglik = object$loglik,
    df = object$df,
    aic = -2 * object$loglik + 2 * object$df,
    deviance = object$deviance,
    df.residual = object$df.residual,
    iter = object$iter,
    converged = object$converged
  )
  class(out) <- "summary.nbreg"
  out
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson, names = FALSE)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCoefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, na.print = "NA")
  cat("\nTheta:", format(x$theta, digits = digits),
      " (SE", paste0(format(x$SE.theta, digits = digits), ")"),
      "   Log(theta):", format(x$log_theta, digits = digits),
      " (SE", paste0(format(x$SE.log_theta, digits = digits), ")"), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = max(digits, 6L)),
      "on", x$df, "df   AIC:", format(x$aic, digits = max(digits, 6L)), "\n")
  cat("Residual deviance:", format(x$deviance, digits = max(digits, 5L)),
      "on", x$df.residual, "degrees of freedom\n")
  cat("Outer iterations (IRLS + Newton-Raphson):", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) {
  object$coefficients
}

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) {
  object$vcov
}

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) {
  object$nobs
}

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) {
  object$fitted.values
}

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object,
                            type = c("response", "pearson", "deviance"),
                            ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  theta <- object$theta
  a <- object$prior.weights
  switch(type,
    # r_i = y_i - mu_i
    response = y - mu,
    # r_i = sqrt(a_i) (y_i - mu_i) / sqrt(mu_i + mu_i^2 / theta)
    pearson = sqrt(a) * (y - mu) / sqrt(mu + mu^2 / theta),
    # r_i = sign(y_i - mu_i) sqrt(d_i), where d_i is observation i's
    # contribution to the deviance (see nb_deviance()).
    deviance = {
      ylog <- ifelse(y > 0, y * log(y / mu), 0)
      d <- 2 * a * (ylog - (y + theta) * log((y + theta) / (mu + theta)))
      sign(y - mu) * sqrt(pmax(d, 0))
    }
  )
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("link", "response"),
                          se.fit = FALSE, ...) {
  type <- match.arg(type)

  if (missing(newdata) || is.null(newdata)) {
    X <- stats::model.matrix(object$terms, object$model,
                             contrasts.arg = object$contrasts)
    eta <- object$linear.predictors
  } else {
    # Rebuild the model matrix for newdata with the same terms, factor
    # levels and contrasts as the fit.
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                             xlev = object$xlevels)
    X <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    # Offsets: offset() terms in the formula, plus the `offset` argument of
    # the original call re-evaluated in newdata.
    o <- stats::model.offset(mf)
    if (is.null(o)) o <- rep(0, nrow(X))
    if (!is.null(object$call$offset)) {
      o <- o + eval(object$call$offset, newdata,
                    environment(object$terms))
    }
    eta <- as.vector(X %*% object$coefficients + o)
    names(eta) <- rownames(mf)
  }

  fit <- if (type == "link") eta else exp(eta)
  if (!se.fit) return(fit)

  # SE of eta_hat = x' beta_hat: sqrt(x' V x).  On the response scale,
  # delta method: SE(mu_hat) = mu_hat * SE(eta_hat).
  se_eta <- sqrt(rowSums((X %*% object$vcov) * X))
  se <- if (type == "link") se_eta else exp(eta) * se_eta
  names(se) <- names(fit)
  list(fit = fit, se.fit = se)
}
