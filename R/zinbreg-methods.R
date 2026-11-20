# S3 methods for objects of class "zinbreg".

#' Methods for `zinbreg` objects
#'
#' Standard extractor, summary, and prediction methods for zero-inflated
#' negative binomial fits returned by [zinbreg()].
#'
#' * `coef()`: \eqn{\hat\beta} and \eqn{\hat\gamma}, either together
#'   (`model = "full"`, named `count_*` and `zero_*` as in
#'   `pscl::zeroinfl()`) or one part at a time.
#' * `vcov()`: the matching block of the inverse observed information
#'   (without \eqn{\log\theta}, whose SE is reported separately).
#' * `logLik()`: the maximized observed-data log-likelihood with attribute
#'   `df` = \eqn{p + q + 1}, so [stats::AIC()] and [stats::BIC()] work.
#' * `fitted()`: the fitted mean \eqn{(1-\hat\pi_i)\hat\mu_i}.
#' * `residuals()`: response residuals \eqn{y_i - (1-\hat\pi_i)\hat\mu_i}, or
#'   Pearson residuals that divide them by \eqn{\sqrt{\mathrm{Var}(Y_i)}},
#'   \eqn{\mathrm{Var}(Y_i) = (1-\pi_i)\mu_i(1 + \pi_i\mu_i + \mu_i/\theta)}.
#' * `predict()`: `type = "response"` gives \eqn{(1-\hat\pi)\hat\mu},
#'   `"count"` gives \eqn{\hat\mu}, `"zero"` gives \eqn{\hat\pi}, and
#'   `"prob"` gives the matrix of \eqn{P(Y_i = k)} for `k` in `at`.
#' * `summary()`: coefficient tables for both parts with Wald z tests; the
#'   count table also shows \eqn{\log\theta}.
#'
#' @param object,x An object of class `"zinbreg"` (or `"summary.zinbreg"`
#'   for the print method of a summary).
#' @param model For `coef()` and `vcov()`: `"full"`, `"count"`, or `"zero"`.
#' @param type For `residuals()`: `"response"` or `"pearson"`. For
#'   `predict()`: `"response"`, `"count"`, `"zero"`, or `"prob"`.
#' @param newdata Optional data frame in which to look for the covariates
#'   (and any count-offset variables). If omitted, the fitted values are
#'   used.
#' @param at For `predict(type = "prob")`: the counts at which to evaluate
#'   the probabilities. Default `0:max(y)`.
#' @param digits Number of significant digits to print.
#' @param signif.stars Logical; show significance stars?
#' @param ... Further arguments (currently ignored).
#'
#' @return `coef()`, `fitted()`, `residuals()` return numeric vectors;
#'   `vcov()` a matrix; `logLik()` an object of class `"logLik"`; `nobs()`
#'   an integer; `predict()` a vector, or a matrix for `type = "prob"`;
#'   `summary()` an object of class `"summary.zinbreg"`. The print methods
#'   return their argument invisibly.
#'
#' @examples
#' set.seed(2)
#' n <- 500
#' x <- rnorm(n)
#' z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(-1 + z)) == 1, 0,
#'             rnbinom(n, mu = exp(1 + 0.5 * x), size = 2))
#' fit <- zinbreg(y ~ x | z)
#' coef(fit)
#' coef(fit, model = "zero")
#' sqrt(diag(vcov(fit)))
#' logLik(fit)
#' AIC(fit)
#' nd <- data.frame(x = c(-1, 0, 1), z = c(0, 0, 2))
#' predict(fit, newdata = nd, type = "response")
#' predict(fit, newdata = nd, type = "zero")
#' # expected number of zeros
#' sum(predict(fit, type = "prob", at = 0))
#' @name zinbreg-methods
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                          ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Count model coefficients (negbin with log link):\n")
  print.default(format(x$coefficients$count, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("Theta =", format(x$theta, digits = digits), "\n\n")
  cat("Zero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nLog-likelihood:", format(x$loglik, digits = max(digits, 6L)),
      "on", x$df, "df\n")
  if (!x$converged) cat("Warning: EM did not converge.\n")
  cat("\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)

  # Wald table: estimate, SE, z = estimate / SE, two-sided normal p-value.
  wald <- function(est, se) {
    z <- est / se
    cbind(Estimate = est, `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  }
  # Count table includes log(theta), as in pscl::zeroinfl().
  count_tab <- wald(c(object$coefficients$count, `Log(theta)` =
                        object$log_theta),
                    se[c(seq_len(p), p + q + 1L)])
  zero_tab <- wald(object$coefficients$zero, se[p + seq_len(q)])
  rownames(count_tab) <- c(names(object$coefficients$count), "Log(theta)")
  rownames(zero_tab) <- names(object$coefficients$zero)

  out <- list(
    call = object$call,
    coefficients = list(count = count_tab, zero = zero_tab),
    pearson = stats::residuals(object, type = "pearson"),
    theta = object$theta,
    SE.theta = object$SE.theta,
    loglik = object$loglik,
    df = object$df,
    aic = -2 * object$loglik + 2 * object$df,
    iter = object$iter,
    converged = object$converged
  )
  class(out) <- "summary.zinbreg"
  out
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars =
                                    getOption("show.signif.stars"),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson, names = FALSE)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCount model coefficients (negbin with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits,
                      signif.stars = signif.stars, na.print = "NA")
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits,
                      signif.stars = signif.stars, na.print = "NA")
  cat("\nTheta =", format(x$theta, digits = digits),
      " (SE", paste0(format(x$SE.theta, digits = digits), ")"), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = max(digits, 6L)),
      "on", x$df, "df   AIC:", format(x$aic, digits = max(digits, 6L)), "\n")
  cat("EM iterations:", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    count = cf$count,
    zero = cf$zero,
    full = stats::setNames(c(cf$count, cf$zero),
                           c(paste0("count_", names(cf$count)),
                             paste0("zero_", names(cf$zero))))
  )
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  idx <- switch(model,
    full = seq_len(p + q),
    count = seq_len(p),
    zero = p + seq_len(q)
  )
  v <- object$vcov[idx, idx, drop = FALSE]
  if (model != "full") {
    nm <- names(object$coefficients[[model]])
    dimnames(v) <- list(nm, nm)
  }
  v
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) {
  object$nobs
}

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) {
  object$fitted.values
}

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"),
                              ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$mu
  p_zero <- object$pi
  theta <- object$theta
  r <- y - (1 - p_zero) * mu
  switch(type,
    response = r,
    # Var(Y_i) = (1 - pi_i) mu_i (1 + pi_i mu_i + mu_i / theta)
    pearson = sqrt(object$prior.weights) * r /
      sqrt((1 - p_zero) * mu * (1 + p_zero * mu + mu / theta))
  )
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero", "prob"),
                            at = NULL, ...) {
  type <- match.arg(type)

  if (missing(newdata) || is.null(newdata)) {
    mu <- object$mu
    p_zero <- object$pi
  } else {
    # Rebuild both model matrices for newdata with the terms, factor levels
    # and contrasts of the fit.
    ttX <- stats::delete.response(object$terms$count)
    ttZ <- stats::delete.response(object$terms$zero)
    mfX <- stats::model.frame(ttX, newdata, na.action = stats::na.pass,
                              xlev = object$xlevels$count)
    mfZ <- stats::model.frame(ttZ, newdata, na.action = stats::na.pass,
                              xlev = object$xlevels$zero)
    X <- stats::model.matrix(ttX, mfX, contrasts.arg = object$contrasts$count)
    Z <- stats::model.matrix(ttZ, mfZ, contrasts.arg = object$contrasts$zero)
    # Count offset: offset() terms plus the original `offset` argument
    # re-evaluated in newdata.
    o <- stats::model.offset(mfX)
    if (is.null(o)) o <- rep(0, nrow(X))
    if (!is.null(object$call$offset)) {
      o <- o + eval(object$call$offset, newdata,
                    environment(object$terms$count))
    }
    mu <- as.vector(exp(X %*% object$coefficients$count + o))
    p_zero <- as.vector(stats::plogis(Z %*% object$coefficients$zero))
    names(mu) <- names(p_zero) <- rownames(mfX)
  }

  switch(type,
    response = (1 - p_zero) * mu,
    count = mu,
    zero = p_zero,
    prob = {
      # P(Y = k) = pi I(k = 0) + (1 - pi) f_NB(k; mu, theta)
      if (is.null(at)) at <- 0:max(object$y)
      prob <- vapply(at, function(k) {
        (1 - p_zero) * stats::dnbinom(k, size = object$theta, mu = mu) +
          p_zero * (k == 0)
      }, numeric(length(mu)))
      prob <- matrix(prob, nrow = length(mu))
      dimnames(prob) <- list(names(mu), at)
      prob
    }
  )
}
