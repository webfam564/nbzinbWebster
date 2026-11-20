# Negative binomial (NB2) regression fitted by IRLS for beta alternated with
# Newton-Raphson for log(theta).  Implements Section 3.1 of the project
# description.  No general-purpose optimizer is used anywhere in this file.
#
# Notation (Section 2.1):
#   y_i       count response
#   x_i       row i of the model matrix X (p columns)
#   o_i       offset
#   a_i       prior weight (1 for an ordinary fit)
#   eta_i     = x_i' beta + o_i            (linear predictor, log link)
#   mu_i      = exp(eta_i)                 (mean)
#   V(mu)     = mu + mu^2 / theta          (NB2 variance function)

# ---------------------------------------------------------------------------
# Control
# ---------------------------------------------------------------------------

#' Control parameters for `nbreg()`
#'
#' Tolerances and iteration limits for the three loops used by [nbreg()]:
#' the outer loop that alternates the two updates, the inner IRLS loop for
#' the regression coefficients, and the inner Newton-Raphson loop for
#' \eqn{\log\theta}.
#'
#' The outer loop stops when the absolute change in the log-likelihood is
#' below `tol`. The IRLS loop uses the relative change
#' \eqn{|\ell_{new} - \ell_{old}| / (|\ell_{new}| + 0.1)} < `irls_tol`. The
#' Newton-Raphson loop stops when the step in \eqn{\log\theta} is smaller
#' than `nr_tol`.
#'
#' @param tol Outer-loop tolerance on the absolute change in the
#'   log-likelihood. Default `1e-8`.
#' @param maxit Maximum number of outer iterations. Default `100`.
#' @param irls_tol IRLS tolerance on the relative change in the
#'   log-likelihood. Default `1e-10`.
#' @param irls_maxit Maximum number of IRLS iterations per outer iteration.
#'   Default `100`.
#' @param nr_tol Newton-Raphson tolerance on the absolute step in
#'   \eqn{\log\theta}. Default `1e-10`.
#' @param nr_maxit Maximum number of Newton-Raphson iterations per outer
#'   iteration. Default `100`.
#' @param theta_max Upper bound for \eqn{\theta}. Data with no
#'   overdispersion push \eqn{\theta} towards infinity (the Poisson limit);
#'   the bound keeps the computations numerically stable. Default `1e8`.
#' @param theta_warn A warning is issued when the final \eqn{\hat\theta}
#'   exceeds this value, because the model is then essentially Poisson.
#'   Default `1e4`.
#' @param trace Logical. If `TRUE`, print the log-likelihood and
#'   \eqn{\theta} after each outer iteration.
#'
#' @return A list of control parameters.
#' @seealso [nbreg()]
#' @examples
#' nbreg_control(tol = 1e-10, trace = TRUE)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100,
                          irls_tol = 1e-10, irls_maxit = 100,
                          nr_tol = 1e-10, nr_maxit = 100,
                          theta_max = 1e8, theta_warn = 1e4,
                          trace = FALSE) {
  check_pos <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(sprintf("'%s' must be a single positive number.", name),
           call. = FALSE)
    }
  }
  check_pos(tol, "tol")
  check_pos(irls_tol, "irls_tol")
  check_pos(nr_tol, "nr_tol")
  check_pos(theta_max, "theta_max")
  check_pos(theta_warn, "theta_warn")
  check_pos(maxit, "maxit")
  check_pos(irls_maxit, "irls_maxit")
  check_pos(nr_maxit, "nr_maxit")
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    stop("'trace' must be TRUE or FALSE.", call. = FALSE)
  }
  list(tol = tol, maxit = as.integer(maxit),
       irls_tol = irls_tol, irls_maxit = as.integer(irls_maxit),
       nr_tol = nr_tol, nr_maxit = as.integer(nr_maxit),
       theta_max = theta_max, theta_warn = theta_warn,
       trace = trace)
}

# ---------------------------------------------------------------------------
# User-facing function
# ---------------------------------------------------------------------------

#' Negative binomial regression
#'
#' Fits a negative binomial (NB2) regression model with log link,
#' \eqn{\log\mu_i = x_i^\top\beta + o_i} and
#' \eqn{\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta}.
#' The coefficients are estimated by iteratively reweighted least squares
#' (IRLS) with \eqn{\theta} held fixed, alternated with Newton-Raphson
#' updates of \eqn{\log\theta} with \eqn{\beta} held fixed, until the
#' log-likelihood stops changing. No general-purpose optimizer is used.
#'
#' @section Algorithm:
#' 1. **Starting values.** \eqn{\beta} from a Poisson fit computed by IRLS
#'    (working weights \eqn{a_i\mu_i}); \eqn{\theta} from the moment
#'    estimate \eqn{\sum a_i\mu_i^2 / \sum a_i\{(y_i-\mu_i)^2-\mu_i\}}.
#' 2. **IRLS for \eqn{\beta}** (\eqn{\theta} fixed): with
#'    \eqn{w_i = a_i\mu_i\theta/(\theta+\mu_i)} and
#'    \eqn{z_i = (\eta_i - o_i) + (y_i-\mu_i)/\mu_i},
#'    update \eqn{\beta \leftarrow (X^\top W X)^{-1}X^\top W z}.
#' 3. **Newton-Raphson for \eqn{\phi=\log\theta}** (\eqn{\beta} fixed):
#'    \eqn{\phi \leftarrow \phi - U_\phi/H_\phi} with
#'    \eqn{U_\phi = \theta U(\theta)} and
#'    \eqn{H_\phi = \theta^2 H(\theta) + \theta U(\theta)}.
#' 4. Repeat 2-3 until the log-likelihood changes by less than
#'    `control$tol`, then do a final IRLS pass at the final \eqn{\theta}.
#'
#' Standard errors: \eqn{\widehat{\mathrm{Var}}(\hat\beta) =
#' (X^\top\hat W X)^{-1}} (treating \eqn{\theta} as known) and
#' \eqn{\mathrm{SE}(\hat\theta) = \{-H(\hat\theta)\}^{-1/2}}.
#'
#' @param formula A model formula, `response ~ terms`. Offsets may be given
#'   with `offset()` terms in the formula.
#' @param data An optional data frame containing the variables in the model.
#' @param offset An optional vector of offsets (e.g. log exposure), added to
#'   any `offset()` terms in the formula.
#' @param control A list of control parameters; see [nbreg_control()].
#' @param weights Optional non-negative prior weights \eqn{a_i}.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values (default: the
#'   `na.action` option, usually [stats::na.omit()]).
#'
#' @return An object of class `"nbreg"`, a list with components
#' \describe{
#'   \item{coefficients}{Named vector of estimated \eqn{\beta}.}
#'   \item{vcov}{Estimated covariance matrix of \eqn{\hat\beta}.}
#'   \item{theta, SE.theta}{\eqn{\hat\theta} and its standard error.}
#'   \item{log_theta, SE.log_theta}{\eqn{\log\hat\theta} and its standard
#'     error (delta method).}
#'   \item{loglik}{Maximized log-likelihood.}
#'   \item{df}{Number of estimated parameters (\eqn{p + 1}).}
#'   \item{deviance}{Residual deviance at \eqn{\hat\theta}.}
#'   \item{iter}{Number of outer iterations.}
#'   \item{converged}{Logical convergence flag.}
#'   \item{fitted.values, linear.predictors, residuals}{\eqn{\hat\mu_i},
#'     \eqn{\hat\eta_i}, and response residuals \eqn{y_i - \hat\mu_i}.}
#'   \item{weights, prior.weights}{Final IRLS working weights and the prior
#'     weights.}
#'   \item{y, offset, nobs, df.residual}{Response, offset, number of
#'     observations with positive weight, and \eqn{n - p}.}
#'   \item{call, formula, terms, model, xlevels, contrasts, na.action,
#'     control}{Information about the fit.}
#' }
#'
#' @references
#' Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics
#' with S* (4th ed.). Springer.
#'
#' @seealso [nbreg_control()]
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- rnorm(n)
#' y <- rnbinom(n, mu = exp(0.5 + 0.8 * x), size = 2)
#' fit <- nbreg(y ~ x)
#' fit$coefficients
#' fit$theta
#' fit$loglik
#'
#' # with an exposure offset
#' t <- runif(n, 1, 5)
#' y2 <- rnbinom(n, mu = t * exp(0.5 + 0.8 * x), size = 2)
#' nbreg(y2 ~ x, offset = log(t))$coefficients
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control(),
                  weights, subset, na.action) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)

  # Build the model frame with the standard model.frame() machinery, so that
  # `data`, `subset`, `weights`, `offset` and `na.action` behave as usual.
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  # Validate the control list by passing it through nbreg_control().
  if (!is.list(control)) {
    stop("'control' must be a list, usually from nbreg_control().",
         call. = FALSE)
  }
  control <- do.call(nbreg_control, control)

  # Response: non-negative integer counts.
  y <- stats::model.response(mf, "numeric")
  if (is.null(y)) stop("The formula must have a response.", call. = FALSE)
  if (is.matrix(y)) stop("The response must be a vector.", call. = FALSE)
  if (any(!is.finite(y)) || any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    stop("The response must contain non-negative integer counts.",
         call. = FALSE)
  }
  y <- as.vector(y)
  n <- length(y)

  # Model matrix X.
  X <- stats::model.matrix(mt, mf)
  p <- ncol(X)
  if (p == 0L) stop("The model must contain at least one coefficient.",
                    call. = FALSE)
  if (qr(X)$rank < p) {
    stop("The model matrix is rank deficient (some coefficients are not ",
         "identifiable). Remove redundant covariates.", call. = FALSE)
  }

  # Prior weights a_i.
  a <- stats::model.weights(mf)
  if (is.null(a)) a <- rep(1, n)
  if (!is.numeric(a) || any(!is.finite(a)) || any(a < 0)) {
    stop("'weights' must be non-negative and finite.", call. = FALSE)
  }
  if (sum(a > 0) <= p) {
    stop("Not enough observations with positive weight to fit the model.",
         call. = FALSE)
  }
  if (all(y[a > 0] == 0)) {
    stop("All responses are zero; the model cannot be fitted.",
         call. = FALSE)
  }

  # Offset o_i: model.offset() sums the offset() terms in the formula and
  # the `offset` argument (stored in the model frame as "(offset)").
  o <- stats::model.offset(mf)
  if (is.null(o)) o <- rep(0, n)
  o <- as.vector(o)
  if (length(o) != n || any(!is.finite(o))) {
    stop("'offset' must be finite and have one value per observation.",
         call. = FALSE)
  }

  fit <- nb_fit(X, y, o, a, control)

  if (!fit$converged) {
    warning(sprintf(paste0("nbreg: algorithm did not converge in %d outer ",
                           "iterations; increase 'maxit' in nbreg_control()."),
                    control$maxit), call. = FALSE)
  }
  if (fit$theta > control$theta_warn) {
    warning(sprintf(paste0("nbreg: theta is very large (%.3g); the data show ",
                           "no overdispersion and the model is essentially ",
                           "Poisson."), fit$theta), call. = FALSE)
  }

  names(fit$coefficients) <- colnames(X)
  dimnames(fit$vcov) <- list(colnames(X), colnames(X))
  rn <- rownames(mf)
  names(fit$fitted.values) <- names(fit$linear.predictors) <- rn
  names(fit$residuals) <- names(fit$weights) <- rn
  names(y) <- rn

  nobs <- sum(a > 0)
  out <- c(fit, list(
    df = p + 1L,
    y = y,
    offset = o,
    prior.weights = a,
    nobs = nobs,
    df.residual = nobs - p,
    call = cl,
    formula = formula,
    terms = mt,
    model = mf,
    xlevels = stats::.getXlevels(mt, mf),
    contrasts = attr(X, "contrasts"),
    na.action = attr(mf, "na.action"),
    control = control
  ))
  class(out) <- "nbreg"
  out
}

# ---------------------------------------------------------------------------
# Internal workhorse (also reused by the ZINB M-step with a_i = 1 - tau_i)
# ---------------------------------------------------------------------------

# Fit an NB2 model to a model matrix.  Section 3.1 (a)-(d).
#   X: n x p model matrix; y: counts; o: offset; a: prior weights.
#   beta_start / theta_start: optional starting values (used by EM).
nb_fit <- function(X, y, o, a, control,
                   beta_start = NULL, theta_start = NULL) {

  # --- (c) Starting values ------------------------------------------------
  # beta: a Poisson fit by IRLS (theta = Inf, so w_i = a_i mu_i).
  if (is.null(beta_start)) {
    beta <- nb_irls(X, y, o, a, theta = Inf, beta = NULL, control)$beta
  } else {
    beta <- beta_start
  }
  mu <- as.vector(exp(X %*% beta + o))

  # theta: moment estimate from Var(Y) = mu + mu^2 / theta, i.e.
  #   theta = sum a mu^2 / sum a {(y - mu)^2 - mu}.
  # If the data look under-dispersed the denominator is <= 0; then start at 1
  # and let Newton-Raphson move theta upward.
  if (is.null(theta_start)) {
    theta <- sum(a * mu^2) / sum(a * ((y - mu)^2 - mu))
    if (!is.finite(theta) || theta <= 0) theta <- 1
    theta <- min(theta, control$theta_max)
  } else {
    theta <- theta_start
  }

  ll_old <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  iter <- 0L

  # --- Outer loop: alternate (a) and (b) ---------------------------------
  for (iter in seq_len(control$maxit)) {
    # (a) IRLS for beta with theta fixed.
    irls <- nb_irls(X, y, o, a, theta, beta, control)
    beta <- irls$beta
    mu <- irls$mu

    # (b) Newton-Raphson for log(theta) with beta (hence mu) fixed.
    theta <- nb_theta_nr(y, mu, a, theta, control)

    ll <- nb_loglik(y, mu, theta, a)
    if (control$trace) {
      cat(sprintf("nbreg iter %3d: logLik = %.10f  theta = %.8g\n",
                  iter, ll, theta))
    }
    # Stop when the log-likelihood changes by less than the tolerance.
    if (abs(ll - ll_old) < control$tol) {
      converged <- TRUE
      break
    }
    ll_old <- ll
  }

  # Final IRLS pass at the final theta.  The loop ends right after a theta
  # update, so beta was last fitted at the previous theta; this pass makes
  # beta the maximizer for the theta we report (the beta score is then ~0).
  irls <- nb_irls(X, y, o, a, theta, beta, control)
  beta <- irls$beta

  # --- (d) Standard errors at convergence --------------------------------
  eta <- as.vector(X %*% beta + o)
  mu <- exp(eta)
  # Working weights w_i = a_i mu_i theta / (theta + mu_i) at the MLE.
  w <- a * mu * theta / (theta + mu)
  # Var(beta_hat) = (X' W X)^{-1}, theta treated as known.
  vcov <- solve(crossprod(X, w * X))
  # SE(theta_hat) = {-H(theta_hat)}^{-1/2}; SE(log theta) by the delta method.
  H <- nb_hess_theta(y, mu, theta, a)
  se_theta <- if (H < 0) 1 / sqrt(-H) else NA_real_

  list(
    coefficients = as.vector(beta),
    vcov = vcov,
    theta = theta,
    SE.theta = se_theta,
    log_theta = log(theta),
    SE.log_theta = se_theta / theta,
    loglik = nb_loglik(y, mu, theta, a),
    deviance = nb_deviance(y, mu, theta, a),
    iter = iter,
    converged = converged,
    fitted.values = mu,
    linear.predictors = eta,
    residuals = y - mu,
    weights = w
  )
}

# (a) IRLS for beta with theta fixed (theta = Inf gives the Poisson fit).
#   w_i = a_i (dmu/deta)^2 / V(mu_i) = a_i mu_i theta / (theta + mu_i)
#   z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
#   beta^(t+1) = (X' W X)^{-1} X' W z
# If beta = NULL, start from mu_i = y_i + 0.1 (so log(mu) is finite).
nb_irls <- function(X, y, o, a, theta, beta, control) {
  if (is.null(beta)) {
    mu <- y + 0.1
    eta <- log(mu)
  } else {
    eta <- as.vector(X %*% beta + o)
    mu <- exp(eta)
  }
  ll_old <- nb_loglik(y, mu, theta, a)
  converged <- FALSE

  for (it in seq_len(control$irls_maxit)) {
    # Working weights.  For the Poisson start (theta = Inf), theta/(theta+mu)
    # equals 1, so w_i = a_i mu_i.
    w <- if (is.infinite(theta)) a * mu else a * mu * theta / (theta + mu)
    # Working response.
    z <- (eta - o) + (y - mu) / mu
    # Weighted least-squares step: solve (X'WX) beta = X'Wz.
    XtW <- t(X * w)
    beta_new <- as.vector(solve(XtW %*% X, XtW %*% z))
    eta_new <- as.vector(X %*% beta_new + o)
    mu_new <- exp(eta_new)
    ll_new <- nb_loglik(y, mu_new, theta, a)

    # Safeguard (step halving): Fisher scoring normally increases the
    # log-likelihood; if a step decreases it (or overflows), halve the step
    # back towards the current beta.  Not used on the very first step from
    # mu = y + 0.1, because that start is not a point of the form X beta.
    if (!is.null(beta)) {
      halve <- 0L
      while ((!is.finite(ll_new) ||
              ll_new < ll_old - 1e-12 * (abs(ll_old) + 0.1)) &&
             halve < 30L) {
        beta_new <- (beta + beta_new) / 2
        eta_new <- as.vector(X %*% beta_new + o)
        mu_new <- exp(eta_new)
        ll_new <- nb_loglik(y, mu_new, theta, a)
        halve <- halve + 1L
      }
    }

    done <- !is.null(beta) &&
      abs(ll_new - ll_old) < control$irls_tol * (abs(ll_new) + 0.1)
    beta <- beta_new
    eta <- eta_new
    mu <- mu_new
    ll_old <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }
  list(beta = beta, mu = mu, iter = it, converged = converged)
}

# (b) Newton-Raphson for phi = log(theta) with mu fixed.
#   U_phi = theta U(theta),   H_phi = theta^2 H(theta) + theta U(theta)
#   phi <- phi - U_phi / H_phi
# Returns the updated theta.
nb_theta_nr <- function(y, mu, a, theta, control) {
  phi <- log(theta)
  phi_max <- log(control$theta_max)
  ll <- nb_loglik(y, mu, theta, a)

  for (it in seq_len(control$nr_maxit)) {
    U <- nb_score_theta(y, mu, theta, a)
    H <- nb_hess_theta(y, mu, theta, a)
    U_phi <- theta * U                  # chain rule: d l / d phi
    H_phi <- theta^2 * H + theta * U    # d^2 l / d phi^2
    # Newton step.  If the log-likelihood is not concave in phi here
    # (H_phi >= 0) the Newton step would go downhill, so take a unit step
    # in the uphill direction instead and let step halving shorten it.
    step <- if (H_phi < 0) -U_phi / H_phi else sign(U_phi)

    # Step halving: never accept a step that lowers the log-likelihood.
    # phi is capped at log(theta_max) (Poisson limit).
    phi_new <- min(phi + step, phi_max)
    ll_new <- nb_loglik(y, mu, exp(phi_new), a)
    halve <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-12 * (abs(ll) + 0.1)) &&
           halve < 30L) {
      step <- step / 2
      phi_new <- min(phi + step, phi_max)
      ll_new <- nb_loglik(y, mu, exp(phi_new), a)
      halve <- halve + 1L
    }
    if (!is.finite(ll_new) || ll_new < ll - 1e-12 * (abs(ll) + 0.1)) break

    change <- abs(phi_new - phi)
    phi <- phi_new
    theta <- exp(phi)
    ll <- ll_new
    if (change < control$nr_tol) break
  }
  theta
}

# ---------------------------------------------------------------------------
# NB log-likelihood and its derivatives in theta
# ---------------------------------------------------------------------------

# Log-likelihood (Section 2.1):
#   l = sum a_i [ lgamma(y_i + theta) - lgamma(theta) - lgamma(y_i + 1)
#                 + theta log(theta) + y_i log(mu_i)
#                 - (theta + y_i) log(theta + mu_i) ]
# theta = Inf gives the Poisson log-likelihood
#   sum a_i [ y_i log(mu_i) - mu_i - lgamma(y_i + 1) ].
# Observations with a_i = 0 are dropped so that they cannot produce NaN.
nb_loglik <- function(y, mu, theta, a) {
  keep <- a > 0
  y <- y[keep]; mu <- mu[keep]; a <- a[keep]
  # y log(mu) is 0 when y = 0 (even if mu underflows to 0).
  ylogmu <- ifelse(y > 0, y * log(mu), 0)
  if (is.infinite(theta)) {
    return(sum(a * (ylogmu - mu - lgamma(y + 1))))
  }
  sum(a * (lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
             theta * log(theta) + ylogmu -
             (theta + y) * log(theta + mu)))
}

# Score in theta (Section 3.1 (b)):
#   U(theta) = sum a_i [ psi(y_i + theta) - psi(theta) + log(theta) + 1
#                        - log(theta + mu_i) - (y_i + theta)/(theta + mu_i) ]
nb_score_theta <- function(y, mu, theta, a) {
  sum(a * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
             log(theta + mu) - (y + theta) / (theta + mu)))
}

# Second derivative in theta (Section 3.1 (b)):
#   H(theta) = sum a_i [ psi'(y_i + theta) - psi'(theta) + 1/theta
#                        - 2/(theta + mu_i) + (y_i + theta)/(theta + mu_i)^2 ]
nb_hess_theta <- function(y, mu, theta, a) {
  sum(a * (trigamma(y + theta) - trigamma(theta) + 1 / theta -
             2 / (theta + mu) + (y + theta) / (theta + mu)^2))
}

# Residual deviance at fixed theta:
#   D = 2 sum a_i [ y_i log(y_i / mu_i)
#                   - (y_i + theta) log((y_i + theta) / (mu_i + theta)) ]
nb_deviance <- function(y, mu, theta, a) {
  ylog <- ifelse(y > 0, y * log(y / mu), 0)
  2 * sum(a * (ylog - (y + theta) * log((y + theta) / (mu + theta))))
}
