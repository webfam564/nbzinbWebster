# Zero-inflated negative binomial (ZINB) regression fitted by the EM
# algorithm with IRLS M-steps.  Implements Section 3.2 of the project
# description.  No general-purpose optimizer is used anywhere in this file.
#
# Notation (Section 2.2):
#   y_i        count response, a_i prior weight (1 for an ordinary fit)
#   x_i, o_i   count-model covariates and offset
#   z_i        zero-model covariates
#   mu_i       = exp(x_i' beta + o_i)             (NB mean, log link)
#   pi_i       = plogis(zeta_i), zeta_i = z_i' gamma  (P(structural zero))
#   f_NB(0; mu, theta) = (theta / (theta + mu))^theta
#   D_i        latent indicator: 1 if observation i is a structural zero
#   tau_i      = P(D_i = 1 | y_i)                  (E-step)

# ---------------------------------------------------------------------------
# Control
# ---------------------------------------------------------------------------

#' Control parameters for `zinbreg()`
#'
#' Tolerances and iteration limits for the EM algorithm used by
#' [zinbreg()] and for the IRLS fits inside its M-steps.
#'
#' The EM algorithm stops when the observed-data log-likelihood changes by
#' less than `tol` (absolute change). Each M-step is fitted to convergence:
#' the zero-part logistic regression by IRLS (relative log-likelihood change
#' below `irls_tol`), and the count part by the NB algorithm of [nbreg()]
#' with outer tolerance `mstep_tol`.
#'
#' @param tol EM tolerance on the absolute change in the observed-data
#'   log-likelihood. Default `1e-8`.
#' @param maxit Maximum number of EM iterations. Default `1000`.
#' @param mstep_tol Outer tolerance (absolute change in the weighted NB
#'   log-likelihood) for the count-part M-step. Default `1e-10`.
#' @param mstep_maxit Maximum number of outer iterations in the count-part
#'   M-step. Default `100`.
#' @param irls_tol,irls_maxit Tolerance (relative change in the
#'   log-likelihood) and iteration limit for the IRLS loops (both parts).
#' @param nr_tol,nr_maxit Tolerance (step in \eqn{\log\theta}) and iteration
#'   limit for the Newton-Raphson update of \eqn{\theta}.
#' @param theta_max,theta_warn Upper bound for \eqn{\theta}, and the value
#'   above which a warning about the Poisson limit is given; see
#'   [nbreg_control()].
#' @param trace Logical. If `TRUE`, print the log-likelihood after each EM
#'   iteration.
#'
#' @return A list of control parameters.
#' @seealso [zinbreg()], [nbreg_control()]
#' @examples
#' zinbreg_control(tol = 1e-10, maxit = 2000)
#' @export
zinbreg_control <- function(tol = 1e-8, maxit = 1000,
                            mstep_tol = 1e-10, mstep_maxit = 100,
                            irls_tol = 1e-10, irls_maxit = 100,
                            nr_tol = 1e-10, nr_maxit = 100,
                            theta_max = 1e8, theta_warn = 1e4,
                            trace = FALSE) {
  check_pos(tol, "tol")
  check_pos(maxit, "maxit")
  check_pos(mstep_tol, "mstep_tol")
  check_pos(mstep_maxit, "mstep_maxit")
  check_pos(irls_tol, "irls_tol")
  check_pos(irls_maxit, "irls_maxit")
  check_pos(nr_tol, "nr_tol")
  check_pos(nr_maxit, "nr_maxit")
  check_pos(theta_max, "theta_max")
  check_pos(theta_warn, "theta_warn")
  check_flag(trace, "trace")
  list(tol = tol, maxit = as.integer(maxit),
       mstep_tol = mstep_tol, mstep_maxit = as.integer(mstep_maxit),
       irls_tol = irls_tol, irls_maxit = as.integer(irls_maxit),
       nr_tol = nr_tol, nr_maxit = as.integer(nr_maxit),
       theta_max = theta_max, theta_warn = theta_warn,
       trace = trace)
}

# ---------------------------------------------------------------------------
# User-facing function
# ---------------------------------------------------------------------------

#' Zero-inflated negative binomial regression
#'
#' Fits a zero-inflated negative binomial (ZINB) regression model. With
#' probability \eqn{\pi_i} observation \eqn{i} is a structural zero; with
#' probability \eqn{1-\pi_i} it comes from an NB(\eqn{\mu_i, \theta})
#' distribution, where \eqn{\log\mu_i = x_i^\top\beta + o_i} and
#' \eqn{\mathrm{logit}\,\pi_i = z_i^\top\gamma}. The model is fitted by an
#' EM algorithm whose M-steps are IRLS fits. No general-purpose optimizer
#' is used.
#'
#' @section Algorithm:
#' Starting values: \eqn{(\beta, \theta)} from an NB fit (as in
#' [nbreg()]); \eqn{\gamma} from a logistic regression of
#' \eqn{I(y_i = 0)} on \eqn{z_i}, fitted by IRLS. Then repeat:
#'
#' 1. **E-step.** Posterior probability of a structural zero,
#'    \eqn{\tau_i = \pi_i / \{\pi_i + (1-\pi_i) f_{NB}(0;\mu_i,\theta)\}}
#'    if \eqn{y_i = 0}, and \eqn{\tau_i = 0} if \eqn{y_i > 0}.
#' 2. **M-step, zero part.** Maximize
#'    \eqn{\sum_i a_i\{\tau_i\log\pi_i + (1-\tau_i)\log(1-\pi_i)\}}: a
#'    logistic regression with fractional response \eqn{\tau_i}, fitted by
#'    IRLS with \eqn{w_i = a_i\pi_i(1-\pi_i)} and working response
#'    \eqn{z_i^\top\gamma + (\tau_i-\pi_i)/\{\pi_i(1-\pi_i)\}}.
#' 3. **M-step, count part.** Maximize
#'    \eqn{\sum_i a_i(1-\tau_i)\log f_{NB}(y_i;\mu_i,\theta)}: the NB fit of
#'    [nbreg()] with prior weights \eqn{a_i(1-\tau_i)}.
#'
#' until the observed-data log-likelihood changes by less than
#' `control$tol`. EM never decreases the log-likelihood; the path is stored
#' in `loglik_path`, and a warning is given if it ever decreases.
#'
#' Standard errors come from the inverse of the observed information,
#' \eqn{-\partial^2\ell/\partial\psi\partial\psi^\top}, with
#' \eqn{\psi = (\beta, \gamma, \log\theta)}, computed from the analytic
#' second derivatives of the observed-data log-likelihood at the MLE.
#' \eqn{\mathrm{SE}(\hat\theta) = \hat\theta\,\mathrm{SE}(\log\hat\theta)}
#' (delta method).
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before the `|`, zero-model covariates after it. Without
#'   `|`, the same covariates are used in both parts. Offsets (count part
#'   only) may be given with `offset()` terms before the `|`.
#' @param data An optional data frame containing the variables in the model.
#' @param offset An optional offset for the count model (e.g. log
#'   exposure), added to any `offset()` terms in the formula.
#' @param control A list of control parameters; see [zinbreg_control()].
#' @param weights Optional non-negative prior weights \eqn{a_i}.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values (default: the
#'   `na.action` option, usually [stats::na.omit()]).
#'
#' @return An object of class `"zinbreg"`, a list with components
#' \describe{
#'   \item{coefficients}{List with elements `count` (\eqn{\beta}) and
#'     `zero` (\eqn{\gamma}).}
#'   \item{vcov}{Covariance matrix of
#'     \eqn{(\hat\beta, \hat\gamma, \log\hat\theta)}; [vcov.zinbreg()]
#'     returns the block for the coefficients.}
#'   \item{theta, SE.theta, log_theta, SE.log_theta}{\eqn{\hat\theta},
#'     \eqn{\log\hat\theta}, and their standard errors.}
#'   \item{loglik, loglik_path}{Maximized observed-data log-likelihood, and
#'     its value at the start and after each EM iteration.}
#'   \item{df}{Number of estimated parameters (\eqn{p + q + 1}).}
#'   \item{iter, converged}{Number of EM iterations and convergence flag.}
#'   \item{fitted.values}{\eqn{(1-\hat\pi_i)\hat\mu_i}, the fitted mean.}
#'   \item{mu, pi, tau}{\eqn{\hat\mu_i}, \eqn{\hat\pi_i}, and the final
#'     posterior probabilities \eqn{\tau_i}.}
#'   \item{residuals}{Response residuals \eqn{y_i - (1-\hat\pi_i)\hat\mu_i}.}
#'   \item{y, offset, prior.weights, nobs, df.residual}{Response, count
#'     offset, prior weights, number of observations with positive weight,
#'     and \eqn{n - p - q - 1}.}
#'   \item{call, formula, terms, model, xlevels, contrasts, na.action,
#'     control}{Information about the fit; `terms`, `xlevels` and
#'     `contrasts` are lists with `count` and `zero` elements.}
#' }
#'
#' @references
#' Lambert, D. (1992). Zero-inflated Poisson regression, with an application
#' to defects in manufacturing. *Technometrics*, 34(1), 1-14.
#'
#' Zeileis, A., Kleiber, C. and Jackman, S. (2008). Regression models for
#' count data in R. *Journal of Statistical Software*, 27(8), 1-25.
#'
#' @seealso [zinbreg_control()], [nbreg()]
#' @examples
#' set.seed(2)
#' n <- 500
#' x <- rnorm(n)
#' z <- rnorm(n)
#' structural_zero <- rbinom(n, 1, plogis(-1 + 1 * z))
#' y <- ifelse(structural_zero == 1, 0,
#'             rnbinom(n, mu = exp(1 + 0.5 * x), size = 2))
#' fit <- zinbreg(y ~ x | z)
#' fit
#' summary(fit)
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control(),
                    weights, subset, na.action) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  if (!is.list(control)) {
    stop("'control' must be a list, usually from zinbreg_control().",
         call. = FALSE)
  }
  control <- do.call(zinbreg_control, control)

  # Split the two-part formula y ~ count | zero.
  parts <- zinb_formulas(formula)

  # One model frame holding every variable of both parts, so that `subset`,
  # `weights`, `offset` and `na.action` apply to both parts in the same way.
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- parts$full
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Terms and model matrices for each part.
  mtX <- stats::terms(parts$count, data = data)
  mtZ <- stats::terms(parts$zero, data = data)
  if (!is.null(attr(mtZ, "offset"))) {
    stop("Offsets in the zero-inflation part are not supported.",
         call. = FALSE)
  }
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)

  y <- check_count_response(stats::model.response(mf, "numeric"))
  n <- length(y)
  a <- check_weights(stats::model.weights(mf), n)
  # The zero part has no offset, so model.offset() of the joint frame is the
  # count offset: offset() terms before the | plus the `offset` argument.
  o <- check_offset(stats::model.offset(mf), n)

  for (M in list(count = X, zero = Z)) {
    if (ncol(M) == 0L) {
      stop("Both parts of the model need at least one coefficient.",
           call. = FALSE)
    }
  }
  if (qr(X)$rank < ncol(X)) {
    stop("The count-model matrix is rank deficient. Remove redundant ",
         "covariates.", call. = FALSE)
  }
  if (qr(Z)$rank < ncol(Z)) {
    stop("The zero-model matrix is rank deficient. Remove redundant ",
         "covariates.", call. = FALSE)
  }
  if (sum(a > 0) <= ncol(X) + ncol(Z) + 1L) {
    stop("Not enough observations with positive weight to fit the model.",
         call. = FALSE)
  }
  if (!any(y[a > 0] == 0)) {
    stop("The response has no zeros; a zero-inflated model cannot be ",
         "fitted.", call. = FALSE)
  }
  if (all(y[a > 0] == 0)) {
    stop("All responses are zero; the model cannot be fitted.",
         call. = FALSE)
  }

  fit <- zinb_fit(X, Z, y, o, a, control)

  if (!fit$converged) {
    warning(sprintf(paste0("zinbreg: EM did not converge in %d iterations; ",
                           "increase 'maxit' in zinbreg_control()."),
                    control$maxit), call. = FALSE)
  }
  if (fit$decreased) {
    warning("zinbreg: the log-likelihood decreased during EM; check the ",
            "results.", call. = FALSE)
  }
  if (fit$theta > control$theta_warn) {
    warning(sprintf(paste0("zinbreg: theta is very large (%.3g); the count ",
                           "part is essentially Poisson."), fit$theta),
            call. = FALSE)
  }
  if (anyNA(fit$vcov)) {
    warning("zinbreg: the observed information matrix is not invertible; ",
            "standard errors are not available.", call. = FALSE)
  }

  # Names.
  names(fit$coefficients$count) <- colnames(X)
  names(fit$coefficients$zero) <- colnames(Z)
  nm <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
          "log(theta)")
  dimnames(fit$vcov) <- list(nm, nm)
  rn <- rownames(mf)
  for (v in c("fitted.values", "mu", "pi", "tau", "residuals")) {
    names(fit[[v]]) <- rn
  }
  names(y) <- rn
  fit$decreased <- NULL

  nobs <- sum(a > 0)
  out <- c(fit, list(
    df = ncol(X) + ncol(Z) + 1L,
    y = y,
    offset = o,
    prior.weights = a,
    nobs = nobs,
    df.residual = nobs - ncol(X) - ncol(Z) - 1L,
    call = cl,
    formula = formula,
    terms = list(count = mtX, zero = mtZ, full = attr(mf, "terms")),
    model = mf,
    xlevels = list(count = stats::.getXlevels(mtX, mf),
                   zero = stats::.getXlevels(mtZ, mf)),
    contrasts = list(count = attr(X, "contrasts"),
                     zero = attr(Z, "contrasts")),
    na.action = attr(mf, "na.action"),
    control = control
  ))
  class(out) <- "zinbreg"
  out
}

# Split y ~ count | zero into a count formula, a one-sided zero formula, and
# a combined formula y ~ count + zero used to build one joint model frame.
# Without |, the same right-hand side is used for both parts.
zinb_formulas <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula such as y ~ x | z.",
         call. = FALSE)
  }
  env <- environment(formula)
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    count_rhs <- rhs[[2L]]
    zero_rhs <- rhs[[3L]]
  } else {
    count_rhs <- zero_rhs <- rhs
  }
  list(
    count = stats::as.formula(call("~", lhs, count_rhs), env = env),
    zero = stats::as.formula(call("~", zero_rhs), env = env),
    full = stats::as.formula(call("~", lhs, call("+", count_rhs, zero_rhs)),
                             env = env)
  )
}

# ---------------------------------------------------------------------------
# EM algorithm (Section 3.2)
# ---------------------------------------------------------------------------

zinb_fit <- function(X, Z, y, o, a, control) {
  # Control list for the count-part M-step, which reuses nb_fit() (the
  # engine of nbreg()).
  nb_control <- list(tol = control$mstep_tol, maxit = control$mstep_maxit,
                     irls_tol = control$irls_tol,
                     irls_maxit = control$irls_maxit,
                     nr_tol = control$nr_tol, nr_maxit = control$nr_maxit,
                     theta_max = control$theta_max,
                     theta_warn = control$theta_warn, trace = FALSE)

  # --- Starting values ----------------------------------------------------
  # (beta, theta) from the NB fit of Section 3.1.
  nb0 <- nb_fit(X, y, o, a, nb_control)
  beta <- nb0$coefficients
  theta <- nb0$theta
  # gamma from an IRLS logistic regression of I(y_i = 0) on z_i.
  gamma <- logit_irls(Z, as.numeric(y == 0), a, rep(0, ncol(Z)),
                      control)$gamma

  mu <- as.vector(exp(X %*% beta + o))
  zeta <- as.vector(Z %*% gamma)
  ll_old <- zinb_loglik(y, mu, theta, zeta, a)
  path <- ll_old
  converged <- FALSE
  decreased <- FALSE
  iter <- 0L

  for (iter in seq_len(control$maxit)) {
    # --- E-step ---------------------------------------------------------
    tau <- zinb_estep(y, mu, theta, zeta)

    # --- M-step, zero part: logistic IRLS with fractional response tau --
    gamma <- logit_irls(Z, tau, a, gamma, control)$gamma

    # --- M-step, count part: NB fit with prior weights a_i (1 - tau_i) --
    nb <- nb_fit(X, y, o, a * (1 - tau), nb_control,
                 beta_start = beta, theta_start = theta)
    beta <- nb$coefficients
    theta <- nb$theta

    # --- Observed-data log-likelihood and convergence --------------------
    mu <- as.vector(exp(X %*% beta + o))
    zeta <- as.vector(Z %*% gamma)
    ll <- zinb_loglik(y, mu, theta, zeta, a)
    path <- c(path, ll)
    if (control$trace) {
      cat(sprintf("zinbreg EM iter %4d: logLik = %.10f  theta = %.8g\n",
                  iter, ll, theta))
    }
    # EM must never decrease l (allow for rounding error).
    if (ll < ll_old - 1e-6) decreased <- TRUE
    if (abs(ll - ll_old) < control$tol) {
      converged <- TRUE
      break
    }
    ll_old <- ll
  }

  # --- Standard errors from the analytic observed information --------------
  H <- zinb_hessian(X, Z, y, mu, theta, zeta, a)
  vcov <- tryCatch(solve(-H), error = function(e) {
    matrix(NA_real_, nrow(H), ncol(H))
  })
  k <- nrow(H)
  se_log_theta <- sqrt(vcov[k, k])
  p_zero <- stats::plogis(zeta)

  list(
    coefficients = list(count = as.vector(beta), zero = as.vector(gamma)),
    vcov = vcov,
    theta = theta,
    SE.theta = theta * se_log_theta,  # delta method
    log_theta = log(theta),
    SE.log_theta = se_log_theta,
    loglik = zinb_loglik(y, mu, theta, zeta, a),
    loglik_path = path,
    iter = iter,
    converged = converged,
    decreased = decreased,
    fitted.values = (1 - p_zero) * mu,  # E(Y_i) = (1 - pi_i) mu_i
    mu = mu,
    pi = p_zero,
    tau = zinb_estep(y, mu, theta, zeta),
    residuals = y - (1 - p_zero) * mu
  )
}

# E-step: tau_i = pi_i / {pi_i + (1 - pi_i) f_NB(0; mu_i, theta)} if y_i = 0,
# and 0 if y_i > 0.  Dividing through by pi_i gives
#   tau_i = 1 / {1 + exp(-zeta_i + log f_NB(0))} = plogis(zeta_i - log f0_i),
# which is the same number but cannot overflow.
zinb_estep <- function(y, mu, theta, zeta) {
  log_f0 <- theta * (log(theta) - log(theta + mu))
  ifelse(y == 0, stats::plogis(zeta - log_f0), 0)
}

# Observed-data log-likelihood (Section 2.2):
#   l = sum_{y_i = 0} a_i log{pi_i + (1 - pi_i) f_NB(0; mu_i, theta)}
#     + sum_{y_i > 0} a_i [log(1 - pi_i) + log f_NB(y_i; mu_i, theta)]
zinb_loglik <- function(y, mu, theta, zeta, a) {
  keep <- a > 0
  sum(a[keep] * zinb_loglik_i(y, mu, theta, zeta)[keep])
}

zinb_loglik_i <- function(y, mu, theta, zeta) {
  # log(pi) and log(1 - pi) on the log scale, so they never underflow.
  log_pi <- stats::plogis(zeta, log.p = TRUE)
  log_1mpi <- stats::plogis(-zeta, log.p = TRUE)
  ll <- numeric(length(y))
  zero <- y == 0
  # y = 0: log{pi + (1 - pi) f0} = log{exp(u1) + exp(u2)} with
  # u1 = log(pi), u2 = log(1 - pi) + log f0, computed by log-sum-exp.
  u1 <- log_pi[zero]
  u2 <- log_1mpi[zero] + theta * (log(theta) - log(theta + mu[zero]))
  um <- pmax(u1, u2)
  ll[zero] <- um + log(exp(u1 - um) + exp(u2 - um))
  # y > 0: log(1 - pi) + log f_NB(y; mu, theta)
  ll[!zero] <- log_1mpi[!zero] + nb_logf(y[!zero], mu[!zero], theta)
  ll
}

# M-step, zero part: logistic regression with fractional response tau by
# IRLS (Section 3.2):
#   w_i = a_i pi_i (1 - pi_i)
#   z_i = zeta_i + (tau_i - pi_i) / {pi_i (1 - pi_i)}
#   gamma^(t+1) = (Z' W Z)^{-1} Z' W z
# Note W z = a_i pi_i (1 - pi_i) zeta_i + a_i (tau_i - pi_i), so we compute
# Z' W z in that form; it is algebraically identical and avoids dividing by
# pi_i (1 - pi_i) when pi_i is close to 0 or 1.
logit_irls <- function(Z, tau, a, gamma, control) {
  # Weighted Bernoulli log-likelihood sum a_i {tau_i log pi_i
  #   + (1 - tau_i) log(1 - pi_i)}, the objective of this M-step.
  obj <- function(zeta) {
    sum(a * (tau * stats::plogis(zeta, log.p = TRUE) +
               (1 - tau) * stats::plogis(-zeta, log.p = TRUE)))
  }
  zeta <- as.vector(Z %*% gamma)
  ll_old <- obj(zeta)
  converged <- FALSE

  for (it in seq_len(control$irls_maxit)) {
    p <- stats::plogis(zeta)
    w <- a * p * (1 - p)
    Wz <- w * zeta + a * (tau - p)
    gamma_new <- as.vector(solve(crossprod(Z, w * Z), crossprod(Z, Wz)))
    zeta_new <- as.vector(Z %*% gamma_new)
    ll_new <- obj(zeta_new)

    # Safeguard (step halving), as in nb_irls(): Fisher scoring for the
    # logistic model normally increases the objective; if not, halve.
    halve <- 0L
    while ((!is.finite(ll_new) ||
            ll_new < ll_old - 1e-12 * (abs(ll_old) + 0.1)) && halve < 30L) {
      gamma_new <- (gamma + gamma_new) / 2
      zeta_new <- as.vector(Z %*% gamma_new)
      ll_new <- obj(zeta_new)
      halve <- halve + 1L
    }

    done <- abs(ll_new - ll_old) < control$irls_tol * (abs(ll_new) + 0.1)
    gamma <- gamma_new
    zeta <- zeta_new
    ll_old <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }
  list(gamma = gamma, iter = it, converged = converged)
}

# ---------------------------------------------------------------------------
# Analytic Hessian of the observed-data log-likelihood
# ---------------------------------------------------------------------------
# Parameters psi = (beta, gamma, phi) with phi = log(theta).  Each l_i depends
# on beta only through eta_i = x_i' beta + o_i, on gamma only through
# zeta_i = z_i' gamma, so by the chain rule
#   d2 l / d beta d beta' = sum a_i h_ee,i x_i x_i'
#   d2 l / d beta d gamma' = sum a_i h_ez,i x_i z_i'
#   d2 l / d beta d phi    = sum a_i h_ep,i x_i        (and so on),
# where h_ee, h_ez, ... are second derivatives of l_i in (eta, zeta, phi).
#
# y_i > 0:  l_i = log(1 - pi) + log f_NB(y; mu, theta)
#   h_ee = -theta mu (theta + y) / (theta + mu)^2
#   h_zz = -pi (1 - pi),         h_ez = h_zp = 0
#   h_ep = theta mu (y - mu) / (theta + mu)^2
#   h_pp = theta^2 H_i(theta) + theta U_i(theta)   (U, H from Section 3.1)
#
# y_i = 0:  l_i = log p0, p0 = pi + (1 - pi) f0, f0 = exp(g),
#   g = theta {log(theta) - log(theta + mu)}.
#   Derivatives of g:
#     g_e = -theta mu / (theta + mu),        g_ee = -theta^2 mu / (theta+mu)^2
#     g_p = theta {log(theta/(theta+mu)) + mu/(theta+mu)}
#     g_pp = mu^2 theta / (theta+mu)^2 + g_p
#     g_ep = -theta mu^2 / (theta + mu)^2
#   Derivatives of p0 (pi' = pi (1 - pi), f0' = f0 g'):
#     p0_e = (1-pi) f0 g_e,   p0_p = (1-pi) f0 g_p,
#     p0_z = pi (1-pi) (1 - f0)
#     p0_ee = (1-pi) f0 (g_ee + g_e^2),  p0_pp = (1-pi) f0 (g_pp + g_p^2)
#     p0_ep = (1-pi) f0 (g_ep + g_e g_p)
#     p0_zz = pi (1-pi) (1 - 2 pi) (1 - f0)
#     p0_ez = -pi (1-pi) f0 g_e,         p0_zp = -pi (1-pi) f0 g_p
#   and h_jk = p0_jk / p0 - p0_j p0_k / p0^2.
zinb_hessian <- function(X, Z, y, mu, theta, zeta, a) {
  p_zero <- stats::plogis(zeta)
  n <- length(y)
  h_ee <- h_zz <- h_pp <- h_ez <- h_ep <- h_zp <- numeric(n)

  # y > 0
  pos <- y > 0
  yp <- y[pos]
  mp <- mu[pos]
  pp <- p_zero[pos]
  U <- digamma(yp + theta) - digamma(theta) + log(theta) + 1 -
    log(theta + mp) - (yp + theta) / (theta + mp)
  H <- trigamma(yp + theta) - trigamma(theta) + 1 / theta -
    2 / (theta + mp) + (yp + theta) / (theta + mp)^2
  h_ee[pos] <- -theta * mp * (theta + yp) / (theta + mp)^2
  h_zz[pos] <- -pp * (1 - pp)
  h_ep[pos] <- theta * mp * (yp - mp) / (theta + mp)^2
  h_pp[pos] <- theta^2 * H + theta * U

  # y = 0
  z0 <- !pos
  m0 <- mu[z0]
  p0pi <- p_zero[z0]
  g <- theta * (log(theta) - log(theta + m0))
  f0 <- exp(g)
  p0 <- p0pi + (1 - p0pi) * f0
  g_e <- -theta * m0 / (theta + m0)
  g_ee <- -theta^2 * m0 / (theta + m0)^2
  g_p <- theta * (log(theta / (theta + m0)) + m0 / (theta + m0))
  g_pp <- theta * m0^2 / (theta + m0)^2 + g_p
  g_ep <- -theta * m0^2 / (theta + m0)^2

  d_e <- (1 - p0pi) * f0 * g_e
  d_p <- (1 - p0pi) * f0 * g_p
  d_z <- p0pi * (1 - p0pi) * (1 - f0)
  dd_ee <- (1 - p0pi) * f0 * (g_ee + g_e^2)
  dd_pp <- (1 - p0pi) * f0 * (g_pp + g_p^2)
  dd_ep <- (1 - p0pi) * f0 * (g_ep + g_e * g_p)
  dd_zz <- p0pi * (1 - p0pi) * (1 - 2 * p0pi) * (1 - f0)
  dd_ez <- -p0pi * (1 - p0pi) * f0 * g_e
  dd_zp <- -p0pi * (1 - p0pi) * f0 * g_p

  h_ee[z0] <- dd_ee / p0 - d_e^2 / p0^2
  h_pp[z0] <- dd_pp / p0 - d_p^2 / p0^2
  h_ep[z0] <- dd_ep / p0 - d_e * d_p / p0^2
  h_zz[z0] <- dd_zz / p0 - d_z^2 / p0^2
  h_ez[z0] <- dd_ez / p0 - d_e * d_z / p0^2
  h_zp[z0] <- dd_zp / p0 - d_z * d_p / p0^2

  # Assemble the (p + q + 1) x (p + q + 1) Hessian.
  H_bb <- crossprod(X, (a * h_ee) * X)
  H_bg <- crossprod(X, (a * h_ez) * Z)
  H_gg <- crossprod(Z, (a * h_zz) * Z)
  H_bp <- crossprod(X, a * h_ep)
  H_gp <- crossprod(Z, a * h_zp)
  H_pp <- sum(a * h_pp)
  rbind(cbind(H_bb, H_bg, H_bp),
        cbind(t(H_bg), H_gg, H_gp),
        cbind(t(H_bp), t(H_gp), H_pp))
}
