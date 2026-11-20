# Simulation study: does zinbreg() recover known parameters?  (Section 5.1)
#
# Design
#   R = 100 data sets, each with n = 500 observations, from the ZINB model
#     count part:  log mu_i     = beta0 + beta1 x_i,   x_i ~ N(0, 1)
#     zero part:   logit pi_i   = gamma0 + gamma1 z_i, z_i ~ N(0, 1)
#     dispersion:  theta
#   y_i = 0 with probability pi_i (structural zero), otherwise
#   y_i ~ NB(mu_i, theta).  The covariates are drawn afresh for every data
#   set, so the results average over the design.
#
# For every parameter we report
#   * bias          = mean(estimate) - true value,
#   * MC SD         = standard deviation of the estimates over replicates,
#   * mean SE       = average of the estimated standard errors,
#   * SE / MC SD    = ratio (close to 1 if the SEs are right),
#   * coverage      = proportion of 95% Wald intervals (est +/- 1.96 SE)
#                     that contain the true value.
# theta is reported on the log scale (the scale of its Wald interval, as in
# zeroinfl) and on the original scale (SE by the delta method).
#
# Run from the package root, with the package installed:
#   Rscript inst/benchmark/simulation.R
# Results are written to inst/benchmark/results/simulation_results.rds.

suppressPackageStartupMessages(library(nbzinbWebster))

if (!dir.exists(file.path("inst", "benchmark"))) {
  stop("Run this script from the package root directory.")
}
out_dir <- file.path("inst", "benchmark", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

seed <- 651
R <- 100
n <- 500
true <- c(beta0 = 1, beta1 = 0.5, gamma0 = -1, gamma1 = 1, theta = 2)
level <- 0.95
zcrit <- qnorm(1 - (1 - level) / 2)

set.seed(seed)

# ---------------------------------------------------------------------------
# Data generation
# ---------------------------------------------------------------------------

simulate_zinb <- function(n, true) {
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(true[["beta0"]] + true[["beta1"]] * x)
  pi0 <- plogis(true[["gamma0"]] + true[["gamma1"]] * z)
  structural_zero <- rbinom(n, 1, pi0)
  y <- ifelse(structural_zero == 1, 0,
              rnbinom(n, size = true[["theta"]], mu = mu))
  data.frame(y = y, x = x, z = z)
}

# ---------------------------------------------------------------------------
# Fit one replicate; record estimates, SEs, convergence and warnings
# ---------------------------------------------------------------------------

fit_one <- function(d) {
  warn <- character(0)
  fit <- withCallingHandlers(
    tryCatch(zinbreg(y ~ x | z, data = d), error = function(e) e),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, error = conditionMessage(fit), warnings = warn))
  }
  se <- sqrt(diag(fit$vcov))  # (beta0, beta1, gamma0, gamma1, log theta)
  list(
    ok = TRUE,
    est = c(beta0 = coef(fit)[[1]], beta1 = coef(fit)[[2]],
            gamma0 = coef(fit)[[3]], gamma1 = coef(fit)[[4]],
            log_theta = fit$log_theta, theta = fit$theta),
    se = c(beta0 = se[[1]], beta1 = se[[2]], gamma0 = se[[3]],
           gamma1 = se[[4]], log_theta = fit$SE.log_theta,
           theta = fit$SE.theta),
    converged = fit$converged,
    iter = fit$iter,
    prop_zero = mean(d$y == 0),
    warnings = warn
  )
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

message(sprintf("Simulating %d data sets of size %d (seed %d) ...",
                R, n, seed))
elapsed <- system.time({
  reps <- vector("list", R)
  for (r in seq_len(R)) {
    reps[[r]] <- fit_one(simulate_zinb(n, true))
    if (r %% 10 == 0) message("  ", r, " / ", R)
  }
})[["elapsed"]]

ok <- vapply(reps, `[[`, logical(1), "ok")
if (any(!ok)) {
  message(sum(!ok), " replicate(s) failed and are excluded:")
  print(unique(unlist(lapply(reps[!ok], `[[`, "error"))))
}
est <- do.call(rbind, lapply(reps[ok], `[[`, "est"))
se <- do.call(rbind, lapply(reps[ok], `[[`, "se"))
converged <- vapply(reps[ok], `[[`, logical(1), "converged")
iter <- vapply(reps[ok], `[[`, numeric(1), "iter")
prop_zero <- vapply(reps[ok], `[[`, numeric(1), "prop_zero")
warnings_all <- unlist(lapply(reps, `[[`, "warnings"))

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

true_ext <- c(true[c("beta0", "beta1", "gamma0", "gamma1")],
              log_theta = log(true[["theta"]]), theta = true[["theta"]])
params <- names(true_ext)

summary_table <- data.frame(
  parameter = params,
  true = unname(true_ext),
  mean_estimate = colMeans(est[, params]),
  bias = colMeans(est[, params]) - unname(true_ext),
  mc_sd = apply(est[, params], 2, sd),
  mean_se = colMeans(se[, params]),
  row.names = NULL
)
summary_table$se_ratio <- summary_table$mean_se / summary_table$mc_sd
# Relative bias in percent of the true value.
summary_table$rel_bias_pct <- 100 * summary_table$bias / summary_table$true
# Monte Carlo SE of the bias, to judge whether a bias is distinguishable
# from simulation noise.
summary_table$mc_se_bias <- summary_table$mc_sd / sqrt(nrow(est))
# Coverage of the 95% Wald interval est +/- 1.96 SE.
lower <- est[, params] - zcrit * se[, params]
upper <- est[, params] + zcrit * se[, params]
covered <- sweep(lower, 2, true_ext, "<=") & sweep(upper, 2, true_ext, ">=")
summary_table$coverage <- colMeans(covered)
summary_table <- summary_table[, c("parameter", "true", "mean_estimate",
                                   "bias", "mc_se_bias", "rel_bias_pct",
                                   "mc_sd", "mean_se", "se_ratio",
                                   "coverage")]

# With R replicates the Monte Carlo SE of an empirical coverage near 0.95 is
# sqrt(0.95 * 0.05 / R), about 0.022 for R = 100.
coverage_mc_se <- sqrt(level * (1 - level) / nrow(est))

# ---------------------------------------------------------------------------
# Print
# ---------------------------------------------------------------------------

op <- options(width = 120)
cat("\n==== ZINB simulation study ====\n")
cat(sprintf("R = %d data sets, n = %d, seed = %d, model y ~ x | z\n",
            R, n, seed))
cat(sprintf("Successful fits: %d / %d; converged: %d; ",
            sum(ok), R, sum(converged)))
cat(sprintf("EM iterations: median %g (range %g-%g)\n",
            median(iter), min(iter), max(iter)))
cat(sprintf("Proportion of zeros per data set: mean %.3f (range %.3f-%.3f)\n",
            mean(prop_zero), min(prop_zero), max(prop_zero)))
cat(sprintf("Total time: %.1f s (%.2f s per fit)\n\n", elapsed, elapsed / R))
print(transform(summary_table,
                mean_estimate = round(mean_estimate, 4),
                bias = round(bias, 4), mc_se_bias = round(mc_se_bias, 4),
                rel_bias_pct = round(rel_bias_pct, 2),
                mc_sd = round(mc_sd, 4), mean_se = round(mean_se, 4),
                se_ratio = round(se_ratio, 3)),
      row.names = FALSE)
cat(sprintf(paste0("\nMonte Carlo SE of a coverage near %.2f: %.3f ",
                   "(so roughly %.2f-%.2f is consistent with %.2f)\n"),
            level, coverage_mc_se, level - 2 * coverage_mc_se,
            min(1, level + 2 * coverage_mc_se), level))
if (length(warnings_all) > 0) {
  cat("\nWarnings during fitting:\n")
  print(table(warnings_all))
}
options(op)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

saveRDS(list(
  summary = summary_table,
  estimates = est,
  std_errors = se,
  covered = covered,
  converged = converged,
  iterations = iter,
  prop_zero = prop_zero,
  n_failed = sum(!ok),
  warnings = warnings_all,
  coverage_mc_se = coverage_mc_se,
  settings = list(R = R, n = n, seed = seed, true = true, level = level,
                  formula = "y ~ x | z"),
  elapsed_seconds = elapsed,
  date = Sys.time(),
  session = utils::sessionInfo()
), file.path(out_dir, "simulation_results.rds"))
message("Results saved to ",
        normalizePath(file.path(out_dir, "simulation_results.rds")))
