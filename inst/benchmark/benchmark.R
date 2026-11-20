# Benchmark: nbzinbWebster vs. the reference implementations
#
#   nbreg()   vs. MASS::glm.nb()
#   zinbreg() vs. pscl::zeroinfl(dist = "negbin")
#
# on pscl::bioChemists and AER::NMES1988 (project Section 5).
#
# For each data set and model this script reports
#   * the estimates and standard errors side by side,
#   * the maximum absolute differences, checked against Table 4,
#   * computation times (median elapsed seconds over repeated fits).
#
# Run from the package root, with the package installed:
#   Rscript inst/benchmark/benchmark.R
# Results are written to inst/benchmark/results/.

suppressPackageStartupMessages({
  library(nbzinbWebster)
  library(MASS)
  library(pscl)
})

if (!dir.exists(file.path("inst", "benchmark"))) {
  stop("Run this script from the package root directory.")
}
out_dir <- file.path("inst", "benchmark", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(651)

# ---------------------------------------------------------------------------
# Data and models
# ---------------------------------------------------------------------------

data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

datasets <- list(
  bioChemists = list(
    data = bioChemists,
    nb = art ~ fem + mar + kid5 + phd + ment,
    zinb = art ~ fem + mar + kid5 + phd + ment |
      fem + mar + kid5 + phd + ment
  ),
  NMES1988 = list(
    data = NMES1988,
    nb = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school +
      insurance | hospital + chronic + gender + school + insurance
  )
)

# Required agreement (Table 4).
tol_nb <- c(coef = 1e-5, se = 1e-5, theta = 1e-4, loglik = 1e-4)
tol_zinb <- c(coef = 1e-3, se = 1e-3, theta = 1e-3, loglik = 1e-3)

# Number of repeated fits used for timing.
reps_nb <- 10
reps_zinb <- 3

# Reference functions with tight tolerances (Section 5 tips).
ref_nb <- function(f, d) {
  glm.nb(f, data = d, control = glm.control(epsilon = 1e-10, maxit = 100))
}
ref_zinb <- function(f, d) {
  zeroinfl(f, data = d, dist = "negbin",
           control = zeroinfl.control(reltol = 1e-12))
}

# Median elapsed time (seconds) of `reps` evaluations of fun().
time_it <- function(fun, reps) {
  median(vapply(seq_len(reps),
                function(i) system.time(fun())[["elapsed"]], numeric(1)))
}

max_abs <- function(a, b) max(abs(unname(a) - unname(b)))

# ---------------------------------------------------------------------------
# NB: nbreg() vs glm.nb()
# ---------------------------------------------------------------------------

bench_nb <- function(name, spec) {
  f <- spec$nb
  d <- spec$data
  mine <- nbreg(f, data = d)
  ref <- ref_nb(f, d)

  se_mine <- sqrt(diag(vcov(mine)))
  se_ref <- summary(ref)$coefficients[, "Std. Error"]

  coef_table <- data.frame(
    term = c(names(coef(mine)), "theta"),
    estimate_nbreg = c(coef(mine), mine$theta),
    estimate_glm.nb = c(coef(ref), ref$theta),
    se_nbreg = c(se_mine, mine$SE.theta),
    se_glm.nb = c(se_ref, ref$SE.theta),
    row.names = NULL, check.names = FALSE
  )
  coef_table$abs_diff_estimate <-
    abs(coef_table$estimate_nbreg - coef_table$estimate_glm.nb)
  coef_table$abs_diff_se <- abs(coef_table$se_nbreg - coef_table$se_glm.nb)

  diffs <- c(coef = max_abs(coef(mine), coef(ref)),
             se = max_abs(se_mine, se_ref),
             theta = abs(mine$theta - ref$theta),
             loglik = abs(mine$loglik - as.numeric(logLik(ref))))
  agreement <- data.frame(
    dataset = name, model = "NB",
    quantity = c("Coefficients", "Standard errors", "theta",
                 "Log-likelihood"),
    max_abs_diff = unname(diffs),
    tolerance = unname(tol_nb),
    pass = unname(diffs < tol_nb)
  )

  timing <- data.frame(
    dataset = name, model = "NB",
    ours = "nbreg", reference = "MASS::glm.nb",
    seconds_ours = time_it(function() nbreg(f, data = d), reps_nb),
    seconds_reference = time_it(function() ref_nb(f, d), reps_nb),
    reps = reps_nb,
    iterations = mine$iter, converged = mine$converged,
    loglik_ours = mine$loglik,
    loglik_reference = as.numeric(logLik(ref))
  )

  list(coef_table = coef_table, agreement = agreement, timing = timing)
}

# ---------------------------------------------------------------------------
# ZINB: zinbreg() vs zeroinfl()
# ---------------------------------------------------------------------------

bench_zinb <- function(name, spec) {
  f <- spec$zinb
  d <- spec$data
  mine <- zinbreg(f, data = d)
  ref <- ref_zinb(f, d)

  se_mine <- sqrt(diag(vcov(mine)))
  se_ref <- sqrt(diag(vcov(ref)))
  stopifnot(identical(names(coef(mine)), names(coef(ref))))

  # log(theta) is reported with its SE, as zeroinfl() does.
  coef_table <- data.frame(
    term = c(names(coef(mine)), "log(theta)"),
    estimate_zinbreg = c(coef(mine), mine$log_theta),
    estimate_zeroinfl = c(coef(ref), log(ref$theta)),
    se_zinbreg = c(se_mine, mine$SE.log_theta),
    se_zeroinfl = c(se_ref, ref$SE.logtheta),
    row.names = NULL, check.names = FALSE
  )
  coef_table$abs_diff_estimate <-
    abs(coef_table$estimate_zinbreg - coef_table$estimate_zeroinfl)
  coef_table$abs_diff_se <- abs(coef_table$se_zinbreg - coef_table$se_zeroinfl)

  diffs <- c(coef = max_abs(coef(mine), coef(ref)),
             se = max_abs(se_mine, se_ref),
             theta = abs(mine$log_theta - log(ref$theta)),
             loglik = abs(mine$loglik - as.numeric(logLik(ref))))
  agreement <- data.frame(
    dataset = name, model = "ZINB",
    quantity = c("Coefficients", "Standard errors", "log(theta)",
                 "Log-likelihood"),
    max_abs_diff = unname(diffs),
    tolerance = unname(tol_zinb),
    pass = unname(diffs < tol_zinb)
  )

  timing <- data.frame(
    dataset = name, model = "ZINB",
    ours = "zinbreg", reference = "pscl::zeroinfl",
    seconds_ours = time_it(function() zinbreg(f, data = d), reps_zinb),
    seconds_reference = time_it(function() ref_zinb(f, d), reps_zinb),
    reps = reps_zinb,
    iterations = mine$iter, converged = mine$converged,
    loglik_ours = mine$loglik,
    loglik_reference = as.numeric(logLik(ref))
  )

  list(coef_table = coef_table, agreement = agreement, timing = timing)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

results <- list()
for (name in names(datasets)) {
  message("Benchmarking ", name, " ...")
  results[[name]] <- list(
    NB = bench_nb(name, datasets[[name]]),
    ZINB = bench_zinb(name, datasets[[name]])
  )
}

agreement <- do.call(rbind, lapply(results, function(r) {
  rbind(r$NB$agreement, r$ZINB$agreement)
}))
timing <- do.call(rbind, lapply(results, function(r) {
  rbind(r$NB$timing, r$ZINB$timing)
}))
rownames(agreement) <- rownames(timing) <- NULL

# ---------------------------------------------------------------------------
# Print
# ---------------------------------------------------------------------------

op <- options(width = 120, digits = 6)
for (name in names(results)) {
  for (model in c("NB", "ZINB")) {
    cat("\n====", name, "-", model, "====\n")
    print(results[[name]][[model]]$coef_table, row.names = FALSE)
  }
}
cat("\n==== Agreement (maximum absolute differences vs. Table 4) ====\n")
print(transform(agreement, max_abs_diff = signif(max_abs_diff, 3)),
      row.names = FALSE)
cat("\n==== Computation time (median elapsed seconds) ====\n")
print(timing[, c("dataset", "model", "ours", "reference", "seconds_ours",
                 "seconds_reference", "reps", "iterations", "converged")],
      row.names = FALSE)
options(op)

if (!all(agreement$pass)) {
  warning("Some quantities are outside the Table 4 tolerances.")
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

saveRDS(list(results = results, agreement = agreement, timing = timing,
             date = Sys.time(),
             session = utils::sessionInfo()),
        file.path(out_dir, "benchmark_results.rds"))
utils::write.csv(agreement, file.path(out_dir, "agreement.csv"),
                 row.names = FALSE)
utils::write.csv(timing, file.path(out_dir, "timing.csv"), row.names = FALSE)
for (name in names(results)) {
  for (model in c("NB", "ZINB")) {
    utils::write.csv(results[[name]][[model]]$coef_table,
                     file.path(out_dir,
                               sprintf("coef_%s_%s.csv", name, model)),
                     row.names = FALSE)
  }
}
writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
message("Results saved to ", normalizePath(out_dir))
