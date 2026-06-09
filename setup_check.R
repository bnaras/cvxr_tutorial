#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## useR! 2026 CVXR tutorial -- environment setup check
##
## Run this BEFORE the tutorial to confirm your machine is ready:
##     Rscript setup_check.R      (terminal)
##     source("setup_check.R")    (R console / RStudio)
##
## One-line install for everything (run once if anything is missing):
##     install.packages(c("CVXR", "scip", "Uno", "sparsediff"))
##
## Base R only -- no external dependencies -- so it runs even if a package
## failed to install. CORE = what every hands-on block needs; DEMO = the two
## instructor demos (MI-SOCP via scip, DNLP via Uno + sparsediff).
## ---------------------------------------------------------------------------

mark <- function(ok) if (isTRUE(ok)) "PASS" else if (is.na(ok)) "WARN" else "FAIL"
line <- function(tag, label, detail = "") cat(sprintf("  [%s] %-28s %s\n", tag, label, detail))
hdr  <- function(x) cat(sprintf("\n== %s ==\n", x))

core_ok <- TRUE   # everything a hands-on block needs
demo_ok <- TRUE   # the two demos

## --- 1. R version --------------------------------------------------------
hdr("R version")
rv_ok <- getRversion() >= "4.3.0"
line(mark(rv_ok), "R >= 4.3.0", R.version.string)
if (!rv_ok) core_ok <- FALSE

## --- 2. Packages present (+ version floors) ------------------------------
hdr("Packages")
check_pkg <- function(pkg, min = NULL, tier = "core") {
  v <- tryCatch(packageVersion(pkg), error = function(e) NULL)
  if (is.null(v)) { line("FAIL", pkg, "not installed"); return(FALSE) }
  ok <- is.null(min) || v >= min
  line(mark(ok), pkg, sprintf("%s%s", v,
       if (!ok) sprintf("  (need >= %s)", min) else ""))
  ok
}
## CVXR + the four solvers it imports (auto-installed with CVXR)
core_ok <- check_pkg("CVXR",     "1.9.1") && core_ok
core_ok <- check_pkg("clarabel")          && core_ok
core_ok <- check_pkg("osqp")              && core_ok
core_ok <- check_pkg("scs")               && core_ok
core_ok <- check_pkg("highs",    "1.14")  && core_ok   # CVXR requires highs >= 1.14
## Demo-only extras (CRAN, but in CVXR Enhances -> separate installs)
demo_ok <- check_pkg("scip",       "1.10", "demo") && demo_ok   # MI-SOCP demo
demo_ok <- check_pkg("Uno",        "2.7.3","demo") && demo_ok   # DNLP backend (pulls rmumps)
demo_ok <- check_pkg("sparsediff", "0.4.0","demo") && demo_ok   # DNLP autodiff

## --- 3. Functional solve test: the four core solvers ---------------------
## Tiny LP: minimize sum(x) s.t. x >= 1, x in R^2  ->  optimal value 2.
hdr("Core solvers (tiny LP, expect value 2)")
if (requireNamespace("CVXR", quietly = TRUE)) {
  suppressMessages(library(CVXR))
  x    <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  for (s in c("CLARABEL", "OSQP", "SCS", "HIGHS")) {
    res <- tryCatch({
      val <- psolve(prob, solver = s)
      list(ok = is.numeric(val) && abs(val - 2) < 1e-4 &&
                status(prob) %in% c("optimal", "optimal_inaccurate"),
           detail = sprintf("value=%.4f status=%s", val, status(prob)))
    }, error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
    line(mark(res$ok), s, res$detail)
    if (!isTRUE(res$ok)) core_ok <- FALSE
  }
} else {
  line("FAIL", "CVXR", "cannot load -- skipping solve tests")
  core_ok <- FALSE
}

## --- 4. Demo: SCIP integer solve (MI-SOCP-capable interface) -------------
## Tiny integer LP: minimize sum(x) s.t. x >= 1.5, integer -> value 4.
hdr("MI-SOCP demo solver (SCIP integer solve, expect value 4)")
if (requireNamespace("CVXR", quietly = TRUE) &&
    requireNamespace("scip", quietly = TRUE)) {
  xi <- Variable(2, integer = TRUE)
  pi <- Problem(Minimize(sum(xi)), list(xi >= 1.5))
  res <- tryCatch({
    val <- psolve(pi, solver = "SCIP")
    list(ok = is.numeric(val) && abs(val - 4) < 1e-4, detail = sprintf("value=%.4f", val))
  }, error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
  line(mark(res$ok), "SCIP", res$detail)
  if (!isTRUE(res$ok)) demo_ok <- FALSE
} else {
  line("WARN", "SCIP", "scip not installed -- MI-SOCP demo only")
  demo_ok <- FALSE
}

## --- 5. Demo: DNLP backend present ---------------------------------------
hdr("DNLP demo backends")
dnlp_ok <- requireNamespace("Uno", quietly = TRUE) &&
           requireNamespace("sparsediff", quietly = TRUE)
line(if (dnlp_ok) "PASS" else "WARN", "Uno + sparsediff",
     if (dnlp_ok) "loadable" else "missing -- DNLP demo only (ipopt NOT needed)")
if (!dnlp_ok) demo_ok <- FALSE

## --- Summary -------------------------------------------------------------
hdr("Summary")
cat(sprintf("  CORE  (every hands-on block): %s\n", if (core_ok) "READY" else "NOT READY"))
cat(sprintf("  DEMOS (MI-SOCP + DNLP)     : %s\n", if (demo_ok) "READY" else "incomplete (demos only)"))
if (!core_ok) {
  cat("\n  -> Fix CORE first. Re-run:\n")
  cat("     install.packages(c(\"CVXR\", \"scip\", \"Uno\", \"sparsediff\"))\n")
  cat("     (or just install.packages(\"CVXR\") for the hands-on blocks)\n")
} else if (!demo_ok) {
  cat("\n  -> You can do every hands-on block. The two demos are run by the\n")
  cat("     instructor; install the extras only if you want to reproduce them:\n")
  cat("     install.packages(c(\"scip\", \"Uno\", \"sparsediff\"))\n")
} else {
  cat("\n  All set -- core and both demos are ready.\n")
}
cat("\n")
