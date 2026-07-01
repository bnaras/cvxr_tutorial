#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## useR! 2026 CVXR tutorial -- environment setup check
##
## Run this BEFORE the tutorial to confirm your machine is ready:
##     Rscript setup_check.R      (terminal)
##     source("setup_check.R")    (R console / RStudio)
##
## One-line install for everything (run once if anything is missing):
##     install.packages(c("CVXR", "scip", "Uno", "sparsediff"))   # solvers
##     install.packages(c("ggplot2","RColorBrewer","tidyr","nnls","glmnet","png","bench"))
##
## Base R only -- no external dependencies -- so it runs even if a package
## failed to install. CORE = the four built-in solvers used by most chapters;
## EXTRA = scip (ch.11, MI-SOCP) and Uno + sparsediff (ch.14, nonlinear) -- both
## are hands-on chapters, but if an EXTRA fails you can still do everything else.
## ---------------------------------------------------------------------------

mark <- function(ok) if (isTRUE(ok)) "PASS" else if (is.na(ok)) "WARN" else "FAIL"
line <- function(tag, label, detail = "") cat(sprintf("  [%s] %-28s %s\n", tag, label, detail))
hdr  <- function(x) cat(sprintf("\n== %s ==\n", x))

core_ok  <- TRUE   # CVXR + four built-in solvers + helpers
extra_ok <- TRUE   # scip (ch.11) and Uno + sparsediff (ch.14)

## --- 1. R version --------------------------------------------------------
hdr("R version")
rv_ok <- getRversion() >= "4.3.0"
line(mark(rv_ok), "R >= 4.3.0", R.version.string)
if (!rv_ok) core_ok <- FALSE

## --- 2. Packages: installed vs. needed vs. version used to build ---------
## REF = the versions this tutorial was developed and tested against. A
## mismatch is usually harmless; it is shown so that if something misbehaves
## during the session you can see at a glance whether your versions differ
## from the author's and factor that in.
REF <- c(CVXR = "1.9.1", clarabel = "0.11.2", osqp = "1.0.0", scs = "3.2.7",
         highs = "1.14.0.2", ggplot2 = "4.0.3", RColorBrewer = "1.1.3",
         tidyr = "1.3.2", nnls = "1.6", glmnet = "5.0", png = "0.1.9",
         bench = "1.1.4", scip = "1.10.0.3", Uno = "2.7.4",
         sparsediff = "0.4.0")

hdr("Packages  (installed = you have | needed = minimum | used = built with)")
cat(sprintf("  %-6s %-12s %-14s %-11s %s\n", "", "package", "installed", "needed", "used"))
check_pkg <- function(pkg, min = NULL) {
  v    <- tryCatch(packageVersion(pkg), error = function(e) NULL)
  ref  <- if (pkg %in% names(REF)) REF[[pkg]] else "--"
  need <- if (is.null(min)) "any" else paste(">=", min)
  if (is.null(v)) {
    cat(sprintf("  [%s] %-12s %-14s %-11s %s\n", "FAIL", pkg, "not installed", need, ref))
    return(FALSE)
  }
  ok <- is.null(min) || v >= min
  cat(sprintf("  [%s] %-12s %-14s %-11s %s\n", mark(ok), pkg, as.character(v), need, ref))
  ok
}
## CVXR + the four solvers it imports (auto-installed with CVXR)
core_ok <- check_pkg("CVXR",  "1.9.1") && core_ok
core_ok <- check_pkg("clarabel")       && core_ok
core_ok <- check_pkg("osqp")           && core_ok
core_ok <- check_pkg("scs")            && core_ok
core_ok <- check_pkg("highs", "1.14")  && core_ok   # CVXR requires highs >= 1.14
## Helper packages the example chapters library() (CRAN installs only; base and
## recommended packages that ship with R -- Matrix, MASS, survival -- are not
## listed here because they are always present).
for (h in c("ggplot2", "RColorBrewer", "tidyr", "nnls", "glmnet", "png", "bench"))
  core_ok <- check_pkg(h) && core_ok
## Extra solvers for chapters 11 and 14 (CRAN, in CVXR Enhances -> separate installs)
extra_ok <- check_pkg("scip",       "1.10")  && extra_ok   # ch.11 MI-SOCP
extra_ok <- check_pkg("Uno",        "2.7.4") && extra_ok   # ch.14 nonlinear (pulls rmumps); 2.7.4 carries the Windows MUMPS fix
extra_ok <- check_pkg("sparsediff", "0.4.0") && extra_ok   # ch.14 autodiff

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

## --- 4. Chapter 11: SCIP integer solve (MI-SOCP-capable interface) -------
## Tiny integer LP: minimize sum(x) s.t. x >= 1.5, integer -> value 4.
hdr("Chapter 11 solver (SCIP integer solve, expect value 4)")
if (requireNamespace("CVXR", quietly = TRUE) &&
    requireNamespace("scip", quietly = TRUE)) {
  xi <- Variable(2, integer = TRUE)
  pi <- Problem(Minimize(sum(xi)), list(xi >= 1.5))
  res <- tryCatch({
    val <- psolve(pi, solver = "SCIP")
    list(ok = is.numeric(val) && abs(val - 4) < 1e-4, detail = sprintf("value=%.4f", val))
  }, error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
  line(mark(res$ok), "SCIP", res$detail)
  if (!isTRUE(res$ok)) extra_ok <- FALSE
} else {
  line("WARN", "SCIP", "scip not installed -- needed for chapter 11")
  extra_ok <- FALSE
}

## --- 5. Chapter 14: DNLP backends present --------------------------------
hdr("Chapter 14 backends (Uno + sparsediff)")
dnlp_ok <- requireNamespace("Uno", quietly = TRUE) &&
           requireNamespace("sparsediff", quietly = TRUE)
line(if (dnlp_ok) "PASS" else "WARN", "Uno + sparsediff",
     if (dnlp_ok) "loadable" else "missing -- needed for chapter 14 (ipopt NOT needed)")
if (!dnlp_ok) extra_ok <- FALSE

## --- Summary -------------------------------------------------------------
hdr("Summary")
cat(sprintf("  CORE   (most chapters)      : %s\n", if (core_ok)  "READY" else "NOT READY"))
cat(sprintf("  EXTRAS (chapters 11 and 14) : %s\n", if (extra_ok) "READY" else "incomplete"))
if (!core_ok) {
  cat("\n  -> Fix CORE first:\n")
  cat("     install.packages(c(\"CVXR\",\"ggplot2\",\"RColorBrewer\",\"tidyr\",\"nnls\",\"glmnet\",\"png\",\"bench\"))\n")
} else if (!extra_ok) {
  cat("\n  -> Core is ready (most chapters work). For chapters 11 and 14 also:\n")
  cat("     install.packages(c(\"scip\", \"Uno\", \"sparsediff\"))\n")
} else {
  cat("\n  All set -- every chapter is ready.\n")
}
cat("\n")
