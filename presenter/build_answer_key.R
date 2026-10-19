#!/usr/bin/env Rscript
## Build the presenter-only ANSWER KEY: a single standalone HTML page holding
## every exercise + solution in the book, in chapter order.
##
##   Regenerate (from the repo root):
##       Rscript presenter/build_answer_key.R
##
## Output: presenter/answer-key.html — standalone, self-contained (one file,
## works offline), math via MathJax (degrades to readable LaTeX offline).
## Decision (PRESENTER_PLAN.md, Artifact 1): code + prose only, NOT executed
## numeric outputs — robust and cannot fail on a solver. The live book covers
## "watch it run".
##
## Presenter-only, but harmless if published: it only duplicates solutions that
## are already click-to-reveal in the book. Not a `.qmd`, so the book render
## ignores it.

source("presenter/extract_book.R")

files <- book_chapter_files(".")
if (!length(files)) stop("No chapter files found — run from the repo root.")

md <- c(
  "> **Presenter answer key.** Every solution below is also a click-to-reveal",
  "> `callout-tip` in the book. Code is shown but **not executed** — run the",
  "> live book to watch each problem solve.",
  ""
)

n_ex <- 0L
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  ex <- extract_exercises(lines)
  if (!length(ex)) next
  md <- c(md, "", paste0("# ", chapter_title(lines)), "")
  for (k in seq_along(ex)) {
    n_ex <- n_ex + 1L
    md <- c(md,
            ex[[k]]$prompt, "",
            "**Solution.**", "",
            clean_code(ex[[k]]$solution),
            "", "---", "")
  }
}

md_file  <- "presenter/answer-key.md"
out_html <- "presenter/answer-key.html"
writeLines(md, md_file)

pandoc <- Sys.which("pandoc")
if (!nzchar(pandoc)) stop("pandoc not found on PATH.")

status <- system2(
  pandoc,
  c(shQuote(md_file), "-o", shQuote(out_html),
    "--standalone", "--toc", "--toc-depth=1",
    "--mathjax", "--embed-resources",
    "--metadata", shQuote("title=Presenter Answer Key — CVXR useR! 2026")),
  stdout = TRUE, stderr = TRUE
)
if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
  cat(status, sep = "\n")
  stop("pandoc failed.")
}

cat(sprintf("Wrote %s (%d exercises across %d chapters).\n",
            out_html, n_ex, sum(vapply(files, function(f)
              length(extract_exercises(readLines(f, warn = FALSE))) > 0L, logical(1)))))
