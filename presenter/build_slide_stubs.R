#!/usr/bin/env Rscript
## Scaffold thin per-chapter revealjs SLIDE STUBS, one deck per book chapter,
## mirroring the book 1:1 (PRESENTER_PLAN.md, Artifact 2).
##
##   Scaffold (from the repo root):
##       Rscript presenter/build_slide_stubs.R          # skip decks that exist
##       Rscript presenter/build_slide_stubs.R --force  # overwrite all
##
## Each deck is deliberately near-empty — a title slide is a valid deck. It is
## seeded with:
##   * a section title slide (chapter title),
##   * a "Roadmap" stub whose bullets are the chapter's `## ` section headings
##     (skeleton to flesh out on rehearsal),
##   * one slide per exercise: the PROMPT is on the slide (audience) and the
##     SOLUTION is in `::: {.notes}` speaker notes (presenter-only) — the
##     elegant complement to the answer key.
##
## Decks live in slides/ (same folder as the opening intro.qmd + slides.scss);
## the book render ignores the folder. Rendered HTML is gitignored; the .qmd
## stubs are the tracked deliverable, meant to be hand-edited — so this script
## SKIPS any deck that already exists unless --force is given.

source("presenter/extract_book.R")

force  <- "--force" %in% commandArgs(trailingOnly = TRUE)
files  <- book_chapter_files(".")
if (!length(files)) stop("No chapter files found — run from the repo root.")

frontmatter <- function(title, ch_no) {
  c("---",
    sprintf('title: "%s"', title),
    sprintf('subtitle: "useR! 2026 · Warsaw · Chapter %d"', ch_no),
    "format:",
    "  revealjs:",
    "    theme: [default, slides.scss]",
    "    slide-number: c/t",
    "    width: 1280",
    "    height: 720",
    "    link-external-newwindow: true",
    '    footer: "[cvxr.rbind.io](https://cvxr.rbind.io)"',
    "execute:",
    "  echo: true",
    "  eval: false        # stub deck: code is shown, not run",
    "  warning: false",
    "  message: false",
    "---",
    "")
}

created <- character(0)
skipped <- character(0)

for (ch_no in seq_along(files)) {
  f      <- files[ch_no]
  base   <- basename(f)
  out    <- file.path("slides", base)
  if (file.exists(out) && !force) { skipped <- c(skipped, out); next }

  lines  <- readLines(f, warn = FALSE)
  title  <- chapter_title(lines)
  heads  <- section_headings(lines)
  ex     <- extract_exercises(lines)

  body <- frontmatter(title, ch_no)

  ## Title / section-break slide.
  body <- c(body,
            sprintf("# %s", title),
            "",
            "<!-- STUB deck — flesh out on rehearsal. Signposts + visuals only;",
            "     do NOT re-teach the book. Many chapters need only this slide. -->",
            "")

  ## Roadmap seeded from the chapter's section headings.
  if (length(heads)) {
    body <- c(body,
              "## Roadmap {.smaller}", "",
              "::: incremental",
              paste0("- ", heads),
              ":::", "")
  }

  ## One slide per exercise: prompt visible, solution in speaker notes.
  for (k in seq_along(ex)) {
    slide_title <- if (length(ex) > 1L) sprintf("## Your turn — Exercise %d {.smaller}", k)
                   else "## Your turn {.smaller}"
    body <- c(body,
              slide_title, "",
              ex[[k]]$prompt, "",
              "::: {.notes}",
              "**Solution (presenter notes).**", "",
              clean_code(ex[[k]]$solution),
              ":::", "")
  }

  writeLines(body, out)
  created <- c(created, out)
}

cat(sprintf("Slide stubs: %d created, %d skipped (already present).\n",
            length(created), length(skipped)))
if (length(created)) cat("  created:\n", paste0("    ", created, collapse = "\n"), "\n", sep = "")
if (length(skipped)) cat("  skipped:\n", paste0("    ", skipped, collapse = "\n"), "\n", sep = "")
