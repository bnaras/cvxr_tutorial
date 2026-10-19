#!/usr/bin/env Rscript
## Generate FULL per-chapter revealjs decks from the book chapters: the ENTIRE
## chapter, same content, same order, rendered as slides ("the book, as slides").
## This is the follow-along presentation companion to the thin stubs in slides/.
##
##   Rscript presenter/build_full_slides.R      # run from the repo root; regenerates all
##
## Output: slides/full/NN-*.qmd  (rendered *.html is gitignored; the .qmd are tracked).
## The book render ignores slides/, so these never interfere with the book.
##
## Faithful mechanical transform (identical for every chapter):
##   1. Swap the book front-matter for a revealjs one (scrollable, big font, echo,
##      execute, footer, theme ../slides.scss).
##   2. Drop the `# Title {#sec-id}` H1  ->  becomes the deck title / subtitle.
##   3. Rewrite cross-chapter `@sec-*` refs to the target chapter's title in italics
##      (they cannot resolve in a standalone deck; every @sec- id here is a chapter).
##   4. Point `bibliography:` at the repo bibs so `@citations` resolve.
##   5. Set knitr root.dir to the repo root so `source("_check_status.R")`, `data/`,
##      and `images/` resolve when chunks execute.
## Nothing else is altered: the body is the chapter verbatim, so executed output,
## exercises (click-to-reveal solutions), and prose match the book exactly.

repo  <- normalizePath(getwd())
files <- sort(list.files(".", pattern = "^[0-9]{2}-.*\\.qmd$"))
if (!length(files)) stop("No chapter files found — run from the repo root.")

## --- pass 1: build the sec-id -> chapter-title map (for cross-ref rewriting) ------
h1_line   <- function(lines) lines[grep("^# ", lines)[1]]
title_of  <- function(lines) trimws(sub("\\s*\\{#[^}]*\\}\\s*$", "", sub("^# ", "", h1_line(lines))))
secid_of  <- function(lines) sub(".*\\{#(sec-[^}]*)\\}.*", "\\1", h1_line(lines))

titles <- secids <- character(length(files))
for (i in seq_along(files)) {
  L <- readLines(files[i], warn = FALSE)
  titles[i] <- title_of(L)
  secids[i] <- secid_of(L)
}
sec_map <- setNames(titles, secids)          # "sec-dcp" -> "Disciplined Convex Programming"

frontmatter <- function(title, ch_no) {
  c("---",
    sprintf('title: "%s"', title),
    sprintf('subtitle: "useR! 2026 · Warsaw · Chapter %d"', ch_no),
    "format:",
    "  revealjs:",
    "    theme: [default, ../slides.scss]",
    "    slide-number: c/t",
    "    width: 1280",
    "    height: 720",
    "    scrollable: true",
    "    code-line-numbers: false",
    "    link-external-newwindow: true",
    '    footer: "[cvxr.rbind.io](https://cvxr.rbind.io)"',
    "execute:",
    "  echo: true",
    "  eval: true",
    "  warning: false",
    "  message: false",
    "  freeze: false",
    "bibliography: [../../references.bib, ../../cvxr.bib]",
    "knitr:",
    "  opts_knit:",
    sprintf('    root.dir: "%s"', repo),
    "---",
    "")
}

## Rewrite every @sec-<id> to *Chapter Title* (italic), matching only whole ids.
rewrite_secrefs <- function(text) {
  for (id in names(sec_map)) {
    pat <- paste0("@", id, "(?![A-Za-z0-9-])")   # id boundary via negative lookahead
    text <- gsub(pat, paste0("*", sec_map[[id]], "*"), text, perl = TRUE)
  }
  text
}

dir.create(file.path("slides", "full"), showWarnings = FALSE, recursive = TRUE)
created <- character(0)

for (ch_no in seq_along(files)) {
  f     <- files[ch_no]
  lines <- readLines(f, warn = FALSE)

  h1   <- grep("^# ", lines)[1]                 # drop the H1 title line only
  body <- lines[-h1]
  body <- rewrite_secrefs(body)

  out  <- file.path("slides", "full", basename(f))
  writeLines(c(frontmatter(titles[ch_no], ch_no), body), out)
  created <- c(created, out)
}

## --- launcher: one page linking every deck (regenerated; gitignored like the decks) --
li <- sprintf(
  '    <li><a href="%s">%02d &middot; %s</a></li>',
  sub("\\.qmd$", ".html", basename(files)), seq_along(files), titles)
index <- c(
  '<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>CVXR Tutorial — Full Slide Decks</title>',
  '<style>',
  ' body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;',
  '   max-width:760px;margin:3rem auto;padding:0 1.2rem;line-height:1.6;color:#222}',
  ' h1{color:#1a5276;font-size:1.6rem;margin-bottom:.2rem}',
  ' p.sub{color:#5d6d7e;margin-top:0}',
  ' ol{padding-left:0;list-style:none;counter-reset:none}',
  ' li{margin:.35rem 0;font-size:1.15rem}',
  ' a{color:#2471a3;text-decoration:none} a:hover{text-decoration:underline}',
  ' .opener{margin:1rem 0 1.6rem;font-size:1.2rem}',
  ' footer{margin-top:2rem;color:#5d6d7e;font-size:.85rem}',
  ' @media(prefers-color-scheme:dark){body{background:#111;color:#ddd}',
  '   h1{color:#5dade2}a{color:#5dade2}}',
  '</style></head><body>',
  '<h1>Convex Optimization in R Using CVXR</h1>',
  '<p class="sub">useR! 2026 &middot; Warsaw &middot; full follow-along slide decks</p>',
  '<p class="opener"><a href="../intro.html">&#9654; Session opener</a> &mdash; start here</p>',
  '<ol>',
  li,
  '</ol>',
  '<footer>Regenerate: <code>Rscript presenter/build_full_slides.R</code>, then render each deck.<br>',
  'Home &amp; all examples: <a href="https://cvxr.rbind.io">cvxr.rbind.io</a></footer>',
  '</body></html>')
writeLines(index, file.path("slides", "full", "index.html"))

cat(sprintf("Full slide decks: %d written to slides/full/ (+ index.html launcher).\n", length(created)))
cat(paste0("    ", created, collapse = "\n"), "\n", sep = "")
