## Shared extraction for the presenter artifacts (answer key + slide stubs).
##
## Pure text parsing of the book `NN-*.qmd` files — NO code execution, so it
## is fast and cannot fail on a solver or a package version. Sourced by both
## `build_answer_key.R` and `build_slide_stubs.R`.
##
## Structural contract (uniform across every chapter, verified 2026-07-01):
##
##   **Exercise ...            <- prompt begins
##   ...prompt lines...
##   ::: {.callout-tip collapse="true"}
##   ## Solution
##   ...solution body (code + prose)...
##   :::                       <- solution ends (fenced-div depth back to 0)

## All chapter files in book order (01-, 02-, ...). `book_dir` is the repo root.
book_chapter_files <- function(book_dir = ".") {
  files <- list.files(book_dir, pattern = "^[0-9]+.*\\.qmd$", full.names = TRUE)
  files[order(basename(files))]
}

## First `# ...` heading, with any `{#sec-...}` attribute stripped.
chapter_title <- function(lines) {
  h <- grep("^#\\s+", lines, value = TRUE)[1]
  if (is.na(h)) return(NA_character_)
  h <- sub("^#\\s+", "", h)
  h <- sub("\\s*\\{[^}]*\\}\\s*$", "", h)
  trimws(h)
}

## Level-2 section headings (`## ...`), attributes stripped. Used to seed slide
## roadmaps. Solution headings inside callouts are excluded by design because
## they live *inside* the fenced div, but we drop the literal "Solution" anyway.
section_headings <- function(lines) {
  h <- grep("^##\\s+", lines, value = TRUE)
  h <- sub("^##\\s+", "", h)
  h <- sub("\\s*\\{[^}]*\\}\\s*$", "", h)
  h <- trimws(h)
  h[h != "Solution"]
}

## Trim leading/trailing all-blank lines.
trim_blank <- function(x) {
  if (!length(x)) return(x)
  nz <- which(nzchar(trimws(x)))
  if (!length(nz)) return(character(0))
  x[seq(min(nz), max(nz))]
}

## Extract every (prompt, solution) pair from a chapter's lines.
## Returns a list of list(prompt = <chr>, solution = <chr>).
extract_exercises <- function(lines) {
  ex <- list()
  i <- 1L
  n <- length(lines)
  while (i <= n) {
    if (grepl("^\\*\\*Exercise", lines[i])) {
      ## Prompt: from the "**Exercise" line up to (not including) the callout.
      prompt <- character(0)
      while (i <= n && !grepl("^:::+\\s*\\{\\.callout-tip", lines[i])) {
        prompt <- c(prompt, lines[i])
        i <- i + 1L
      }
      if (i > n) break
      i <- i + 1L                       # step past the callout-open line
      ## Solution: collect until the fenced-div depth returns to 0.
      depth <- 1L
      sol <- character(0)
      while (i <= n && depth > 0L) {
        ln <- lines[i]
        if (grepl("^:::+\\s*\\{", ln)) {
          depth <- depth + 1L
        } else if (grepl("^:::+\\s*$", ln)) {
          depth <- depth - 1L
          if (depth == 0L) { i <- i + 1L; break }
        }
        sol <- c(sol, ln)
        i <- i + 1L
      }
      sol <- sol[!grepl("^##\\s+Solution\\s*$", sol)]
      ex[[length(ex) + 1L]] <- list(prompt = trim_blank(prompt),
                                    solution = trim_blank(sol))
    } else {
      i <- i + 1L
    }
  }
  ex
}

## Rewrite executable Quarto fences to plain highlighted fences and drop `#|`
## chunk-option lines: ```{r label} -> ```r, ```{python} -> ```python, etc.
## The result is syntax-highlighted but never executed.
clean_code <- function(lines) {
  lines <- lines[!grepl("^\\s*#\\|", lines)]
  lines <- sub("^(\\s*)```\\{([A-Za-z][A-Za-z0-9]*)[^}]*\\}\\s*$", "\\1```\\2", lines)
  lines
}
