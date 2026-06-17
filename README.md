# Convex Optimization in R Using CVXR

A hands-on tutorial on disciplined convex optimization in R with
[CVXR](https://cvxr.rbind.io), prepared for **useR! 2026** (Warsaw).

Authors: **Anqi Fu, Balasubramanian Narasimhan, Stephen Boyd.**

This is a [Quarto book](https://quarto.org/docs/books/). It targets **CVXR 1.9.1**
and uses **CRAN solvers only**.

> The earlier (2019, CVXR 1.x, bookdown) tutorial is preserved on the
> [`useR2019`](https://github.com/bnaras/cvxr_tutorial/tree/useR2019) branch.

## For participants

1. **Install R ≥ 4.3.0** and a recent IDE ---
   [RStudio](https://posit.co/download/rstudio-desktop/) or
   [Positron](https://positron.posit.co).
2. **Get the materials** (this repository):

   - **RStudio:** *File → New Project → Version Control → Git*, paste
     `https://github.com/bnaras/cvxr_tutorial.git`, then open `cvxr-tutorial-2026.Rproj`.
   - **Positron:** Command Palette → *Git: Clone* → paste the same URL → open the folder.
   - **Terminal:** `git clone https://github.com/bnaras/cvxr_tutorial.git`
   - **No Git:** *Code → Download ZIP* on the repo page and unzip.

   You also need the [Quarto CLI](https://quarto.org/docs/get-started/) (≥ 1.4) to
   render chapters or the whole book. RStudio bundles Quarto; Positron and terminal
   users install it once from that link. (Running individual code chunks needs no
   render step.)
3. Install the packages (CRAN only):

   ```r
   # Solvers (CVXR auto-pulls clarabel, osqp, scs, highs)
   install.packages(c("CVXR", "scip", "Uno", "sparsediff"))
   # Helpers used by the examples
   install.packages(c("ggplot2", "tidyr", "nnls", "glmnet",
                      "boot", "png", "bench"))
   ```

4. **Check your setup:**

   ```r
   source("setup_check.R")     # or:  Rscript setup_check.R
   ```

   It verifies your R version and packages and runs a tiny solve on each solver,
   reporting `PASS`/`FAIL`. The **core hands-on chapters need only `CVXR` + the
   helpers**; `scip`/`Uno`/`sparsediff` power two end-of-book demonstrations.

5. Work through the chapters by **running the code chunks** in RStudio.

## Building the book (authors)

```r
# install.packages("quarto")  # the R wrapper, if needed
quarto::quarto_render()        # or, in a terminal:  quarto render
```

Requires the [Quarto CLI](https://quarto.org/docs/get-started/) (≥ 1.4). For PDF
output you also need a LaTeX install (e.g. `tinytex::install_tinytex()`); HTML
needs no LaTeX. Rendered output goes to `_book/`.

## Repository layout

| Path | What |
|------|------|
| `index.qmd` | Preface |
| `01-…`–`16-…qmd` | Chapters (see `_quarto.yml` for the structure) |
| `setup_check.R` | Pre-tutorial environment check |
| `_check_status.R` | Helper sourced by example chunks |
| `references.bib`, `cvxr.bib` | Bibliography |
| `custom.scss` | Theme (shared with the CVXR website) |
| `figures/` | Reused figures |

## License

See [`LICENSE`](LICENSE).
