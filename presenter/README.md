# Presenter materials

Presenter-facing artifacts for the useR! 2026 CVXR tutorial. Both are
**auto-derived from the book chapters** (`../NN-*.qmd`), so they never drift:
re-run the generator after any chapter edit. Neither is part of the book render.

## 1. Answer key — `answer-key.html`

A single standalone page with **every exercise + solution** in chapter order —
the presenter safety net for "I forgot the solution." The same solutions are
already click-to-reveal in the book; this just gathers them, order-independent.

```sh
Rscript presenter/build_answer_key.R      # run from the repo root
```

- One self-contained HTML (`--embed-resources`), works offline; print to PDF
  from the browser for a paper backup.
- Math via MathJax (degrades to readable LaTeX offline).
- **Code + prose only — not executed.** Robust; cannot fail on a solver. The
  live book covers "watch it run."
- Intermediate `answer-key.md` is gitignored; only the HTML + script are tracked.

## 2. Slide stubs — `../slides/NN-*.qmd`

One thin `revealjs` deck per chapter, mirroring the book 1:1. Deliberately
near-empty — a title slide is a valid deck. Each is seeded with a title slide, a
**Roadmap** stub (bullets = the chapter's `##` section headings), and one slide
per exercise: **prompt on the slide, solution in `::: {.notes}` speaker notes.**

```sh
Rscript presenter/build_slide_stubs.R           # scaffold; SKIPS existing decks
Rscript presenter/build_slide_stubs.R --force   # overwrite all
```

- **Scaffold-once, then hand-edit.** The generator skips any deck that already
  exists so your rehearsal edits are safe; `--force` regenerates from scratch.
- Decks live in `slides/` alongside the opening `intro.qmd` and `slides.scss`.
  The **book render ignores the folder**; rendered `*.html` is gitignored — the
  `.qmd` stubs are the tracked deliverable.
- Render one: `quarto render slides/02-gentle-intro.qmd --to revealjs`.
- Note: `slides/intro.qmd` is the separate **session-opening** deck (why-here →
  Hello World); `slides/01-intro.qmd` is the chapter-1 stub. Merge or drop one on
  rehearsal as you prefer.

## 3. Full follow-along decks — `../slides/full/NN-*.qmd`

The **whole book, as slides** — one `revealjs` deck per chapter carrying the
*entire* chapter (same content, same order, code executed), for presenting
from a single medium instead of flipping between the web pages and the thin
stubs. Big font; long slides scroll (`scrollable`).

```sh
Rscript presenter/build_full_slides.R      # regenerate all 16 .qmd + index.html
```

- **Fully mechanical, re-runnable.** Each deck is the chapter verbatim with a
  faithful transform: reveal front-matter swapped in, the `# Title {#sec-id}` H1
  lifted to the deck title, cross-chapter `@sec-*` refs rewritten to the target
  chapter's title, `bibliography:` pointed at the repo bibs, and knitr
  `root.dir` set to the repo root so `data/`, `images/`, `_check_status.R`
  resolve. Nothing else is edited — output, exercises, and prose match the book.
- **Render one:** `quarto render slides/full/02-gentle-intro.qmd --to revealjs`.
  Freeze is format-specific, so these re-execute against CRAN CVXR 1.9.1.
- **Launcher:** `slides/full/index.html` links the session opener + all 16 decks.
- Rendered `*.html` is gitignored (regenerated); the `.qmd` are the tracked
  deliverable. The book render ignores `slides/`, so no interference.
- This complements — does not replace — the thin stubs in `slides/` and the
  opening `slides/intro.qmd`; pick whichever fits how you want to present.

## Shared

`extract_book.R` — the text-only parser both stub/answer-key generators source
(chapter titles, section headings, and `(prompt, solution)` pairs). No code
execution. (`build_full_slides.R` is self-contained — it does not use it.)
