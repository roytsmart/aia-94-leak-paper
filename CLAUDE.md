# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`aia-94-leak-paper` is a **reproducible scientific article**, not a conventional
software library. The `aia_94_leak_paper` Python package programmatically generates
a complete AAS-journal LaTeX article (text, figures, tables, and numeric values)
reporting a visible light leak in the 94 Å channel of SDO/AIA.
Calling `aia_94_leak_paper.pdf()` produces the final `aia-94-leak.pdf`.

This is one package within the larger Kankelborg-Group workspace (see the parent
`../CLAUDE.md`). It follows the layout of `esis2-science-paper`, without that
repo's Sphinx docs, and with ruff for both linting and formatting.

## The article

The leak is measured four independent ways, and each gets its own section:

- **Photon noise.** On the disk, part of the 94 Å signal carries no photon
  noise; its map is a limb-darkened image of the photosphere.
- **Earth eclipses.** During ingress the atmosphere absorbs the EUV well before
  the visible disk sets, so the leak is seen alone. One ingress per season
  gives its history since 2010.
- **Independent checks.** The thick focal-plane filter (1000 Å more Zr)
  attenuates it only as much as the EUV, so it passes that filter through
  openings rather than the metal; the 304 Å channel (same telescope and CCD)
  has none; SUVI on GOES-16, -18, and -19 gives the same history.
- **Consequences and a correction.** DEMs (hot quiet Sun and coronal holes),
  Fe XVIII maps, the degradation calibration, and published results that show
  its symptoms.

Roy T. Smart is the corresponding author, then Charles C. Kankelborg.

## Rules

- **Every number in the prose is an `aastex.Variable`** defined in `_variables.py`.
  A number the analysis cannot supply yet is declared with `_pending()`, which
  renders as a bold `XX`, so the draft shows plainly what is still owed.
- Sections not yet written are drafted as an `outline` environment (defined in
  `_preamble.py`), which renders as a gray bulleted list, so the plan cannot be
  mistaken for the article.
- **Where the analysis lives.** Analysis specific to this article (the noise
  fits, the eclipse forward model, the SUVI cross-calibration) lives here as
  thin modules. Anything reusable may move down the stack (`sdo`, `utu`,
  `named-arrays`) only after the user agrees to the new library API.
- **Data.** AIA level 1 images come from JSOC (`aia.lev1_euv_12s`, by SUMS path
  through `drms`), the synoptic level 1.5 images for the DEMs from
  `jsoc1.stanford.edu/data/aia/synoptic`, AIA flat fields from
  `aia.flatfield`, and SUVI level 2 composites from NOAA NCEI. Cache
  downloads, and keep the number of parallel requests to JSOC small (it
  answers 12 at once with HTTP 503).
- **The archives.** Each of `eclipses` (the ingresses, binned 16 x 16),
  `noise` (the noon runs, reduced to statistics per region), `filters` (a run
  through each focal-plane filter on 2026-03-07, reduced the same way), and
  `dems` (the six coronal channels at noon, binned to 4.8 arcsec, with the
  94 Å flat field) is built from JSOC by `python -m aia_94_leak_paper.<name> build` (needs
  `.[data]`), bundled by `... pack`, and published as an asset of the release
  named by `<name>.version`; the article reads them through `<name>.load`,
  which downloads the asset once and checks its SHA-256. A change to what a
  builder writes is a new release and a new `version`/`sha256`, never a
  replaced asset.

## Commands

Run from this package directory:

```bash
pip install -e .[test]          # install for development
pytest                          # run tests; test_pdf compiles the LaTeX → PDF
ruff format                     # format (CI enforces --check)
ruff check                      # lint (CI enforces)
```

Building the PDF requires a **LaTeX installation** with `latexmk` (TeX Live:
`texlive-publishers texlive-science cm-super latexmk`). Without `latexmk`,
`pylatex` falls back to a single `pdflatex` pass and the bibliography is
silently missing; `test_pdf` checks for this.

## Architecture

The whole article is assembled in `_document.py`: `document()` builds an
`aastex.Document`, appending acronyms, variables, title, authors, abstract,
keywords, each section, and finally the bibliography (`sources.bib`). `pdf()`
renders it.

- **`sections/`**: each module returns an `aastex.Section` (or `aastex.Abstract`).
  Prose lives here as raw LaTeX strings (with `\cite`, `\ref`, equations) and
  embeds figures/tables by calling the corresponding factory.
- **`figures/`**: each module builds a matplotlib figure and wraps it in an
  `aastex.Figure`/`FigureStar` with a caption.
- **`tables/`**: factories returning `aastex.Table` objects.
- **`_variables.py`**: `aastex.Variable` LaTeX macros for every numeric value
  cited in the prose.
- **`_acronyms.py`**: `aastex.Acronym` definitions; prose uses `\ACRONYM` macros,
  and `\ACRONYMCapital` where the acronym opens a sentence.
- **`_preamble.py`**: custom commands: ion macros such as `\FeXVIII`, and a
  colored comment macro per author (`\roy`, `\charles`).

## Conventions

- Every module declares an explicit `__all__` and exposes functionality through small
  factory functions, re-exported up the package via `__init__.py`. To add a new
  figure/table/section: create `_name.py`, then add it to the subpackage `__init__.py`.
- Modules import the top-level package as `import aia_94_leak_paper` and reach back
  into it (`aia_94_leak_paper.figures.x()`) rather than deep relative imports.
- House style: `Figure~\ref`, `Section~\ref`, `Table~\ref`; `\eg,` with a comma;
  a tie `~` before every `\citep`; no em dashes in the prose.
- Cite only entries in `sources.bib`. Add new ones from the ADS BibTeX export
  (`/v1/export/bibtex`, token in `~/.ads/dev_key.txt`, never printed), keyed as
  first author and year (`Lemen2012`), rather than typing them by hand.
- Generated LaTeX/PDF output and caches are gitignored; never commit build artifacts.
- Stage explicit files with `git add`; never `git add` whole directories.
