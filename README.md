# aia-94-leak-paper

[![tests](https://github.com/roytsmart/aia-94-leak-paper/actions/workflows/tests.yml/badge.svg)](https://github.com/roytsmart/aia-94-leak-paper/actions/workflows/tests.yml)
[![ruff](https://github.com/roytsmart/aia-94-leak-paper/actions/workflows/ruff.yml/badge.svg)](https://github.com/roytsmart/aia-94-leak-paper/actions/workflows/ruff.yml)

An AAS journal article reporting a visible light leak in the 94 Å channel
of the Atmospheric Imaging Assembly (AIA) on the Solar Dynamics Observatory,
which has grown exponentially over the mission and is now most of the
quiet-Sun signal in the channel.

Roy T. Smart and Charles C. Kankelborg

📄 **[Read the article (pdf)](https://roytsmart.github.io/aia-94-leak-paper/aia-94-leak.pdf)**

## What this is

The article measures the leak four independent ways:

- **Photon noise.** Part of the signal on the disk carries no photon noise,
  and its map is a limb-darkened image of the photosphere.
- **Earth eclipses.** As SDO passes into the shadow of the Earth, the
  atmosphere absorbs the EUV long before the visible disk sets, so the leak
  can be seen alone. Two eclipse seasons a year give its history from 2010.
- **The focal-plane filters**, which attenuate it more than the EUV, and the
  304 Å channel, which shares the telescope and detector but has no leak.
- **Cross-calibration with SUVI** on GOES-16, -18, and -19.

It then follows the leak into differential emission measures, Fe XVIII maps,
and the degradation calibration of the channel, and gives a correction.

## This repository *is* the article

Every figure, table, and numeric value in the prose is computed at build time
by the `aia_94_leak_paper` package, so the text cannot drift away from the
analysis. Numbers reach the prose as `aastex.Variable` macros rather than as
literals; a quantity the analysis has not produced yet is rendered as a bold
**XX**. Sections not yet written appear in the PDF as a gray outline.

```python
import aia_94_leak_paper

aia_94_leak_paper.pdf()
```

| path | contents |
|---|---|
| `_document.py` | assembles the article and renders the pdf |
| `sections/` | the prose, as LaTeX strings, one module per section |
| `figures/` | one module per figure, each returning an `aastex.Figure` |
| `tables/` | one module per table |
| `_variables.py` | the `aastex.Variable` macros quoted in the prose |
| `_acronyms.py` | acronym definitions used as `\ACRONYM` macros |
| `_preamble.py` | custom commands, including a comment macro per author |
| `sources.bib` | the bibliography, exported from ADS |

## Building

Requires Python 3.12 or newer and a LaTeX installation. On Ubuntu:

```bash
sudo apt-get install latexmk texlive-publishers texlive-science cm-super
```

`latexmk` is not optional. Without it `pylatex` falls back to a single
`pdflatex` pass, which silently produces an article with no bibliography.

```bash
pip install -e .[test]
pytest              # builds and validates the pdf
ruff check          # linting, enforced in CI
ruff format         # formatting, enforced in CI
```

## Continuous integration

Every push runs the test suite, which builds the article and checks that no
citation or cross-reference is left unresolved. Pushes to `main` publish the
article to the link above, and pull requests get their own preview at
`…/pr/<number>/aia-94-leak.pdf`, linked in a comment on the pull request and
removed when it closes.

## Commenting on the draft

Each author has a colored comment macro for notes in the prose:
`\roy{...}` and `\charles{...}`.
