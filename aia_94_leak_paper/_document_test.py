import pathlib

import pylatex
import pymupdf

import aia_94_leak_paper


def test_document():
    doc = aia_94_leak_paper.document()
    assert isinstance(doc, pylatex.Document)


def test_pdf(capsys):
    with capsys.disabled():
        pdf = aia_94_leak_paper.pdf()
    assert isinstance(pdf, pathlib.Path)
    assert pdf.exists()

    with pymupdf.open(pdf) as document:
        text = "".join(page.get_text() for page in document)

    # The bibliography and the cross-references are only resolved if the
    # article is compiled repeatedly with BibTeX in between, which `latexmk`
    # does and a bare `pdflatex` does not. Without this check the article
    # builds "successfully" with every citation rendered as `(?)` and no
    # reference list at all.
    assert "(?)" not in text
    assert "??" not in text

    bbl = pdf.with_suffix(".bbl")
    assert bbl.exists()
    assert "bibitem" in bbl.read_text(encoding="utf-8", errors="replace")
