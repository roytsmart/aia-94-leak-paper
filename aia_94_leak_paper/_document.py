import pathlib

import aastex
import matplotlib.pyplot as plt

import aia_94_leak_paper

__all__ = [
    "document",
    "pdf",
]


def document() -> aastex.Document:
    """
    An :mod:`aastex` representation of the article.
    """

    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 9
    plt.rcParams["lines.linewidth"] = 1

    doc = aastex.Document(
        document_options=[
            "twocolumn",
        ],
        lmodern=False,
        textcomp=False,
    )

    doc.packages.append(aastex.Package("amsmath"))
    doc.packages.append(aastex.Package("hyperref"))
    doc.packages.append(aastex.Package("siunitx"))

    doc.preamble += aia_94_leak_paper.preamble()
    doc.preamble += aia_94_leak_paper.acronyms()

    doc.variables += aia_94_leak_paper.variables()

    doc.append(
        aastex.Title(
            aastex.NoEscape(
                r"A Growing Visible Light Leak in the 94\,\AA\ Channel "
                r"of the Atmospheric Imaging Assembly"
            )
        )
    )
    doc.append(
        aastex.ShortTitle(aastex.NoEscape(r"A Visible Light Leak in AIA 94\,\AA"))
    )

    doc += aia_94_leak_paper.authors()

    doc.append(aia_94_leak_paper.sections.abstract())
    doc.append(aia_94_leak_paper.keywords())
    doc.append(aia_94_leak_paper.sections.introduction())
    doc.append(aia_94_leak_paper.sections.noise())
    doc.append(aia_94_leak_paper.sections.eclipses())
    doc.append(aia_94_leak_paper.sections.checks())
    doc.append(aia_94_leak_paper.sections.consequences())
    doc.append(aia_94_leak_paper.sections.correction())
    doc.append(aia_94_leak_paper.sections.conclusions())

    doc.append(aastex.Bibliography("sources"))

    return doc


def pdf() -> pathlib.Path:
    """
    Build a pdf version of :func:`document` and return the path of the document.
    """

    doc = document()

    path = pathlib.Path(__file__).parent / "aia-94-leak"
    doc.generate_pdf(
        filepath=path,
        clean_tex=False,
    )

    return path.with_suffix(".pdf")
