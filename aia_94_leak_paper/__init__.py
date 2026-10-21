"""
Create the figures and compile the LaTeX files for this article.
"""

from . import dems, eclipses, figures, noise, sections, tables
from ._acronyms import acronyms
from ._authors import authors
from ._document import document, pdf
from ._keywords import keywords
from ._preamble import preamble
from ._variables import variables

__all__ = [
    "acronyms",
    "authors",
    "dems",
    "document",
    "eclipses",
    "figures",
    "keywords",
    "noise",
    "pdf",
    "preamble",
    "sections",
    "tables",
    "variables",
]
