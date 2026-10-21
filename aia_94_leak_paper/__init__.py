"""
Create the figures and compile the LaTeX files for this article.
"""

from . import eclipses, figures, sections, tables
from ._acronyms import acronyms
from ._authors import authors
from ._document import document, pdf
from ._keywords import keywords
from ._preamble import preamble
from ._variables import variables

__all__ = [
    "acronyms",
    "authors",
    "document",
    "eclipses",
    "figures",
    "keywords",
    "pdf",
    "preamble",
    "sections",
    "tables",
    "variables",
]
