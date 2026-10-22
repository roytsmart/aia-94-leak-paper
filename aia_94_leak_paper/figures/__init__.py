"""
The figures of this article, each a factory function returning an
:class:`aastex.Figure`.
"""

from ._dems import dems
from ._filters import filters
from ._ingress import date_ingress, ingress, times_ingress
from ._leak_history import leak_history
from ._noise import noise

__all__ = [
    "date_ingress",
    "dems",
    "filters",
    "ingress",
    "leak_history",
    "noise",
    "times_ingress",
]
