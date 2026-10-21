"""
The figures of this article, each a factory function returning an
:class:`aastex.Figure`.
"""

from ._ingress import date_ingress, ingress, times_ingress
from ._leak_history import leak_history

__all__ = [
    "date_ingress",
    "ingress",
    "leak_history",
    "times_ingress",
]
