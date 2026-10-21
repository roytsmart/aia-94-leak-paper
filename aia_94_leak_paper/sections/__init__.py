"""
The sections of this article, each a factory function returning an
:class:`aastex.Section`.
"""

from ._s0_abstract import abstract
from ._s1_introduction import introduction
from ._s2_noise import noise
from ._s3_eclipses import eclipses
from ._s4_checks import checks
from ._s5_consequences import consequences
from ._s6_correction import correction
from ._s7_conclusions import conclusions

__all__ = [
    "abstract",
    "checks",
    "conclusions",
    "consequences",
    "correction",
    "eclipses",
    "introduction",
    "noise",
]
