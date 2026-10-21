import aastex

__all__ = [
    "keywords",
]


def keywords() -> aastex.Keywords:
    """
    The Unified Astronomy Thesaurus concepts describing this article.
    """
    return aastex.Keywords(
        [
            aastex.UAT("Solar instruments", 1499),
            aastex.UAT("Calibration", 2179),
            aastex.UAT("Solar extreme ultraviolet emission", 1493),
            aastex.UAT("Solar corona", 1483),
            aastex.UAT("Quiet sun", 1322),
            aastex.UAT("Solar coronal holes", 1484),
            aastex.UAT("Eclipses", 442),
            aastex.UAT("Limb darkening", 922),
        ]
    )
