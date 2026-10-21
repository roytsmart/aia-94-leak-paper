import aastex

import aia_94_leak_paper

__all__ = [
    "variables",
]


def _pending(name: str) -> aastex.Variable:
    """
    A quantity the analysis has not produced yet, marked where it is cited.

    Every number in this article is computed by the analysis, and these are
    the ones it has no answer for yet. They are written into the text as a
    bold ``XX`` rather than guessed at or left out, so that a draft shows
    plainly what is still owed and no placeholder can be mistaken for a
    measurement. ``??`` is avoided since it is what LaTeX prints for an
    unresolved cross-reference, which the tests look for.
    """
    return aastex.Variable(name=name, value=aastex.NoEscape(r"\textbf{XX}"))


def variables() -> list[aastex.Variable]:
    """
    A list of LaTeX variables for every numeric quantity cited in the prose.

    Reference these macros in section strings instead of hardcoding numbers so
    that the text stays in sync with the analysis.
    """
    dates = aia_94_leak_paper.eclipses.dates
    return [
        # the eclipses the leak is measured in, and the years they span
        aastex.Variable(name="numEclipses", value=len(dates)),
        aastex.Variable(name="yearFirst", value=min(dates).year),
        aastex.Variable(name="yearLast", value=max(dates).year),
        # the leak at disk center in the first and the last eclipse
        _pending("leakFirst"),
        _pending("leakLast"),
        # the year the leak began to grow exponentially, and how fast
        _pending("yearOnset"),
        _pending("doublingTime"),
        # how long after the EUV is absorbed the visible disk sets
        _pending("sunsetDelay"),
        # the fraction of the quiet-Sun signal at disk center that is leak
        _pending("leakFractionQuietSun"),
    ]
