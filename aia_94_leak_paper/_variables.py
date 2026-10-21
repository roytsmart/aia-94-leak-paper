import aastex
import astropy.units as u
import numpy as np

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


def _significant(value: u.Quantity, digits: int = 2) -> u.Quantity:
    """A quantity rounded to some significant figures, which is all a measurement here supports."""
    if value.value == 0:
        return value
    decimals = digits - 1 - int(np.floor(np.log10(abs(value.value))))
    return np.round(value, decimals)


def variables() -> list[aastex.Variable]:
    """
    A list of LaTeX variables for every numeric quantity cited in the prose.

    Reference these macros in section strings instead of hardcoding numbers so
    that the text stays in sync with the analysis.
    """
    eclipses = aia_94_leak_paper.eclipses
    date_ingress = aia_94_leak_paper.figures.date_ingress
    dates = eclipses.dates
    leaks = eclipses.history()
    growth = eclipses.growth(leaks)
    rate = (u.DN, u.s**-1)
    first, last = leaks[0], leaks[-1]
    fraction = (last.center / last.quiet).to_value(u.percent)
    noise = aia_94_leak_paper.noise.fit(date_ingress)
    disk = aia_94_leak_paper.noise.disk(date_ingress)
    # the noise method measures the leak once the leak is three times the
    # typical spread of the method among the quadrants of the detector
    history = aia_94_leak_paper.noise.history()
    spread = aia_94_leak_paper.noise.spread()
    measured = [
        (n, k)
        for n, k, m in zip(
            history, leaks, aia_94_leak_paper.noise.measured(history, leaks)
        )
        if m
    ]
    ratios = [(n.rate / k.center).to_value(u.percent) for n, k in measured]
    return [
        # the eclipses the leak is measured in, and the years they span
        aastex.Variable(name="numEclipses", value=len(dates)),
        aastex.Variable(name="yearFirst", value=min(dates).year),
        aastex.Variable(name="yearLast", value=max(dates).year),
        # the leak at disk center in the first and the last eclipse
        aastex.Variable(name="leakFirst", value=_significant(first.center), unit=rate),
        aastex.Variable(name="leakLast", value=_significant(last.center), unit=rate),
        # the year the leak began to grow exponentially, and how fast
        aastex.Variable(name="yearOnset", value=growth.since),
        aastex.Variable(name="doublingTime", value=_significant(growth.doubling)),
        # how long after the EUV edge the visible disk is half set, the median over the
        # eclipses since the onset (before it the setting is held at this median)
        aastex.Variable(
            name="sunsetDelay",
            value=_significant(
                np.median(
                    u.Quantity(
                        [k.sunset for k in leaks if k.date.year >= eclipses.onset]
                    )
                )
            ),
        ),
        # the eclipse whose ingress is shown, and when its visible disk set after the EUV edge
        aastex.Variable(
            name="ingressDate",
            value=aastex.NoEscape(
                f"{date_ingress:%Y} {date_ingress:%B} {date_ingress.day}"
            ),
        ),
        aastex.Variable(
            name="ingressSunset",
            value=_significant(eclipses.fit(date_ingress).sunset),
        ),
        # the fraction of the signal at disk center that is leak, in the last eclipse
        aastex.Variable(
            name="leakFractionQuietSun",
            value=aastex.NoEscape(rf"{_significant(fraction * u.one).value:g}\%"),
        ),
        # the noise of the run on the day of the ingress shown
        aastex.Variable(name="noiseFrames", value=aia_94_leak_paper.noise.frames),
        aastex.Variable(name="readNoise", value=_significant(np.nanmean(noise.read))),
        aastex.Variable(
            name="noiseSlope", value=_significant(np.nanmean(noise.slope) * u.DN)
        ),
        aastex.Variable(name="noiselessCenter", value=_significant(noise.noiseless)),
        aastex.Variable(
            name="noiselessSpread", value=_significant(noise.spread, digits=1)
        ),
        # the map of the noiseless signal on that day: its limb darkening, how
        # much it grows per DN of the 94 Å signal of a region, and the part of
        # that signal it is at disk center
        aastex.Variable(
            name="noiseLimbDarkening",
            value=aastex.NoEscape(
                rf"\ensuremath{{{_significant(disk.limb_darkening * u.one).value:g}}}"
            ),
        ),
        aastex.Variable(
            name="noisePerSignal",
            value=aastex.NoEscape(
                rf"\ensuremath{{{_significant(disk.per_signal * u.one, digits=1).value:g}}}"
            ),
        ),
        aastex.Variable(
            name="noiseShare",
            value=aastex.NoEscape(
                rf"{_significant((disk.amplitude / disk.signal).to(u.percent)).value:g}\%"
            ),
        ),
        # the typical spread of the noise method among the quadrants, and how
        # its noiseless signal compares with the eclipses once the leak is
        # three times that
        aastex.Variable(name="noiseSpread", value=_significant(spread, digits=1)),
        aastex.Variable(name="yearNoise", value=min(k.date.year for _, k in measured)),
        aastex.Variable(
            name="noiseRatioMin",
            value=aastex.NoEscape(rf"{_significant(min(ratios) * u.one).value:g}\%"),
        ),
        aastex.Variable(
            name="noiseRatioMax",
            value=aastex.NoEscape(rf"{_significant(max(ratios) * u.one).value:g}\%"),
        ),
    ]
