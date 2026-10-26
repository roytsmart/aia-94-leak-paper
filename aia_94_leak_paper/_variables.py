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


def _figures(value: float, digits: int = 2) -> str:
    """
    A number written to some significant figures, which is all a measurement
    here supports.

    The zeros among the figures are written too: a number rounded and then
    printed as a float loses them, and 0.90 written as 0.9 claims one figure
    rather than two.
    """
    mantissa, exponent = f"{value:.{digits - 1}e}".split("e")
    decimals = max(digits - 1 - int(exponent), 0)
    return f"{float(f'{mantissa}e{exponent}'):.{decimals}f}"


def _measured(value: u.Quantity, digits: int = 2) -> aastex.NoEscape:
    """A quantity written to some significant figures, with its unit."""
    text = value.to_string(
        format="latex_inline",
        formatter=lambda v: _figures(v, digits),
    )
    return aastex.NoEscape(rf"\ensuremath{{{text[1:~0]}}}")


def _percent(fraction: float) -> aastex.NoEscape:
    """A fraction written as a percentage, to two significant figures."""
    return aastex.NoEscape(rf"{_figures(fraction * 100)}\%")


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
    rate = u.DN / u.s
    first, last = leaks[0], leaks[-1]
    fraction = (last.center / last.quiet).to_value(u.percent)
    gradient = eclipses.residuals(last.date).gradient
    run = aia_94_leak_paper.noise.load(date_ingress)
    noise = aia_94_leak_paper.noise.fit(run)
    disk = aia_94_leak_paper.noise.disk(run)
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
    filters = aia_94_leak_paper.filters
    attenuation = filters.attenuation()
    gap = filters.load("thin").time - filters.load("thick").time
    entrance = aia_94_leak_paper.entrance
    monitor = entrance.load()
    bright = entrance.brightness(monitor)
    comparison = entrance.compare()
    first_dominant = int(np.flatnonzero(comparison.dominant)[0])
    dems = aia_94_leak_paper.dems
    date_dem = dems.date_dem
    shown = dems.summary(date_dem)
    summaries = dems.history()
    early = [s for s in summaries if s.date.year < eclipses.onset]
    return [
        # the eclipses the leak is measured in, and the years they span
        aastex.Variable(name="numEclipses", value=len(dates)),
        aastex.Variable(name="yearFirst", value=min(dates).year),
        aastex.Variable(name="yearLast", value=max(dates).year),
        # the leak at disk center in the first and the last eclipse
        aastex.Variable(name="leakFirst", value=_measured(first.center.to(rate))),
        aastex.Variable(name="leakLast", value=_measured(last.center.to(rate))),
        # the year the leak began to grow exponentially, and how fast
        aastex.Variable(name="yearOnset", value=growth.since),
        aastex.Variable(name="doublingTime", value=_measured(growth.doubling)),
        # how long after the EUV edge the visible disk is half set, the median over the
        # eclipses since the onset (before it the setting is held at this median)
        aastex.Variable(
            name="sunsetDelay",
            value=_measured(
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
            value=_measured(eclipses.fit(date_ingress).sunset),
        ),
        # how much brighter than the fit the plateau is toward the right of the
        # detector, per solar radius, in the last eclipse
        aastex.Variable(
            name="plateauGradient",
            value=_percent(
                (np.hypot(gradient.x, gradient.y) / last.center).to_value(u.one)
            ),
        ),
        # the fraction of the signal at disk center that is leak, in the last eclipse
        aastex.Variable(
            name="leakFractionQuietSun",
            value=aastex.NoEscape(rf"{_figures(fraction)}\%"),
        ),
        # the noise of the run on the day of the ingress shown
        aastex.Variable(name="noiseFrames", value=aia_94_leak_paper.noise.frames),
        aastex.Variable(name="readNoise", value=_measured(np.nanmean(noise.read))),
        aastex.Variable(
            name="noiseSlope", value=_measured(np.nanmean(noise.slope) * u.DN)
        ),
        aastex.Variable(name="noiselessCenter", value=_measured(noise.noiseless)),
        aastex.Variable(
            name="noiselessSpread", value=_measured(noise.spread, digits=1)
        ),
        # the map of the noiseless signal on that day: its limb darkening, how
        # much it grows per DN of the 94 Å signal of a region, and the part of
        # that signal it is at disk center
        aastex.Variable(
            name="noiseLimbDarkening",
            value=aastex.NoEscape(rf"\ensuremath{{{_figures(disk.limb_darkening)}}}"),
        ),
        aastex.Variable(
            name="noisePerSignal",
            value=aastex.NoEscape(
                rf"\ensuremath{{{_figures(disk.per_signal, digits=1)}}}"
            ),
        ),
        aastex.Variable(
            name="noiseShare",
            value=aastex.NoEscape(
                rf"{_figures((disk.amplitude / disk.signal).to_value(u.percent))}\%"
            ),
        ),
        # the typical spread of the noise method among the quadrants, and how
        # its noiseless signal compares with the eclipses once the leak is
        # three times that
        aastex.Variable(name="noiseSpread", value=_measured(spread, digits=1)),
        aastex.Variable(name="yearNoise", value=min(k.date.year for _, k in measured)),
        aastex.Variable(
            name="noiseRatioMin",
            value=aastex.NoEscape(rf"{_figures(min(ratios))}\%"),
        ),
        aastex.Variable(
            name="noiseRatioMax",
            value=aastex.NoEscape(rf"{_figures(max(ratios))}\%"),
        ),
        # the day the channel observed through its thick focal-plane filter,
        # how far apart the runs through each filter began, how much more
        # zirconium the thick filter has, and the band they are compared in
        aastex.Variable(
            name="filterDate",
            value=aastex.NoEscape(
                f"{filters.switch:%Y} {filters.switch:%B} {filters.switch.day}"
            ),
        ),
        aastex.Variable(name="filterGap", value=_measured(gap.to(u.min))),
        aastex.Variable(
            name="filterZirconium",
            value=filters.zirconium["thick"] - filters.zirconium["thin"],
        ),
        aastex.Variable(
            name="filterBand",
            value=aastex.NoEscape(rf"\ensuremath{{{filters.band:g}\,R_\odot}}"),
        ),
        # how much the thick filter passes, relative to the thin one, of the
        # EUV (from the structure of the corona, and as the calibration
        # expects) and of the leak (from the noise, with its spread among the
        # quadrants, and from the signals alone)
        aastex.Variable(
            name="filterEUV",
            value=aastex.NoEscape(rf"\ensuremath{{{_figures(attenuation.euv)}}}"),
        ),
        aastex.Variable(
            name="filterTable",
            value=aastex.NoEscape(rf"\ensuremath{{{_figures(filters.euv_table)}}}"),
        ),
        aastex.Variable(
            name="filterLeak",
            value=aastex.NoEscape(rf"\ensuremath{{{_figures(attenuation.leak)}}}"),
        ),
        aastex.Variable(
            name="filterLeakSpread",
            value=aastex.NoEscape(
                rf"\ensuremath{{{_figures(attenuation.leak_spread, digits=1)}}}"
            ),
        ),
        aastex.Variable(
            name="filterLeakSignals",
            value=aastex.NoEscape(
                rf"\ensuremath{{{_figures(attenuation.leak_signals)}}}"
            ),
        ),
        # the images without a focal-plane filter: the month of the first, the
        # ring they are compared with the leak in, and how much they have
        # brightened there since the first
        aastex.Variable(
            name="entranceFirst",
            value=aastex.NoEscape(f"{monitor.time[0].datetime:%Y %B}"),
        ),
        aastex.Variable(
            name="entranceRing",
            value=aastex.NoEscape(
                rf"{entrance.ring[0]:g}--{entrance.ring[1]:g}\,\ensuremath{{R_\odot}}"
            ),
        ),
        aastex.Variable(
            name="entranceGrowth",
            value=aastex.NoEscape(_figures((bright[-1] / bright[0]).to_value(u.one))),
        ),
        # the leak against them: since when it dominates the glow before the
        # onset, its mean ratio to them since then and the scatter of that,
        # how much the leak grew meanwhile, and how much higher the ratio was
        # in the first eclipse with an image
        aastex.Variable(name="entranceDominant", value=entrance.dominance),
        aastex.Variable(
            name="entranceSince", value=comparison.dates[first_dominant].year
        ),
        aastex.Variable(
            name="entranceRatio",
            value=aastex.NoEscape(
                rf"\ensuremath{{{_figures(comparison.ratio_mean * 1e7)}"
                r" \times 10^{-7}}"
            ),
        ),
        aastex.Variable(name="entranceScatter", value=_percent(comparison.scatter)),
        aastex.Variable(
            name="entranceGrowthSince",
            value=aastex.NoEscape(
                _figures(
                    (comparison.leak[-1] / comparison.leak[first_dominant]).to_value(
                        u.one
                    )
                )
            ),
        ),
        aastex.Variable(
            name="entranceRatioFirst",
            value=aastex.NoEscape(
                _figures(comparison.ratio[0] / comparison.ratio_mean)
            ),
        ),
        aastex.Variable(name="entranceRatioFirstYear", value=comparison.dates[0].year),
        # the day whose DEMs are shown, the temperature above which plasma
        # counts as hot, and how much of the emission measure of its coronal
        # holes and quiet Sun is hot, as observed and with the leak subtracted
        aastex.Variable(
            name="demDate",
            value=aastex.NoEscape(f"{date_dem:%Y} {date_dem:%B} {date_dem.day}"),
        ),
        aastex.Variable(
            name="hotTemperature",
            value=_measured((10**dems.hot * u.K).to(u.MK)),
        ),
        aastex.Variable(name="hotQuietObserved", value=_percent(shown.quiet_observed)),
        aastex.Variable(
            name="hotQuietCorrected", value=_percent(shown.quiet_corrected)
        ),
        aastex.Variable(name="hotHolesObserved", value=_percent(shown.holes_observed)),
        aastex.Variable(
            name="hotHolesCorrected", value=_percent(shown.holes_corrected)
        ),
        # the most of it before the onset as observed, and the most of it on any
        # day with the leak subtracted
        aastex.Variable(
            name="hotQuietEarly",
            value=_percent(max(s.quiet_observed for s in early)),
        ),
        aastex.Variable(
            name="hotHolesEarly",
            value=_percent(
                max(s.holes_observed for s in early if np.isfinite(s.holes_observed))
            ),
        ),
        aastex.Variable(
            name="hotQuietCorrectedMax",
            value=_percent(max(s.quiet_corrected for s in summaries)),
        ),
        aastex.Variable(
            name="hotHolesCorrectedMax",
            value=_percent(
                max(
                    s.holes_corrected
                    for s in summaries
                    if np.isfinite(s.holes_corrected)
                )
            ),
        ),
        aastex.Variable(name="holesFewest", value=dems.fewest),
        # the side of the bins the DEMs are inverted in
        aastex.Variable(name="demBin", value=_measured(dems.load(date_dem).scale)),
        # how bright the quiet Sun is in 94 Å with the leak subtracted, against
        # what its DEM predicts, in the first and the last eclipse
        aastex.Variable(name="ratioFirst", value=_percent(summaries[0].ratio_94)),
        aastex.Variable(name="ratioLast", value=_percent(summaries[-1].ratio_94)),
    ]
