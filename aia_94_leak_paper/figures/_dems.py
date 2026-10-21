import datetime

import aastex
import astropy.units as u
import matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import aia_94_leak_paper

__all__ = [
    "dems",
]


def _year(date: datetime.date) -> float:
    """The decimal year of noon on a day."""
    return date.year + (date.timetuple().tm_yday - 0.5) / 365.25


def dems() -> aastex.FigureStar:
    """
    The fraction of the emission measure above 3 MK on one day, as observed
    and with the leak subtracted, the DEMs of its coronal holes and quiet
    Sun, and the history of both over the mission.
    """
    module = aia_94_leak_paper.dems
    date = module.date_dem
    synoptic = module.load(date)
    chosen = module.regions(date)
    observed = module.invert(date, corrected=False)
    corrected = module.invert(date, corrected=True)

    fig = plt.figure(figsize=(aastex.text_width_inches, 4.6), constrained_layout=True)
    outer = fig.add_gridspec(2, 1, height_ratios=[1.15, 1])
    top = outer[0].subgridspec(1, 4, width_ratios=[1, 1, 0.06, 1.2])
    left = fig.add_subplot(top[0])
    middle = fig.add_subplot(top[1])
    cax = fig.add_subplot(top[2])
    right = fig.add_subplot(top[3])
    bottom = fig.add_subplot(outer[1])

    # the fraction above 3 MK on the disk, where the fit succeeded
    axes = ("y", "x")
    r = synoptic.r.ndarray_aligned(axes)
    position = synoptic.position
    x, y = position.x.ndarray, position.y.ndarray
    holes = chosen.holes.ndarray_aligned(axes)
    for ax, dem, title in (
        (left, observed, "as observed"),
        (middle, corrected, "leak subtracted"),
    ):
        fraction = dem.fraction_hot.ndarray_aligned(axes)
        good = dem.chi2.ndarray_aligned(axes) >= 0
        image = np.where((r < 1) & good, fraction, np.nan)
        im = ax.imshow(
            image,
            origin="lower",
            extent=[x[0], x[-1], y[0], y[-1]],
            cmap="inferno",
            vmin=0,
            vmax=1,
            interpolation="none",
        )
        ax.contour(x, y, holes, levels=[0.5], colors="white", linewidths=0.5, alpha=0.8)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(r"$x$ ($R_\odot$)")
    left.set_ylabel(r"$y$ ($R_\odot$)")
    middle.tick_params(labelleft=False)
    bar = fig.colorbar(im, cax=cax)
    label = f"above {10**module.hot / 1e6:.1f} MK"
    bar.set_label(label)
    bar.ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))

    # the median DEM of the coronal holes and the quiet Sun
    logt = np.log10(observed.temperature.ndarray.to_value(u.K))
    unit = 1e26 / u.cm**5
    colors = {"coronal holes": "tab:purple", "quiet Sun": "tab:green"}
    for name, where in (("coronal holes", chosen.holes), ("quiet Sun", chosen.quiet)):
        for dem, style in ((observed, "dashed"), (corrected, "solid")):
            use = (where & (dem.chi2 >= 0)).ndarray_aligned(axes)
            values = dem.values.to_value(unit).ndarray_aligned(axes + ("temperature",))
            right.plot(
                logt,
                np.median(values[use], axis=0),
                color=colors[name],
                linestyle=style,
                linewidth=1,
            )
    right.axvline(module.hot, color="gray", linestyle="dotted", linewidth=0.8)
    right.set_xlim(5.6, 7.3)
    # headroom above the curves for the key to the styles
    right.set_ylim(0, 1.4 * right.get_ylim()[1])
    right.set_xlabel(r"$\log_{10} T$ (K)")
    right.set_ylabel(r"DEM ($10^{26}\,\mathrm{cm^{-5}}$)")
    # the regions are labeled on their curves, and the styles in the corner
    # the curves leave empty
    for name, (x0, y0, align) in (
        ("quiet Sun", (6.2, 0.6, "right")),
        ("coronal holes", (6.78, 0.44, "center")),
    ):
        right.text(
            x0,
            y0,
            name,
            color=colors[name],
            fontsize=7,
            ha=align,
            transform=right.get_xaxis_transform(),
        )
    handles = [
        matplotlib.lines.Line2D(
            [], [], color="black", linestyle=s, linewidth=1, label=n
        )
        for n, s in (("as observed", "dashed"), ("leak subtracted", "solid"))
    ]
    right.legend(handles=handles, fontsize=7, loc="upper right", handlelength=1.5)

    # the history of both, as observed (open) and with the leak subtracted (filled)
    history = module.history()
    year = np.array([_year(s.date) for s in history])
    for name, observed_, corrected_ in (
        ("quiet Sun", "quiet_observed", "quiet_corrected"),
        ("coronal holes", "holes_observed", "holes_corrected"),
    ):
        color = colors[name]
        bottom.plot(
            year,
            [getattr(s, observed_) for s in history],
            "o",
            color=color,
            markerfacecolor="none",
            markersize=4,
        )
        bottom.plot(
            year,
            [getattr(s, corrected_) for s in history],
            "o",
            color=color,
            markersize=4,
            label=name,
        )
    bottom.set_ylim(0, 0.8)
    bottom.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    bottom.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    bottom.set_xlabel("year")
    bottom.set_ylabel(label)
    handles, _ = bottom.get_legend_handles_labels()
    handles += [
        matplotlib.lines.Line2D(
            [], [], color="gray", marker="o", linestyle="", markersize=4, label=n, **k
        )
        for n, k in (
            ("as observed", {"markerfacecolor": "none"}),
            ("leak subtracted", {}),
        )
    ]
    bottom.legend(handles=handles, fontsize=7, loc="upper left", ncols=2)

    result = aastex.FigureStar("fig:dems", position="!htb")
    result.append(aastex.NoEscape(r"\centering"))
    result.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    result.add_caption(
        aastex.NoEscape(r"""
The \DEM\ of each bin, \demBin\ on a side, inverted from the six coronal channels at noon of each
eclipse day, as observed and with the leak measured in the eclipse of the same day subtracted
from 94\,\AA.
Top left and middle: the fraction of the emission measure above \hotTemperature\ on \demDate, as
observed and with the leak subtracted, on the disk where the inversion succeeded; the coronal holes
are outlined in white.
Top right: the median \DEM\ of the coronal holes and of the quiet Sun that day, as observed (dashed)
and with the leak subtracted (solid); the dotted line is \hotTemperature.
Bottom: the median fraction above \hotTemperature\ in the coronal holes and the quiet Sun of every
eclipse day, as observed (open) and with the leak subtracted (filled), with the coronal holes shown
where they cover at least \holesFewest\ bins.""")
    )
    return result
