import datetime

import aastex
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import aia_94_leak_paper

__all__ = [
    "date_ingress",
    "ingress",
    "times_ingress",
]

date_ingress: datetime.date = datetime.date(2026, 7, 15)
"""The eclipse shown, the latest, in which the leak is brightest."""

times_ingress: u.Quantity = [-20, 60, 110, 150, 185, 215] * u.s
"""When the frames shown were observed, after the 171 Å disk mean began to fall."""


def ingress() -> aastex.FigureStar:
    """
    Frames of 94 Å and 171 Å through the ingress of one eclipse, each less
    the umbra and divided by its exposure time.
    """
    eclipse = aia_94_leak_paper.eclipses.load(date_ingress)
    position = eclipse.position
    x = position.x.ndarray
    y = position.y.ndarray
    extent = [x[0], x[-1], y[0], y[-1]]

    fig, axs = plt.subplots(
        nrows=2,
        ncols=len(times_ingress),
        figsize=(aastex.text_width_inches, 2.75),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    rows = [
        # the leak is a few tenths of a DN per second, the corona of 171 a few
        # thousand, so one is shown linearly and the other logarithmically
        (94, {"cmap": "inferno", "vmin": 0, "vmax": 1.2}, lambda a: a),
        (
            171,
            {"cmap": "gray", "vmin": -1, "vmax": 3.3},
            lambda a: np.log10(np.maximum(a, 0.01)),
        ),
    ]
    for row, (channel, style, scale) in enumerate(rows):
        frames = (
            eclipse.images[f"INGRESS{channel}"] / eclipse.exposures[f"INGRESS{channel}"]
        )
        dark = np.median(
            eclipse.images[f"UMBRA{channel}"] / eclipse.exposures[f"UMBRA{channel}"],
            axis="time",
        )
        rate = (frames - dark).to(u.DN / u.s)
        if channel == 94:
            # as the leak is fit, in the units of the detector at disk center
            rate = rate * eclipse.flat_center
        seconds = (eclipse.times[f"INGRESS{channel}"].ndarray - eclipse.start).to(u.s)
        for col, t in enumerate(times_ingress):
            k = int(np.argmin(np.abs(seconds - t)))
            image = rate[{"time": k}].ndarray_aligned(("y", "x")).value
            ax = axs[row, col]
            im = ax.imshow(
                scale(image),
                origin="lower",
                extent=extent,
                interpolation="none",
                **style,
            )
            if row == 0:
                ax.set_title(f"${seconds[k].value:+.0f}$\\,s", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(rf"{channel}\,\AA", fontsize=8)
        label = r"DN\,s$^{-1}$" if channel == 94 else r"$\log_{10}$ DN\,s$^{-1}$"
        bar = fig.colorbar(im, ax=axs[row, :], shrink=0.9, pad=0.01)
        bar.set_label(label, fontsize=7)
        bar.ax.tick_params(labelsize=7)

    result = aastex.FigureStar("fig:ingress", position="!htb")
    result.append(aastex.NoEscape(r"\centering"))
    result.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    result.add_caption(
        aastex.NoEscape(r"""
The ingress of the eclipse of \ingressDate\ in the 94\,\AA\ (top) and 171\,\AA\ (bottom) channels,
at the times marked after the 171\,\AA\ disk mean began to fall.
Each frame is level~1, averaged over $16 \times 16$ pixels, less the median of the umbra just after
the ingress, and divided by its
exposure time; 94\,\AA\ is multiplied by the flat field at disk center, so that it is in the units of
the detector there.
The atmosphere absorbs the \EUV\ from the lower right, and the corona of both channels goes with it,
but in 94\,\AA\ a smooth, limb-darkened disk is left behind.
It sets \ingressSunset\ after the \EUV\ edge has passed, flattened by refraction, as the visible
Sun does.""")
    )
    return result
