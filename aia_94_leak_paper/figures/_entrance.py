import aastex
import astropy.units as u
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import aia_94_leak_paper

__all__ = [
    "entrance",
]

_spacing = 5
"""The years between the examples of the images, the last of which is the last image."""


def entrance() -> aastex.FigureStar:
    """
    Examples of the images without a focal-plane filter, their brightness
    against the leak in the eclipses near the limb, and the ratio of the two.
    """
    module = aia_94_leak_paper.entrance
    monitor = module.load()
    brightness = module.brightness(monitor)
    comparison = module.compare()
    rate = u.DN / u.s
    year = comparison.time.decimalyear
    dominant = comparison.dominant

    fig = plt.figure(
        figsize=(aastex.text_width_inches, 4.7),
        constrained_layout=True,
    )
    top, bottom = fig.subfigures(nrows=2, height_ratios=[2.1, 2.6])
    examples = top.subplots(ncols=4, sharex=True, sharey=True)
    left, right = bottom.subplots(ncols=2)

    # examples of the images, each divided by its brightness in the ring
    taken = monitor.time.decimalyear
    targets = taken[0] + _spacing * np.arange(len(examples) - 1)
    index = [int(np.argmin(np.abs(taken - t))) for t in targets] + [taken.size - 1]
    images = monitor.images.to_value(rate).ndarray_aligned(("time", "y", "x"))
    position = monitor.position
    for ax, i in zip(examples, index):
        x = position.x[{"time": i}].ndarray
        y = position.y[{"time": i}].ndarray
        im = ax.imshow(
            images[i] / brightness[i].to_value(rate),
            origin="lower",
            extent=[x[0], x[-1], y[0], y[-1]],
            cmap="gray",
            vmin=0,
            vmax=2,
            interpolation="none",
        )
        for radius in module.ring:
            ax.add_patch(
                matplotlib.patches.Circle(
                    (0, 0),
                    radius,
                    fill=False,
                    color="white",
                    linewidth=0.5,
                    linestyle="dashed",
                )
            )
        day = monitor.time[i].datetime
        ax.set_title(f"{day.year} {day:%B} {day.day}", fontsize=8)
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])
        ax.set_xlabel(r"$x$ ($R_\odot$)")
    examples[0].set_ylabel(r"$y$ ($R_\odot$)")
    bar = top.colorbar(im, ax=examples, shrink=0.9)
    bar.set_label("image / brightness")

    # the images on the right axis, scaled by the mean ratio, so that the
    # leak on the left coincides with them where the ratio holds
    ratio = comparison.ratio_mean
    left.plot(
        monitor.time.decimalyear,
        ratio * brightness.to_value(rate),
        "-",
        color="0.5",
        linewidth=0.8,
        label="without a focal-plane filter",
    )
    leak = comparison.leak.to_value(rate)
    for where, face, label in (
        (dominant, "tab:red", "leak, eclipses"),
        (~dominant, "none", None),
    ):
        left.plot(
            year[where],
            leak[where],
            "o",
            color="tab:red",
            markerfacecolor=face,
            markersize=4,
            label=label,
        )
    left.set_yscale("log")
    left.set_xlabel("year")
    left.set_ylabel(r"leak (DN\,s$^{-1}$)")
    left.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(4))
    twin = left.twinx()
    twin.set_yscale("log")
    low, high = left.get_ylim()
    twin.set_ylim(low / ratio, high / ratio)
    twin.set_ylabel(r"without filter (DN\,s$^{-1}$)")
    left.legend(fontsize=7, loc="upper left")

    # their ratio, and its mean and scatter where the leak dominates
    scale = 1e-7
    for where, face in ((dominant, "tab:red"), (~dominant, "none")):
        right.plot(
            year[where],
            comparison.ratio[where] / scale,
            "o",
            color="tab:red",
            markerfacecolor=face,
            markersize=4,
        )
    span = ratio * np.exp(comparison.scatter * np.array([-1, 1])) / scale
    right.axhspan(*span, color="0.85", zorder=0)
    right.axhline(ratio / scale, color="0.4", linewidth=0.8, zorder=0)
    right.set_yscale("log")
    right.set_xlim(left.get_xlim())
    right.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(4))
    right.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    right.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    right.set_yticks([2, 3, 5, 10, 20])
    right.set_xlabel("year")
    right.set_ylabel(r"leak / without filter ($10^{-7}$)")

    figure = aastex.FigureStar("fig:entrance", position="!htb")
    figure.append(aastex.NoEscape(r"\centering"))
    figure.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    figure.add_caption(
        aastex.NoEscape(r"""
The images of the 94 and 304\,\AA\ telescope taken without a focal-plane filter, one a month since
\entranceFirst, against the leak.
Top: four of the images, from the first to the last, in bins of $16 \times 16$ pixels, each divided
by its brightness, the median of its bins in the ring at \entranceRing\ (dashed).
Bottom left: the leak in each eclipse (red), through the point spread function fit to the eclipses,
and the brightness of the images (gray, right axis), scaled by the mean ratio of the two, both in
the ring and in the units of the detector.
Bottom right: the ratio of the leak to the brightness of the images at the time of each eclipse,
with its mean and scatter (gray) over the eclipses since the leak in the ring was
\entranceDominant\ times its mean before \yearOnset\ (filled).""")
    )
    return figure
