import aastex
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import aia_94_leak_paper

__all__ = [
    "entrance",
]


def entrance() -> aastex.FigureStar:
    """
    The brightness of the images without a focal-plane filter against the
    leak in the eclipses, near the limb, and the ratio of the two.
    """
    module = aia_94_leak_paper.entrance
    monitor = module.load()
    brightness = module.brightness(monitor)
    comparison = module.compare()
    rate = u.DN / u.s
    year = comparison.time.decimalyear
    dominant = comparison.dominant

    fig, (left, right) = plt.subplots(
        ncols=2,
        figsize=(aastex.text_width_inches, 2.5),
        constrained_layout=True,
    )

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
\entranceFirst, against the leak, both in the ring at \entranceRing, in the units of the detector.
Left: the leak in each eclipse (red), through the point spread function fit to the eclipses, and
the brightness of the images (gray, right axis), scaled by the mean ratio of the two.
Right: the ratio of the leak to the brightness of the images at the time of each eclipse, with its
mean and scatter (gray) over the eclipses since the leak in the ring was \entranceDominant\ times
its mean before \yearOnset\ (filled).""")
    )
    return figure
