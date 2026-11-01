import itertools

import aastex
import astropy.units as u
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np

import aia_94_leak_paper

__all__ = [
    "filters",
]


def filters() -> aastex.FigureStar:
    """
    The signal of the last frame through the thick focal-plane filter
    relative to the run through the thin one, and how much of the leak the
    thick filter passes.
    """
    module = aia_94_leak_paper.filters
    noise = aia_94_leak_paper.noise
    thin = module.run()
    thick = module.frame("clean")
    axes = ("y", "x")
    r = thin.r.ndarray_aligned(axes)
    position = thin.position
    x, y = position.x.ndarray, position.y.ndarray
    dx, dy = x[1] - x[0], y[1] - y[0]
    extent = [x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2]

    fig, (left, right) = plt.subplots(
        ncols=2,
        figsize=(aastex.text_width_inches, 2.4),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1, 1.2]},
    )

    # the frame through the thick filter over the run through the thin one,
    # where there is signal enough to see it
    signal_thick = thick.image.to_value(u.DN).ndarray_aligned(axes)
    signal_thin = thin.image.to_value(u.DN).ndarray_aligned(axes)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(r < 1.25, signal_thick / signal_thin, np.nan)
    im = left.imshow(
        ratio,
        origin="lower",
        extent=extent,
        cmap="viridis",
        vmin=0,
        vmax=1,
        interpolation="none",
    )
    left.add_patch(
        matplotlib.patches.Circle(
            (0, 0), 1, fill=False, color="white", linewidth=0.5, linestyle="dotted"
        )
    )
    left.set_xlim(-1.3, 1.3)
    left.set_ylim(-1.3, 1.3)
    left.set_xlabel(r"$x$ ($R_\odot$)")
    left.set_ylabel(r"$y$ ($R_\odot$)")
    bar = fig.colorbar(im, ax=left, shrink=0.85)
    bar.set_label("thick / thin")

    # how much of the leak the thick filter passes in each region: the EUV
    # through it, plus what is left of the frame through it once that EUV is
    # taken away, over the leak through the thin filter, by each measure of
    # that leak, with the EUV fit with it; the median of each ring, and the
    # spread of the medians among the quadrants
    quadrants = [quadrant.ndarray_aligned(axes) for quadrant in noise.quadrants(thin)]
    styles = {
        "eclipses": {"color": "tab:blue", "label": "leak from the eclipses"},
        "noise": {"color": "tab:gray", "mfc": "white", "label": "leak from the noise"},
    }
    offset = {"eclipses": 0.008, "noise": -0.008}
    edges = np.arange(0, module.radius + 1e-9, 0.1)
    right.axhline(1, color="black", linestyle="dashed", linewidth=0.8, label="passed")
    right.axhline(0, color="black", linewidth=0.8, label="stopped")
    for measure in module.measures:
        tau = module.attenuation(measure).euv
        leak = module.leak_thin(measure).to_value(u.DN).ndarray_aligned(axes)
        # the leak from the noise is nothing beyond the limb
        with np.errstate(divide="ignore", invalid="ignore"):
            passed = tau + (signal_thick - tau * signal_thin) / leak
        median, spread = [], []
        for lo, hi in itertools.pairwise(edges):
            ring = (r >= lo) & (r < hi)
            median.append(np.nanmedian(passed[ring]))
            spread.append(np.std([np.nanmedian(passed[ring & q]) for q in quadrants]))
        right.errorbar(
            (edges[:-1] + edges[1:]) / 2 + offset[measure],
            median,
            spread,
            fmt="o",
            markersize=3,
            elinewidth=0.6,
            **styles[measure],
        )
    right.set_xlim(0, module.radius)
    right.set_ylim(-0.4, 1.2)
    right.set_xlabel(r"radius ($R_\odot$)")
    right.set_ylabel("leak, thick / thin")
    right.legend(
        fontsize=7,
        loc="center right",
        handlelength=1.5,
        labelspacing=0.3,
        frameon=False,
    )

    figure = aastex.FigureStar("fig:filters", position="!htb")
    figure.append(aastex.NoEscape(r"\centering"))
    figure.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    figure.add_caption(
        aastex.NoEscape(r"""
The two focal-plane filters of 94\,\AA\ on \filterDate: the last frame through the thick filter,
and the first \filterFrames\ frames through the thin one, which began \filterGap\ later.
Left: the frame through the thick filter over the frames through the thin one, in each region of
$64 \times 64$ pixels.
Right: how much of the leak the thick filter passes relative to the thin one, within
\filterRadius\ of disk center: $\tau$, what it passes of the \EUV, plus the signal through the
thick filter less $\tau$ times that through the thin one, over the leak through the thin filter, as
the eclipses that bracket the day measure it (filled; Section~\ref{sec:EarthEclipses}) and as the
noise of the frames through the thin filter does (open; Section~\ref{sec:PhotonNoise}), with
$\tau$ fit with each.
The median of each ring of $0.1\,R_\odot$ (points), the spread of the medians among the quadrants
of the detector (bars), and what it would be if the thick filter stopped the leak (solid) or passed
all of it (dashed).""")
    )
    return figure
