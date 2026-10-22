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
    The signal through the thick focal-plane filter relative to the thin
    one, and the noiseless signal through each in the band where the thick
    filter transmits evenly.
    """
    module = aia_94_leak_paper.filters
    noise = aia_94_leak_paper.noise
    runs = {name: module.load(name) for name in module.filters}
    thin = runs["thin"]
    result = module.attenuation()
    axes = ("y", "x")
    r = thin.r.ndarray_aligned(axes)
    use = module.within(thin)
    position = thin.position
    x, y = position.x.ndarray, position.y.ndarray

    fig, (left, right) = plt.subplots(
        ncols=2,
        figsize=(aastex.text_width_inches, 2.6),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1, 1.25]},
    )

    # the signal through the thick filter over that through the thin one,
    # where there is signal enough to see it
    signal = {
        n: run.image.to_value(u.DN).ndarray_aligned(axes) for n, run in runs.items()
    }
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(r < 1.25, signal["thick"] / signal["thin"], np.nan)
    im = left.imshow(
        ratio,
        origin="lower",
        extent=[x[0], x[-1], y[0], y[-1]],
        cmap="viridis",
        vmin=0.3,
        vmax=1,
        interpolation="none",
    )
    left.contour(x, y, use, levels=[0.5], colors="white", linewidths=0.8)
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

    # the noiseless signal through each filter, the median of each ring of
    # the band, and what the thick filter would pass of it if it attenuated
    # the leak as it does the EUV
    edges = np.arange(0, 0.9001, 0.1)
    middle = (edges[:-1] + edges[1:]) / 2
    rings = {}
    for name, run in runs.items():
        p = noise.noiseless(run).to_value(u.DN).ndarray_aligned(axes)
        rings[name] = np.array(
            [
                np.nanmedian(p[use & (r >= lo) & (r < hi)])
                for lo, hi in itertools.pairwise(edges)
            ]
        )
    colors = {"thin": "tab:red", "thick": "tab:blue"}
    for name in module.filters:
        right.plot(
            middle,
            rings[name],
            "o",
            color=colors[name],
            markersize=4,
            label=f"{name} filter",
        )
    right.plot(
        middle,
        result.euv * rings["thin"],
        color="tab:blue",
        linestyle="dashed",
        linewidth=0.8,
        label="thin, attenuated as the EUV",
    )
    right.set_xlim(0, 0.9)
    right.set_ylim(0, 2.5)
    right.set_xlabel(r"radius ($R_\odot$)")
    right.set_ylabel("$p$ (DN)")
    right.legend(fontsize=7, loc="lower left", handlelength=1.5, labelspacing=0.3)

    figure = aastex.FigureStar("fig:filters", position="!htb")
    figure.append(aastex.NoEscape(r"\centering"))
    figure.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    figure.add_caption(
        aastex.NoEscape(r"""
The two focal-plane filters of 94\,\AA\ on \filterDate: the last run of \noiseFrames\ frames through
the thick filter and the first through the thin one, which began \filterGap\ apart.
Left: the signal through the thick filter over that through the thin one, in each region of
$64 \times 64$ pixels; the thick filter passes less toward the top and bottom of the field, and the
two are compared within the white outline, \filterBand\ to either side of the middle of the disk and
within $0.9\,R_\odot$ of its center.
Right: the noiseless signal $p$ through each filter, the median of each ring of $0.1\,R_\odot$
within the outline, and that through the thin filter times \filterEUV, the attenuation of the \EUV\
measured from the structure of the corona (dashed).""")
    )
    return figure
