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
    The signal of each frame through the thick focal-plane filter relative
    to the run through the thin one, and what the thick filter passes of the
    leak.
    """
    module = aia_94_leak_paper.filters
    noise = aia_94_leak_paper.noise
    thin = module.run()
    frames = {name: module.frame(name) for name in module.frames}
    result = module.attenuation()
    disk = noise.disk(thin)
    axes = ("y", "x")
    r = thin.r.ndarray_aligned(axes)
    position = thin.position
    x, y = position.x.ndarray, position.y.ndarray
    dx, dy = x[1] - x[0], y[1] - y[0]
    extent = [x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2]

    fig, (banded, clean, right) = plt.subplots(
        ncols=3,
        figsize=(aastex.text_width_inches, 2.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1, 1, 1.15]},
    )

    # each frame through the thick filter over the run through the thin one,
    # where there is signal enough to see it
    signal = {
        name: frame.image.to_value(u.DN).ndarray_aligned(axes)
        for name, frame in frames.items()
    }
    signal_thin = thin.image.to_value(u.DN).ndarray_aligned(axes)
    for ax, name in ((banded, "banded"), (clean, "clean")):
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(r < 1.25, signal[name] / signal_thin, np.nan)
        im = ax.imshow(
            ratio,
            origin="lower",
            extent=extent,
            cmap="viridis",
            vmin=0,
            vmax=1,
            interpolation="none",
        )
        ax.add_patch(
            matplotlib.patches.Circle(
                (0, 0), 1, fill=False, color="white", linewidth=0.5, linestyle="dotted"
            )
        )
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_xlabel(r"$x$ ($R_\odot$)")
    banded.set_ylabel(r"$y$ ($R_\odot$)")
    clean.tick_params(labelleft=False)
    bar = fig.colorbar(im, ax=[banded, clean], shrink=0.85)
    bar.set_label("thick / thin")

    # the last frame through the thick filter less the EUV through it, the
    # median of each ring of the disk, and what it would be if the thick
    # filter stopped the leak or attenuated it as it does the EUV
    edges = np.arange(0, module.radius + 1e-9, 0.1)
    middle = (edges[:-1] + edges[1:]) / 2
    excess = signal["clean"] - result.euv * signal_thin
    rings = [
        np.nanmedian(excess[(r >= lo) & (r < hi)])
        for lo, hi in itertools.pairwise(edges)
    ]
    right.plot(middle, rings, "o", color="tab:blue", markersize=4, label="measured")
    radius = np.linspace(0, module.radius, 101)
    stopped = -result.euv * disk(radius).to_value(u.DN)
    right.plot(radius, stopped, color="black", linewidth=0.8, label="leak stopped")
    right.axhline(
        0,
        color="black",
        linestyle="dashed",
        linewidth=0.8,
        label="leak attenuated as the EUV",
    )
    right.set_xlim(0, module.radius)
    right.set_ylim(-2, 1)
    right.set_xlabel(r"radius ($R_\odot$)")
    right.set_ylabel(r"$S_\mathrm{thick} - \tau S_\mathrm{thin}$ (DN)")
    right.legend(fontsize=7, loc="upper right", handlelength=1.5, labelspacing=0.3)

    figure = aastex.FigureStar("fig:filters", position="!htb")
    figure.append(aastex.NoEscape(r"\centering"))
    figure.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    figure.add_caption(
        aastex.NoEscape(r"""
The two focal-plane filters of 94\,\AA\ on \filterDate: the last two frames through the thick
filter, as the channel returned to the thin one, and the first \filterFrames\ frames through the
thin filter, which began \filterGap\ after the last.
Left: the first of the two frames through the thick filter, read out as the filter wheel moved to
the thick filter of 304\,\AA, over the frames through the thin one, in each region of
$64 \times 64$ pixels; it carries a band across the middle rows of the detector.
Middle: the same for the last frame, read out as the wheel moved to the thin filter of 304\,\AA.
Right: the signal of the last frame through the thick filter less $\tau = \filterEUV$ times that
through the thin one, the median of each ring of $0.1\,R_\odot$ (points), and what it would be if
the thick filter stopped the leak (solid) or attenuated it as much as the \EUV\ (dashed).""")
    )
    return figure
