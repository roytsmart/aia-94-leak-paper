import itertools

import aastex
import astropy.units as u
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np

import aia_94_leak_paper

__all__ = [
    "noise",
]


def noise() -> aastex.FigureStar:
    """
    The variance of 94 Å against its signal near disk center and off the
    limb in one run, and the map of the noiseless signal and its radial
    profile.
    """
    date = aia_94_leak_paper.figures.date_ingress
    module = aia_94_leak_paper.noise
    run = module.load(date)
    result = module.fit(run)
    r = run.r

    fig, (left, middle, right) = plt.subplots(
        ncols=3,
        figsize=(aastex.text_width_inches, 2.4),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.2, 1, 1]},
    )

    # variance against signal, near disk center and off the limb
    regions = [
        ("near disk center", r < 0.5, "tab:red"),
        (r"off the limb", (r > 1.05) & (r < 1.3), "tab:blue"),
    ]
    for label, where, color in regions:
        table = module.pooled(run, where)
        count = table["count"].ndarray
        signal = table["signal"].ndarray.to_value(u.DN)
        variance = table["variance"].ndarray.to_value(u.DN**2)
        use = (count >= 3000) & (signal > 0.3) & (signal < 30)
        left.plot(
            signal[use], variance[use], "o", color=color, markersize=3, label=label
        )
    # the line and the read noise shown are the means of the quadrants'
    s = np.linspace(0, 20, 50)
    read = np.nanmean(result.read).to_value(u.DN**2)
    left.plot(
        s,
        np.nanmean(result.intercept).to_value(u.DN**2) + np.nanmean(result.slope) * s,
        color="tab:red",
        linewidth=0.8,
    )
    left.axhline(
        read, color="black", linestyle="dashed", linewidth=0.8, label="read noise"
    )
    # what the variance would be if all of the signal carried photon noise:
    # the slope near disk center, started from the read noise at zero signal
    left.plot(
        s,
        read + np.nanmean(result.slope) * s,
        color="black",
        linestyle="dotted",
        linewidth=0.8,
        label="photon noise alone",
    )
    p = result.noiseless.to_value(u.DN)
    left.annotate(
        "",
        xy=(p, read),
        xytext=(0, read),
        arrowprops={
            "arrowstyle": "<->",
            "color": "tab:red",
            "linewidth": 0.8,
            "mutation_scale": 6,
            "shrinkA": 0,
            "shrinkB": 0,
        },
    )
    left.text(
        p / 2, read * 0.85, "$p$", color="tab:red", ha="center", va="top", fontsize=8
    )
    left.set_xlim(0, 20)
    left.set_ylim(0, 16)
    left.set_xlabel("signal (DN)")
    left.set_ylabel(r"variance (DN$^2$)")
    left.legend(fontsize=7, loc="upper left", handlelength=1.5, labelspacing=0.3)

    # the map of the noiseless signal on the disk; off the limb the slope of
    # the disk does not hold
    pmap = module.noiseless(run).to_value(u.DN)
    image = np.where(r < 1, pmap, np.nan).ndarray_aligned(("y", "x"))
    position = run.position
    x, y = position.x.ndarray, position.y.ndarray
    im = middle.imshow(
        image,
        origin="lower",
        extent=[x[0], x[-1], y[0], y[-1]],
        cmap="inferno",
        vmin=0,
        vmax=2.5,
        interpolation="none",
    )
    middle.add_patch(
        matplotlib.patches.Circle(
            (0, 0), 1, fill=False, color="white", linewidth=0.5, linestyle="dotted"
        )
    )
    middle.set_xlim(-1.3, 1.3)
    middle.set_ylim(-1.3, 1.3)
    middle.set_xlabel(r"$x$ ($R_\odot$)")
    middle.set_ylabel(r"$y$ ($R_\odot$)")
    bar = fig.colorbar(im, ax=middle, shrink=0.85)
    bar.set_label("$p$ (DN)")

    # its radial profile, each region and the median of each ring, against
    # the limb-darkened disk fit to the map
    radius = r.ndarray_aligned(("y", "x"))
    right.plot(
        radius.ravel(),
        image.ravel(),
        ".",
        color="0.7",
        markersize=1,
        label="regions",
    )
    edges = np.arange(0, 1.001, 0.05)
    rings = []
    for lo, hi in itertools.pairwise(edges):
        ring = image[(radius >= lo) & (radius < hi)]
        ring = ring[np.isfinite(ring)]
        rings.append(np.median(ring) if ring.size else np.nan)
    right.plot(
        (edges[:-1] + edges[1:]) / 2,
        rings,
        "o",
        color="tab:red",
        markersize=3,
        label="median of each ring",
    )
    fitted = module.disk(run)
    t = np.linspace(0, 1, 200)
    right.plot(
        t,
        fitted(t).to_value(u.DN),
        color="black",
        linewidth=0.8,
        label=r"photosphere, $1 - u (1 - \mu)$",
    )
    right.axhline(0, color="gray", linewidth=0.5)
    right.set_xlim(0, 1)
    right.set_ylim(-0.5, 3.5)
    right.set_xlabel(r"radius ($R_\odot$)")
    right.set_ylabel("$p$ (DN)")
    right.legend(fontsize=7, loc="lower left", handlelength=1.2, borderaxespad=0.3)

    result_figure = aastex.FigureStar("fig:noise", position="!htb")
    result_figure.append(aastex.NoEscape(r"\centering"))
    result_figure.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    result_figure.add_caption(
        aastex.NoEscape(r"""
The noise of 94\,\AA\ in a run of \noiseFrames\ frames at noon on \ingressDate, in the units of the
detector.
Left: the variance of one frame against its signal, near disk center (within half a radius) and off
the limb (1.05 to 1.3 radii), and the read noise from the corners beyond the field stop.
Off the limb the variance rises from the read noise; near disk center it reaches the read noise
only at a signal of $p = \noiselessCenter \pm \noiselessSpread$, which carries no photon noise.
The dotted line is the variance if all of the signal carried photon noise, the slope near disk
center started from the read noise at zero signal: the points off the limb follow it, and those near
disk center lie $p$ to its right.
Each quadrant of the detector, read out by its own amplifier, is fit with its own read noise and
slope; the line and the read noise shown are their means, and the uncertainty of $p$ is its spread
among the quadrants.
Middle: $p$ in each region of $64 \times 64$ pixels on the disk, with the read noise and the slope of
its quadrant; off the limb the slope of the disk does not hold.
Right: its radial profile, each region (gray) and the median of each ring of $0.05\,R_\odot$ (red),
with the limb darkening of a photosphere, $1 - u (1 - \mu)$, where $\mu = \sqrt{1 - r^2}$ is the cosine
of the angle between the line of sight and the vertical, fit to the regions within $0.95\,R_\odot$:
$u = \noiseLimbDarkening$.""")
    )
    return result_figure
