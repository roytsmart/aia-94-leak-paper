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
    limb in one run, and the map of the noiseless signal.
    """
    date = aia_94_leak_paper.figures.date_ingress
    module = aia_94_leak_paper.noise
    run = module.load(date)
    result = module.fit(date)
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
    s = np.linspace(0, 20, 50)
    read = result.read.to_value(u.DN**2)
    left.plot(
        s,
        result.intercept.to_value(u.DN**2) + result.slope * s,
        color="tab:red",
        linewidth=0.8,
    )
    left.axhline(
        read, color="black", linestyle="dashed", linewidth=0.8, label="read noise"
    )
    p = result.noiseless.to_value(u.DN)
    left.annotate(
        "",
        xy=(p, read),
        xytext=(0, read),
        arrowprops={"arrowstyle": "<->", "color": "tab:red", "linewidth": 0.8},
    )
    left.text(
        p / 2, read * 0.85, "$p$", color="tab:red", ha="center", va="top", fontsize=8
    )
    left.set_xlim(0, 20)
    left.set_ylim(0, 12)
    left.set_xlabel("signal (DN)")
    left.set_ylabel(r"variance (DN$^2$)")
    left.legend(fontsize=7, loc="upper left")

    # the map of the noiseless signal
    pmap = module.noiseless(date).to_value(u.DN)
    image = pmap.ndarray_aligned(("y", "x"))
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

    # its radial profile, against the leak measured in the eclipse of the same day,
    # both per second
    edges = np.arange(0, 1.41, 0.05)
    mids = (edges[:-1] + edges[1:]) / 2
    rate = u.DN / u.s

    def profile(values: np.ndarray, radius: np.ndarray) -> list[float]:
        return [
            np.nanmedian(values[(radius >= lo) & (radius < hi)])
            for lo, hi in itertools.pairwise(edges)
        ]

    exposure = result.exposure.to_value(u.s)
    radius = r.ndarray_aligned(("y", "x"))
    right.plot(
        mids,
        np.array(profile(image, radius)) / exposure,
        "o",
        color="tab:red",
        markersize=3,
        label="noiseless signal",
    )
    eclipse = aia_94_leak_paper.eclipses.load(date)
    position = eclipse.position
    radius = np.sqrt(np.square(position.x) + np.square(position.y))
    seen = aia_94_leak_paper.eclipses.seen(date).to_value(rate)
    right.plot(
        mids,
        profile(seen.ndarray_aligned(("y", "x")), radius.ndarray_aligned(("y", "x"))),
        color="black",
        linewidth=0.8,
        label="eclipse",
    )
    right.axhline(0, color="gray", linewidth=0.5)
    right.set_xlabel(r"radius ($R_\odot$)")
    right.set_ylabel(f"leak ({rate:latex_inline})")
    right.legend(fontsize=7)

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
only at a signal of $p = \noiselessCenter$, which carries no photon noise.
Middle: $p$ in each region of $64 \times 64$ pixels, with the slope of the line held at its value
near disk center.
Right: its radial profile per second, against the leak measured in the eclipse of the same day.""")
    )
    return result_figure
