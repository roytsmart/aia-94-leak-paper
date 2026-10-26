import datetime
import itertools

import aastex
import astropy.units as u
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage

import aia_94_leak_paper

__all__ = [
    "dates_residuals",
    "residuals",
]

dates_residuals: tuple[datetime.date, ...] = (
    datetime.date(2017, 8, 25),
    datetime.date(2020, 2, 13),
    datetime.date(2023, 8, 1),
    datetime.date(2026, 1, 19),
)
"""
The eclipses shown, one every two and a half years since soon after the
onset, alternating between the seasons, in which the edge sweeps the disk in
different directions.
"""

smoothing: int = 5
"""The side of the squares of bins the map of the plateau is averaged over."""

regions: dict[str, tuple[float, float]] = {
    "disk": (0, 0.95),
    "off the limb": (1.05, 1.4),
}
"""The distances from the center of the disk the frames are pooled over, in solar radii."""


def _smooth(image: np.ndarray) -> np.ndarray:
    """An image averaged over squares of :data:`smoothing` bins, the bins left out ignored."""
    good = np.isfinite(image)
    total = scipy.ndimage.uniform_filter(np.where(good, image, 0), smoothing)
    count = scipy.ndimage.uniform_filter(good.astype(float), smoothing)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(good, total / count, np.nan)


def _clipped_mean(values: np.ndarray) -> float:
    """
    The mean of values within five robust standard deviations of their
    median, against the particle hits the despiker missed.
    """
    median = np.median(values)
    spread = 1.4826 * np.median(np.abs(values - median))
    return float(np.mean(values[np.abs(values - median) < 5 * spread]))


def residuals() -> aastex.FigureStar:
    """
    The fit to the ingress of a few eclipses: the mean residual over the
    plateau of each bin, and the frames, the fit, and the residual against
    the time after the EUV edge of each bin.
    """
    eclipses = aia_94_leak_paper.eclipses
    fig, axs = plt.subplots(
        nrows=3,
        ncols=len(dates_residuals),
        figsize=(aastex.text_width_inches, 5.2),
        sharex="row",
        sharey="row",
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.3, 0.8, 0.8]},
    )
    axes = ("time", "y", "x")
    seconds = np.arange(35, 201)
    styles = {
        "disk": ("black", "tab:red"),
        "off the limb": ("0.55", "tab:orange"),
    }
    for col, date in enumerate(dates_residuals):
        center = eclipses.fit(date).center
        result = eclipses.residuals(date)

        # everything in percent of the leak at disk center, so that the
        # eclipses can be compared
        data, model = (
            (a / center).to_value(u.percent).ndarray_aligned(axes)
            for a in (result.data, result.model)
        )
        residual = data - model
        s = result.s.ndarray_aligned(axes).to_value(u.s)
        use = result.use.ndarray_aligned(axes)
        r = np.broadcast_to(result.r.ndarray_aligned(("y", "x")), s.shape)
        setting = result.setting.to_value(u.s)

        # the mean residual over the plateau of each bin
        image = (
            (result.plateau / center).to_value(u.percent).ndarray_aligned(("y", "x"))
        )
        x, y = result.position.x.ndarray, result.position.y.ndarray
        top = axs[0, col]
        top.set_facecolor("0.85")
        im = top.imshow(
            _smooth(image),
            origin="lower",
            extent=[x[0], x[-1], y[0], y[-1]],
            cmap="RdBu_r",
            vmin=-10,
            vmax=10,
            interpolation="none",
        )
        top.add_patch(
            matplotlib.patches.Circle((0, 0), 1, fill=False, color="0.3", linewidth=0.4)
        )
        # the edges of the quadrants of the detector, each read out by its own amplifier
        n = x.size // 2
        top.axvline(
            (x[n - 1] + x[n]) / 2, color="0.2", linewidth=0.4, linestyle="dashed"
        )
        top.axhline(
            (y[n - 1] + y[n]) / 2, color="0.2", linewidth=0.4, linestyle="dashed"
        )
        # the direction the edge sweeps the disk: the clock of each bin
        # starts later along it
        xx, yy = np.meshgrid(x, y)
        disk = r[0] < 0.95
        design = np.stack([np.ones(disk.sum()), xx[disk], yy[disk]], axis=-1)
        (_, sx, sy), *_ = np.linalg.lstsq(design, -s[0][disk], rcond=None)
        dx, dy = np.array([sx, sy]) / np.hypot(sx, sy)
        top.annotate(
            "",
            xy=(0.6 * dx, 0.6 * dy),
            xytext=(-0.6 * dx, -0.6 * dy),
            arrowprops={"arrowstyle": "->", "linewidth": 0.6, "color": "0.2"},
        )
        top.set_xlim(-1.1, 1.1)
        top.set_ylim(-1.1, 1.1)
        top.set_xticks([])
        top.set_yticks([])
        top.set_title(f"{date:%Y} {date:%B} {date.day}", fontsize=8)

        # the frames, the fit, and the residual against the time after the
        # EUV edge of each bin, pooled over the bins of each region
        middle, bottom = axs[1, col], axs[2, col]
        for region, (inner, outer) in regions.items():
            where = use & (r >= inner) & (r < outer)
            t, d, m, e = [], [], [], []
            for lo, hi in itertools.pairwise(seconds):
                k = where & (s >= lo) & (s < hi)
                if k.sum() < 200:
                    continue
                t.append((lo + hi) / 2)
                d.append(_clipped_mean(data[k]))
                m.append(np.mean(model[k]))
                e.append(_clipped_mean(residual[k]))
            color_data, color_model = styles[region]
            middle.plot(t, d, ".", color=color_data, markersize=1.5, label=f"{region}")
            middle.plot(t, m, color=color_model, linewidth=0.8, label="fit")
            bottom.plot(t, e, color=color_model, linewidth=0.8)
        bottom.axhspan(-1, 1, color="0.92", linewidth=0, zorder=0)
        bottom.axhline(0, color="0.6", linewidth=0.5)
        for ax in (middle, bottom):
            ax.axvline(setting, color="0.4", linewidth=0.5, linestyle="dotted")
            ax.set_xlim(seconds[0], seconds[-1])
            ax.tick_params(labelsize=7)
        middle.set_ylim(-5, 105)
        bottom.set_ylim(-8, 8)

    axs[1, 0].set_ylabel(r"\% of the leak" "\n" "at disk center", fontsize=8)
    axs[2, 0].set_ylabel(r"residual (\%)", fontsize=8)
    axs[1, -1].legend(
        fontsize=6, loc="center right", ncols=2, handlelength=1, columnspacing=0.8
    )
    fig.supxlabel("time after the EUV edge of each bin (s)", fontsize=8)
    bar = fig.colorbar(im, ax=axs[0, :], shrink=0.85, pad=0.01)
    bar.set_label(r"residual (\%)", fontsize=8)
    bar.ax.tick_params(labelsize=7)

    plateau = eclipses.before_setting.to_value(u.s)
    (_, disk), (inner, outer) = regions["disk"], regions["off the limb"]
    result = aastex.FigureStar("fig:residuals", position="!htb")
    result.append(aastex.NoEscape(r"\centering"))
    result.add_fig(fig, width=aastex.NoEscape(r"\textwidth"))
    result.add_caption(
        aastex.NoEscape(rf"""
The fit to the ingress of four eclipses, in percent of the leak at disk center.
Top: the mean residual of each bin over its plateau, from when its \EUV\ is gone until
{plateau:g}\,s before the visible disk sets there, averaged over squares of
${smoothing} \times {smoothing}$ bins; gray where no frame is fit, dashed at the edges of the
quadrants of the detector, and the arrow the direction the \EUV\ edge sweeps the disk.
Middle: the frames (points) and the fit (lines) against the time after the \EUV\ edge of each
bin, the mean of the bins on the disk (within {disk:g}\,$R_\odot$) and off the limb ({inner:g} to
{outer:g}\,$R_\odot$) in each second.
Bottom: their difference; the gray band is $\pm 1\%$.
The dotted line is when the visible disk is half set.""")
    )
    return result
