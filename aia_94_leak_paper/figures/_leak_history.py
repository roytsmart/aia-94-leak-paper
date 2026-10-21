import aastex
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import aia_94_leak_paper

__all__ = [
    "leak_history",
]


def leak_history() -> aastex.Figure:
    """
    The leak at disk center and just off the limb in every eclipse, against
    time, with the exponential fit to the eclipses since the onset.
    """
    eclipses = aia_94_leak_paper.eclipses
    leaks = eclipses.history()
    growth = eclipses.growth(leaks)

    year = np.array([leak.time.decimalyear for leak in leaks])
    rate = u.DN / u.s
    center = u.Quantity([leak.center for leak in leaks]).to_value(rate)
    glow = u.Quantity([leak.glow for leak in leaks]).to_value(rate)

    fig, ax = plt.subplots(
        figsize=(aastex.column_width_inches, 2.4),
        constrained_layout=True,
    )
    ax.plot(year, center, "o", color="tab:red", markersize=4, label="disk center")
    ax.plot(
        year,
        glow,
        "o",
        color="tab:red",
        markerfacecolor="none",
        markersize=4,
        label=r"1.05--1.15\,$R_\odot$",
    )
    noise = aia_94_leak_paper.noise.history()
    ax.plot(
        [n.time.decimalyear for n in noise],
        u.Quantity([n.rate for n in noise]).to_value(rate),
        "^",
        color="tab:blue",
        markersize=4,
        label="noiseless, disk center",
    )
    t = np.linspace(growth.since, year.max(), 100)
    ax.plot(
        t, growth(t).to_value(rate), color="black", linewidth=0.8, label="exponential"
    )

    ax.set_yscale("log")
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    ax.set_xlabel("year")
    ax.set_ylabel(f"leak ({rate:latex_inline})")
    ax.legend(fontsize=7, loc="upper left")

    result = aastex.Figure("fig:leakHistory", position="!htb")
    result.append(aastex.NoEscape(r"\centering"))
    result.add_fig(fig, width=aastex.NoEscape(r"\columnwidth"))
    result.add_caption(
        aastex.NoEscape(r"""
The leak in the \numEclipses\ eclipse ingresses, in the units of the detector, at disk center
(filled circles) and just off the limb (open circles), each with the whole disk up, and the
noiseless signal near disk center in a run at noon of each day (triangles).
The line is an exponential fit to the eclipses since \yearOnset, which doubles every
\doublingTime.""")
    )
    return result
