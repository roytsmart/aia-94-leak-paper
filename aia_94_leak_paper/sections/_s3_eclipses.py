import aastex

import aia_94_leak_paper

__all__ = [
    "eclipses",
]


def eclipses() -> aastex.Section:
    result = aastex.Section("Earth Eclipses")
    result.append(r"""
\begin{outline}
\item \textbf{Eclipse seasons.}
    \SDOCapital\ passes behind the Earth once a day for a few weeks twice a year, and the seasons drift
    earlier from year to year.
\item \textbf{Ingress.}
    The atmosphere absorbs the \EUV\ high above the surface, while the visible disk sets only near
    it, so for a while each part of the disk is seen in visible light alone.
    The 171\,\AA\ channel marks where the \EUV\ has gone.
\item \textbf{Method.}
    Each bin of the disk on a clock started by its own 171\,\AA\ edge; a forward model of a
    limb-darkened disk setting behind the atmosphere, seen through a \PSF\ with a sharp core and
    broad wings, plus whatever \EUV\ is left, fit to the frames in which the \EUV\ is gone.
    The frames in the umbra are the dark.
\item \textbf{Results.}
    The leak at disk center in \numEclipses\ eclipses from \yearFirst\ to \yearLast: from
    \leakFirst\ to \leakLast, with a doubling time of \doublingTime\ since \yearOnset\
    (Figure~\ref{fig:leakHistory}).
    Its limb darkening, the sharpness of its limb, and the halo beyond it.
\item \textbf{Caveats.}
    The \PSF\ is fit once to the brightest eclipses and held fixed, since in a single eclipse its
    core and wing trade off against each other, though the light they deliver to the disk does not.
    Before \yearOnset\ the leak is too faint to time its setting, which is held at its median since,
    and it lies above the exponential the later eclipses follow.
    One ingress fell during a flare.
\item \textbf{Figures.}
    An ingress in 94\,\AA\ and 171\,\AA; the leak at disk center against time; its radial profile.
\end{outline}""")
    result.append(aia_94_leak_paper.figures.leak_history())
    return result
