import aastex

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
    \leakFirst\ to \leakLast, with a doubling time of \doublingTime\ since \yearOnset.
    Its limb darkening, the sharpness of its limb, and the halo beyond it.
\item \textbf{Caveats.}
    The core and the wings of the \PSF\ are degenerate, though the light they deliver to the disk
    is not; one ingress fell during a flare; before \yearOnset\ there is only a faint diffuse glow.
\item \textbf{Figures.}
    An ingress in 94\,\AA\ and 171\,\AA; the leak at disk center against time; its radial profile.
\end{outline}""")
    return result
