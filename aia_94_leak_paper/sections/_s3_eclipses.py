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
    The 171\,\AA\ channel marks where the \EUV\ has gone (Figure~\ref{fig:ingress}).
\item \textbf{Method.}
    Each bin of the disk on a clock started by its own 171\,\AA\ edge; a forward model of a
    limb-darkened disk setting behind the atmosphere, seen through a \PSF\ with a sharp core and
    broad wings, plus whatever \EUV\ is left, fit to the frames in which the \EUV\ is gone.
    The leak reaches the same detector as the \EUV, so the 94\,\AA\ frames are flat fielded as the
    \EUV\ is, and multiplied by the flat field at disk center to keep the units of the detector
    there; in the units of the detector the leak carries the pattern of the flat field, with a step
    at the edge of each quadrant from the gain of its amplifier.
    The dark is the median of the umbra just after the whole disk has set, rather than minutes
    later, since the offset of each quadrant drifts.
\item \textbf{Results.}
    The leak at disk center in \numEclipses\ eclipses from \yearFirst\ to \yearLast: from
    \leakFirst\ to \leakLast, with a doubling time of \doublingTime\ since \yearOnset\
    (Figure~\ref{fig:leakHistory}).
    Its limb darkening, the sharpness of its limb, and the halo beyond it.
\item \textbf{Against the noise.}
    The noiseless signal near disk center in the run at noon of each day follows the leak
    (Figure~\ref{fig:leakHistory}).
    Since \yearNoise, when the leak first exceeded three times the spread of the noise method among
    the quadrants (\noiseSpread), the noiseless signal is \noiseRatioMin\ to \noiseRatioMax\ of the
    leak, the fraction falling as the leak grows; before, the leak is too faint for the noise to
    measure.
    Part of the shortfall is expected: a visible photon frees one electron, so the leak does carry
    a little photon noise, a twentieth of what an \EUV\ signal of the same size would.
    \roy{Account for the rest: the rise of the variance per DN with signal, the step in the
    noiseless signal between the quadrants, and the change of the read noise in 2018.}
\item \textbf{Caveats.}
    The \PSF\ is fit once to the brightest eclipses and held fixed, since in a single eclipse its
    core and wing trade off against each other, though the light they deliver to the disk does not.
    Before \yearOnset\ the leak is too faint to time its setting, which is held at its median since,
    and it lies above the exponential the later eclipses follow.
    One ingress fell during a flare.
\item \textbf{Figures.}
    Still to come: the radial profile of the leak, and the sharpness of its limb.
\end{outline}""")
    result.append(aia_94_leak_paper.figures.ingress())
    result.append(aia_94_leak_paper.figures.leak_history())
    return result
