import aastex

import aia_94_leak_paper

__all__ = [
    "noise",
]


def noise() -> aastex.Section:
    result = aastex.Section("Photon Noise")
    result.append(r"""
\begin{outline}
\item \textbf{Data.}
    A run of \noiseFrames\ consecutive level~1 frames at noon on each day of an eclipse, returned to
    the units of the detector by multiplying each by the flat field it was divided by.
\item \textbf{Method.}
    The residual of each frame from the mean of its two neighbors cancels the Sun and anything
    that changes linearly with time, and has 3/2 the variance of one frame.
    It is sorted by the signal of the other frames of the run, whose noise is independent of it.
    The read noise, \readNoise, comes from the corners of the \CCD\ beyond the field stop.
\item \textbf{The noiseless signal.}
    Off the limb the variance rises from the read noise with the signal, as photon noise should.
    Near disk center the line of slope \noiseSlope\ meets the read noise only at a signal of
    \noiselessCenter: that much of the signal carries no photon noise (Figure~\ref{fig:noise}).
    Its map is a limb-darkened disk which vanishes off the limb and does not follow the structure
    of the corona.
\item \textbf{Its history.}
    The noiseless signal follows the leak measured in the eclipses of the same days
    (Figure~\ref{fig:leakHistory}), at \noiseRatio\ of it since \yearOnset.
    Part of the shortfall is expected: a visible photon frees one electron, so the leak does carry
    a little photon noise, a twentieth of what an \EUV\ signal of the same size would.
    \roy{Account for the rest: the rise of the variance per DN with signal, and the wing.}
\end{outline}""")
    result.append(aia_94_leak_paper.figures.noise())
    return result
