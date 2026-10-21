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
    The read noise, about \readNoise, comes from the corners of the \CCD\ beyond the field stop;
    each quadrant of the \CCD\ is read out by its own amplifier, and is fit on its own.
\item \textbf{The noiseless signal.}
    Off the limb the variance rises from the read noise with the signal, as photon noise alone
    would, at the slope found near disk center.
    Near disk center the line of slope \noiseSlope\ meets the read noise only at a signal of
    $\noiselessCenter \pm \noiselessSpread$: that much of the signal carries no photon noise
    (Figure~\ref{fig:noise}, left).
    The photon noise of a signal grows with the electrons each photon frees in the silicon:
    dozens for a 94\,\AA\ photon, one for a visible photon, so a signal of visible light carries
    almost none.
\item \textbf{Its map.}
    The noiseless signal fills a disk the size of the photosphere (Figure~\ref{fig:noise}, middle),
    and does not follow the corona: across the disk it grows by \noisePerSignal\ per DN of the
    94\,\AA\ signal of a region, though it is \noiseShare\ of that signal.
    It darkens toward the limb as $1 - u (1 - \mu)$ with $u = \noiseLimbDarkening$
    (Figure~\ref{fig:noise}, right), where $\mu = \sqrt{1 - r^2}$ is the cosine of the angle between
    the line of sight and the vertical.
    This is limb darkening, the signature of an optically thick photosphere: we see into it only to
    a fixed optical depth along the line of sight, which near the limb, at a grazing angle, lies in
    higher and cooler layers, so the disk dims toward its edge.
    The optically thin corona does the opposite, brightening toward the limb, where the line of sight
    through it is longest.
    The channel sees the photosphere, in light that frees few electrons per photon: a visible light
    leak.
    \roy{Cite a measured limb darkening of the photosphere in visible light (\eg, Neckel and Labs
    1994) to compare $u$ with.}
\end{outline}""")
    result.append(aia_94_leak_paper.figures.noise())
    return result
