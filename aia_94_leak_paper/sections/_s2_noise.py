import aastex

__all__ = [
    "noise",
]


def noise() -> aastex.Section:
    result = aastex.Section("Photon Noise")
    result.append(r"""
\begin{outline}
\item \textbf{Data.}
    Runs of consecutive level~1 images at the full cadence of the channel, returned to the units
    of the detector by multiplying by the flat field they were divided by.
\item \textbf{Method.}
    The variance of the second difference of three consecutive frames, which cancels the linear
    evolution of the Sun, against the signal estimated from frames independent of it.
    The read noise comes from the corners of the \CCD\ beyond the field stop.
\item \textbf{The noiseless signal.}
    Off the limb the variance is the read noise plus a constant times the signal, as expected.
    On the disk it is offset: part of the signal carries no photon noise.
    Its map is a limb-darkened disk which vanishes off the limb and does not follow the structure
    of the corona in any channel.
\item \textbf{Its history.}
    The noiseless signal at disk center, from runs spread over the mission.
\item \textbf{Figures.}
    Variance against signal on and off the disk; the map of the noiseless signal; its history.
\end{outline}""")
    return result
