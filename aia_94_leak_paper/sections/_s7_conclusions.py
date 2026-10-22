import aastex

__all__ = [
    "conclusions",
]


def conclusions() -> aastex.Section:
    result = aastex.Section("Conclusions")
    result.append(r"""
\begin{outline}
\item \textbf{The leak.}
    The 94\,\AA\ channel records visible light from the photosphere: a limb-darkened image of the
    disk with a broad halo, which adds almost no photon noise.
    At disk center it grew from \leakFirst\ in \yearFirst\ to \leakLast\ in \yearLast, doubling
    every \doublingTime\ since \yearOnset, and it is now \leakFractionQuietSun\ of the quiet-Sun
    signal in the channel.
\item \textbf{Two measurements.}
    The noise of the channel and the eclipses find it independently.
    Since \yearNoise\ the noiseless signal is \noiseRatioMin\ to \noiseRatioMax\ of the leak in the
    eclipses.
\item \textbf{Its path.}
    The zirconium entrance filters have always passed some visible light, and they pass more and
    more: the images taken without a focal-plane filter brightened by a factor of
    \entranceGrowth\ since \entranceFirst, and since \entranceSince\ the leak has been a constant
    \entranceRatio\ of them.
    The focal-plane filter passes a fixed part of that light through openings in or around it,
    since the thick one attenuates the leak only as much as the \EUV, and 304\,\AA, behind an
    aluminum focal-plane filter, has none.
\item \textbf{Consequences.}
    Left in the data, the leak appears in \DEMs\ as hot plasma outside active regions: on
    \demDate, \hotHolesObserved\ of the emission measure of coronal holes and \hotQuietObserved\ of
    the quiet Sun is above \hotTemperature, and \hotHolesCorrected\ and \hotQuietCorrected\ with the
    leak subtracted.
    Results on hot plasma in the quiet Sun and coronal holes since \yearOnset\ should be revisited.
    With the leak subtracted, 94\,\AA\ in the quiet Sun falls from \ratioFirst\ to \ratioLast\ of
    what its \DEM\ predicts, so the channel may have degraded more than its calibration says.
\item \textbf{What to do.}
    Subtract the model of Section~\ref{sec:ACorrection}, kept current by the eclipse seasons and,
    between them, by the images taken without a focal-plane filter.
    Observing through the thick focal-plane filter would not help, since it attenuates the leak only
    as much as the \EUV.
\item \textbf{For the \AIA\ team.}
    Why the entrance filters pass more and more visible light, where the openings in the
    focal-plane filter are, and where the arcs and streaks in the images taken without a
    focal-plane filter come from.
\end{outline}""")
    return result
