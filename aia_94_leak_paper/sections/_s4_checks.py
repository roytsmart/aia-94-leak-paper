import aastex

import aia_94_leak_paper

__all__ = [
    "checks",
]


def checks() -> aastex.Section:
    result = aastex.Section("Independent Checks")
    result.append(r"""
\begin{outline}
\item \textbf{The focal-plane filters.}
    The channel has two focal-plane filters: the thin one it observes through, and a spare with
    \filterZirconium\ more zirconium, in case the thin one develops
    pinholes~\citep{Lemen2012,Boerner2012}.
    In the routine checks of the spare, each frame through it is followed by a frame of
    304\,\AA\ through that channel's thick filter, and each carries a band of extra signal across
    the middle rows of the detector (Figure~\ref{fig:filters}, left).
    The band is added as the frame is read out, while the filter wheel moves between the two thick
    filters: a frame through the spare followed by one through the thin filter of 304\,\AA\ has
    none (Figure~\ref{fig:filters}, middle).
    On \filterDate, the day after a reset of the instrument, the channel observed through the
    spare all day, and its last frame through it was followed by one through the thin filter of
    304\,\AA.
    In the frame before, the band reaches \filterBand, about as much as the leak through the thin
    filter at disk center, \filterNoiseless.
    The last frame through the spare and the first \filterFrames\ frames through the thin filter,
    which began \filterGap\ later, are compared within \filterRadius\ of disk center.
    There the spare passes \filterEUV\ of the \EUV, evenly across the field, measured from the
    structure of the corona, where the calibration of the channel expects \filterTable.
    Of the leak, measured from the limb-darkened disk that the frame through the spare lacks, it
    passes $\filterLeak \pm \filterLeakSpread$, the spread among the quadrants of the detector
    (Figure~\ref{fig:filters}, right): none of it.
    \roy{The leak through the spare comes out negative; the noiseless signal through the thin
    filter may read low, as it does against the eclipses (Section~\ref{sec:EarthEclipses}).}
    Light that crossed the metal of the filter would be attenuated by the extra zirconium by
    $e^{-4\pi k\,\Delta d/\lambda}$, where $k$ is the extinction coefficient of zirconium at the
    wavelength $\lambda$ and $\Delta d$ is the extra thickness.
    For the \EUV\ this is \filterMetalEUV\ with the optical constants of \citet{Henke1993}, as the
    calibration expects.
    From \filterVisibleShort\ to \filterVisibleLong, beyond which the silicon of the detector is
    transparent, it is at most \filterMetalQuerry\ with the optical constants of
    \citet{Querry1987}, \filterMetalLynch\ with those of \citet{Lynch1975}, and \filterMetalPalm\
    with those of the film of \citet{Palm2018}, the least absorbing of them.
    So the leak is visible light that crosses the metal of the thin filter, and that the extra
    zirconium of the spare all but stops.
    \roy{Redo this comparison if the \AIA\ team can take frames through the spare outside an
    anomaly of the instrument, each followed by a frame through the thin filter of 304\,\AA: it
    rests on a single frame, taken the day after a reset.}
\item \textbf{The entrance filters.}
    The zirconium entrance filters of the channel pass enough visible light to image the
    photosphere, and \AIA\ monitors them with images taken without a focal-plane
    filter~\citep{Boerner2012}.
    Without one, the light of both halves of the aperture of the 94 and 304\,\AA\ telescope reaches
    the detector: such an image records the visible light that reaches the focal-plane filter.
    The images show the photosphere, darkened toward the limb and with its sunspots, and also
    arcs and streaks of light that are not on the Sun and change over the years
    (Figure~\ref{fig:entrance}, top).
    One a month since the first, in \entranceFirst, they brightened by a factor of \entranceGrowth\
    in the ring at \entranceRing, measured by the median of its bins, which stays in the linear
    range of the camera (Figure~\ref{fig:entrance}, bottom left).
    In every eclipse since \entranceSince, once the leak in that ring was \entranceDominant\ times
    its mean before \yearOnset, it has been \entranceRatio\ of the brightness of these images there,
    with a scatter of \entranceScatter, while it grew by a factor of \entranceGrowthSince\
    (Figure~\ref{fig:entrance}, bottom right).
    Before, the ratio was higher, \entranceRatioFirst\ times as high in \entranceRatioFirstYear,
    when the leak was mostly the faint glow seen before the onset.
    So the leak grows because the entrance filters pass more and more visible light, of which the
    focal-plane filter passes a constant part.
    \roy{The images look more limb-darkened than the model of the leak; compare the radial
    profiles.}
    \roy{The arcs and streaks cross the ring; show that the ratio does not depend on the
    statistic of the ring, and find where they come from.}
\item \textbf{The other channels.}
    The 304\,\AA\ channel shares the telescope, the \CCD, and the visible light that reaches the
    focal-plane filter, but has no leak, so the leak passes the zirconium focal-plane filter of
    94\,\AA\ and not the aluminum one of 304\,\AA.
    What the noise says, and cannot say, about the 131, 171, 193, 211, and 335\,\AA\ channels.
\end{outline}""")
    result.append(aia_94_leak_paper.figures.filters())
    result.append(aia_94_leak_paper.figures.entrance())
    return result
