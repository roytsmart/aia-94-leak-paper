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
    On \filterDate\ the channel observed through the spare until late in the day; the last run of
    frames through it and the first through the thin filter began \filterGap\ apart.
    The spare transmits unevenly, passing less toward the top and bottom of the field, and its
    flat field is that of the thin filter, so the two are compared only within \filterBand\ of the
    middle of the disk (Figure~\ref{fig:filters}, left).
    There the spare passes \filterEUV\ of the \EUV, measured from the structure of the corona,
    where the calibration of the channel expects \filterTable.
    It passes $\filterLeak \pm \filterLeakSpread$ of the noiseless signal, and \filterLeakSignals\
    of it from the signals alone (Figure~\ref{fig:filters}, right): it attenuates the leak about as
    much as the \EUV.
    So the leak reaches the detector through the focal-plane filter.
    Visible light passing through the metal of the filter would be attenuated by orders of
    magnitude more by the extra zirconium, so the leak passes the filter through openings in it or
    around it rather than through the metal.
    \roy{Cite the optical constants of zirconium in the visible.}
\item \textbf{The entrance filters.}
    The zirconium entrance filters of the channel pass enough visible light to image the
    photosphere, and \AIA\ monitors them with images taken without a focal-plane
    filter~\citep{Boerner2012}.
    Without one, the light of both halves of the aperture of the 94 and 304\,\AA\ telescope reaches
    the detector: such an image records the visible light that reaches the focal-plane filter.
    One a month since the first, in \entranceFirst, they brightened by a factor of \entranceGrowth\
    in the ring at \entranceRing, where the detector stays in its linear range
    (Figure~\ref{fig:entrance}, left).
    In every eclipse since \entranceSince, once the leak in that ring was \entranceDominant\ times
    its mean before \yearOnset, it has been \entranceRatio\ of the brightness of these images there,
    with a scatter of \entranceScatter, while it grew by a factor of \entranceGrowthSince\
    (Figure~\ref{fig:entrance}, right).
    Before, the ratio was higher, \entranceRatioFirst\ times as high in \entranceRatioFirstYear,
    when the leak was mostly the faint glow seen before the onset.
    So the leak grows because the entrance filters pass more and more visible light, of which the
    focal-plane filter passes a constant part.
    \roy{The images look more limb-darkened than the model of the leak; compare the radial
    profiles.}
\item \textbf{The other channels.}
    The 304\,\AA\ channel shares the telescope, the \CCD, and the visible light that reaches the
    focal-plane filter, but has no leak, so the leak passes the zirconium focal-plane filter of
    94\,\AA\ and not the aluminum one of 304\,\AA.
    What the noise says, and cannot say, about the 131, 171, 193, 211, and 335\,\AA\ channels.
\end{outline}""")
    result.append(aia_94_leak_paper.figures.filters())
    result.append(aia_94_leak_paper.figures.entrance())
    return result
