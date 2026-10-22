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
\item \textbf{The other channels.}
    The 304\,\AA\ channel shares the telescope and the \CCD\ but has no leak, so the leak is in the
    filters of the 94\,\AA\ half of the aperture.
    What the noise says, and cannot say, about the 131, 171, 193, 211, and 335\,\AA\ channels.
\item \textbf{\SUVI.}
    A joint fit of \AIA\ to \SUVI\ on \GOESShort-16, \GOESShort-18, and \GOESShort-19 with a
    limb-darkened term on the disk and a constant gives the same history, the constant being the
    halo.
\end{outline}""")
    result.append(aia_94_leak_paper.figures.filters())
    return result
