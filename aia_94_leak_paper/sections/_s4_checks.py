import aastex

__all__ = [
    "checks",
]


def checks() -> aastex.Section:
    result = aastex.Section("Independent Checks")
    result.append(r"""
\begin{outline}
\item \textbf{The focal-plane filters.}
    On the one day since 2024 that the channel observed through its thick focal-plane filter, the
    leak fell by more than the \EUV\ did, so it passes through that filter and is light rather
    than an artifact of the detector or its calibration.
\item \textbf{The other channels.}
    The 304\,\AA\ channel shares the telescope and the \CCD\ but has no leak, so the leak is in the
    filters of the 94\,\AA\ half of the aperture.
    What the noise says, and cannot say, about the 131, 171, 193, 211, and 335\,\AA\ channels.
\item \textbf{\SUVI.}
    A joint fit of \AIA\ to \SUVI\ on \GOESShort-16, \GOESShort-18, and \GOESShort-19 with a
    limb-darkened term on the disk and a constant gives the same history, the constant being the
    halo.
\end{outline}""")
    return result
