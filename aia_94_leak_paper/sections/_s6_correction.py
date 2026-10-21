import aastex

__all__ = [
    "correction",
]


def correction() -> aastex.Section:
    result = aastex.Section("A Correction")
    result.append(r"""
\begin{outline}
\item \textbf{Model.}
    A limb-darkened disk seen through the \PSF\ of the leak, scaled by an amplitude interpolated
    between eclipses, in the units of the detector and divided by the flat field.
\item \textbf{Uncertainty.}
    From the scatter about the exponential and between the methods.
\item \textbf{Keeping it current.}
    Two eclipse seasons a year measure the leak as it grows.
    \roy{Decide where the correction lives: \texttt{aiapy}, \texttt{sdo}, or this package.}
\end{outline}""")
    return result
