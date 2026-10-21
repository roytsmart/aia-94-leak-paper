import aastex

__all__ = [
    "consequences",
]


def consequences() -> aastex.Section:
    result = aastex.Section("Consequences")
    result.append(r"""
\begin{outline}
\item \textbf{Differential emission measures.}
    The fraction of the \DEM\ of coronal holes and the quiet Sun above 3\,MK, early in the mission
    and now, with and without the leak subtracted.
    Subtracting it restores the quiet Sun, but coronal holes stay hotter than they were early in
    the mission.
\item \textbf{\FeXVIII.}
    In quiet regions the leak survives the subtraction of the cool part of the channel.
\item \textbf{Calibration.}
    In full-disk comparisons a growing leak reads as rising sensitivity.
    The correction tables now hold the sensitivity of the channel constant.
\item \textbf{Published results.}
    Which of the results of Section~\ref{sec:Introduction} the leak accounts for, and how much.
\end{outline}""")
    return result
