import aastex

import aia_94_leak_paper

__all__ = [
    "consequences",
]


def consequences() -> aastex.Section:
    result = aastex.Section("Consequences")
    result.append(r"""
\begin{outline}
\item \textbf{Differential emission measures.}
    The \DEM\ of each bin of the six coronal channels at noon of each eclipse day, inverted by the
    method of \citet{Plowman2020} with the response of SolarSoft, corrected for degradation and for
    the lines missing from 94\,\AA~\citep{Boerner2014}, as observed and with the leak measured in the
    eclipse of the same day subtracted from 94\,\AA\ (Figure~\ref{fig:dems}).
    On \demDate, \hotHolesObserved\ of the emission measure of the coronal holes and
    \hotQuietObserved\ of the quiet Sun is above \hotTemperature\ as observed, and
    \hotHolesCorrected\ and \hotQuietCorrected\ with the leak subtracted.
    As observed, both grow with the leak, from at most \hotHolesEarly\ and \hotQuietEarly\ before
    \yearOnset; with the leak subtracted, neither is more than \hotHolesCorrectedMax\ and
    \hotQuietCorrectedMax\ on any day.
    The quiet Sun is chosen by its brightness in 193\,\AA, and near solar maximum it is brighter in
    94\,\AA\ too, so that the leak is a smaller part of it.
\item \textbf{\FeXVIII.}
    In quiet regions the leak survives the subtraction of the cool part of the channel.
\item \textbf{Calibration.}
    In full-disk comparisons a growing leak reads as rising sensitivity.
    The correction tables now hold the sensitivity of the channel constant.
    With the leak subtracted, 94\,\AA\ in the quiet Sun is \ratioFirst\ of what its \DEM\ predicts in
    \yearFirst\ and \ratioLast\ in \yearLast, falling steadily between: the sensitivity of the
    channel to the \EUV\ may have fallen while the leak hid it from full-disk comparisons.
    \roy{Check against the degradation of the other channels and against SUVI.}
\item \textbf{Published results.}
    Which of the results of Section~\ref{sec:Introduction} the leak accounts for, and how much.
\end{outline}""")
    result.append(aia_94_leak_paper.figures.dems())
    return result
