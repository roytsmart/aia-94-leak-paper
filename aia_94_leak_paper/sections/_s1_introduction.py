import aastex

__all__ = [
    "introduction",
]


def introduction() -> aastex.Section:
    result = aastex.Section("Introduction")
    result.append(r"""
\begin{outline}
\item \textbf{The channel.}
    \AIACapital~\citep{Lemen2012} images the corona in seven \EUV\ channels from \SDO~\citep{Pesnell2012}.
    The 94\,\AA\ channel sees \FeXVIII\ from plasma near 7\,MK, and cooler lines of \FeIX\ and \FeX\
    which atomic databases underestimate~\citep{Testa2012}, so its temperature response was
    corrected empirically against \EVE~\citep{Woods2012, Boerner2014}.
\item \textbf{What depends on it.}
    It is the main constraint on hot plasma in \DEM\ inversions~\citep{Hannah2012, Plowman2013,
    Cheung2015, Morgan2019, Plowman2020}, and maps of \FeXVIII\ are made by subtracting the cool
    part of it predicted from the other channels~\citep{Warren2012, DelZanna2013}.
\item \textbf{Calibration over the mission.}
    The degradation of each channel is tracked against \EVE\ and its successors~\citep{Boerner2012,
    Boerner2014, DosSantos2021} and distributed in \texttt{aiapy}~\citep{Barnes2020}.
    A signal that grows with time reads as a channel which degrades more slowly than it does.
\item \textbf{A precedent.}
    The zirconium filters of the 94\,\AA\ channels transmit a little visible
    light~\citep{Tadikonda2019}, and pinholes in the entrance filters of \SUVI\ on
    \GOESShort-17 leak visible light, most in its 94\,\AA\ channel~\citep{Darnel2022}.
    \roy{Check that Darnel2022 describes the GOES-17 leak; otherwise cite the SUVI L1b readme.}
\item \textbf{Symptoms already in the literature.}
    In coronal holes the 94\,\AA\ signal matches what remains behind an occulter, and was taken
    for stray light~\citep{Saqri2020, Heinemann2021}.
    A survey of quiet-Sun \DEMs\ ended in 2018 because the reconstructions after it were
    poor~\citep{Purkhart2022}, and some inversions drop the channel~\citep{Hofmeister2025}.
    Coronal holes lose their contrast in it~\citep{Lezzi2023}, and it is brighter than a
    \DEM\ constrained in X-rays allows~\citep{Athiray2024}.
\item \textbf{Roadmap.}
    Section~\ref{sec:PhotonNoise} finds a noiseless signal in the channel,
    Section~\ref{sec:EarthEclipses} measures it directly in eclipses over the mission,
    Section~\ref{sec:IndependentChecks} confirms it with the focal-plane filters, the other
    channels, and \SUVI, Section~\ref{sec:Consequences} follows it into thermal diagnostics and
    calibration, Section~\ref{sec:ACorrection} gives a correction, and
    Section~\ref{sec:Conclusions} concludes.
\end{outline}""")
    return result
