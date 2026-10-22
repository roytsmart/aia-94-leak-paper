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
\item \textbf{Zirconium and visible light.}
    The 94\,\AA\ channel shares a telescope and a \CCD\ with 304\,\AA.
    Each has half of the aperture, behind its own entrance filter, and a filter wheel in front of
    the detector chooses between them with a focal-plane filter: zirconium for 94\,\AA\ and aluminum
    for 304\,\AA~\citep{Lemen2012}.
    Zirconium has always passed some visible light.
    The zirconium entrance filters of the 94 and 131\,\AA\ channels pass enough of it to image the
    photosphere, and \AIA\ monitors them routinely with images taken without a focal-plane
    filter~\citep{Boerner2012}.
    It is the focal-plane filter that keeps this light from the detector, and the filter wheel
    holds a thicker spare in case the thin filter develops pinholes~\citep{Lemen2012}.
\item \textbf{Symptoms already in the literature.}
    In coronal holes the 94\,\AA\ signal matches what remains behind an occulter, and was taken
    for stray light~\citep{Saqri2020, Heinemann2021}.
    A survey of quiet-Sun \DEMs\ ended in 2018 because the reconstructions after it were
    poor~\citep{Purkhart2022}, and some inversions drop the channel~\citep{Hofmeister2025}.
    Coronal holes lose their contrast in it~\citep{Lezzi2023}, and it is brighter than a
    \DEM\ constrained in X-rays allows~\citep{Athiray2024}.
\item \textbf{Roadmap.}
    Section~\ref{sec:PhotonNoise} finds a signal in the channel that carries no photon noise,
    Section~\ref{sec:EarthEclipses} sees it alone in eclipses and measures its history over the
    mission, Section~\ref{sec:IndependentChecks} follows its path through the instrument with the
    two focal-plane filters, the images taken without one, and the other channels,
    Section~\ref{sec:Consequences} follows it into \DEMs, maps of \FeXVIII, and the calibration of
    the channel, Section~\ref{sec:ACorrection} gives a correction, and
    Section~\ref{sec:Conclusions} concludes.
\end{outline}""")
    return result
