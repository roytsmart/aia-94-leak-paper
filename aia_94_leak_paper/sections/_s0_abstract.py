import aastex

__all__ = [
    "abstract",
]


def abstract() -> aastex.Abstract:
    result = aastex.Abstract()
    result.append(r"""
The 94\,\AA\ channel of \AIA\ on \SDO\ is the instrument's main view of coronal plasma hotter than
a few megakelvin, and with the other \EUV\ channels it underlies most \DEM\ analyses of the corona.
We find that the channel also records visible light from the photosphere.
The leak is a sharp, limb-darkened image of the solar disk with a broad halo, and since each
visible photon frees a single electron, it adds almost none of the photon noise that an \EUV\
signal of the same size would.
It is absent from the 304\,\AA\ channel, which shares the telescope and the detector,
it is attenuated more than the \EUV\ by the thicker of the two focal-plane filters,
and it can be seen directly as \SDO\ passes into the shadow of the Earth, where it sets
\sunsetDelay\ after the \EUV\ has been absorbed, flattened by refraction like a sunset.
In \numEclipses\ such eclipses between \yearFirst\ and \yearLast, the leak at disk center grew
from \leakFirst\ to \leakLast, exponentially since \yearOnset\ with a doubling time of
\doublingTime, and it is now \leakFractionQuietSun\ of the quiet-Sun signal in the channel.
Cross-calibration with \SUVI\ on \GOESShort-18 and \GOESShort-19 and the noise of the channel
itself give the same history.
Left in the data, the leak appears in \DEMs\ as hot plasma in the quiet Sun and in coronal holes,
it biases maps of \FeXVIII, and it hides part of the degradation of the channel from its
calibration.
We give a model of the leak and its history from which it can be subtracted.
\acresetall""")
    return result
