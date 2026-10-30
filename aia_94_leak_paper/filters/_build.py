"""
Build the archive of AIA 94 Å through each of its focal-plane filters from
the level-1 images at JSOC.

The channel has two focal-plane filters, the thin one it normally observes
through and a thick spare. In the routine checks of the spare, every frame
through it is followed by a frame of 304 Å through that channel's own thick
filter, and every such frame carries a band of extra signal across the
middle rows of the detector, about as large as the leak. On 2026 March 7,
the channel observed through the spare all day while it recovered from an
anomaly of the instrument. The last frame through the spare, just before
:data:`switch`, was followed by a frame of 304 Å through its thin filter
instead, and has no band.

The archive holds that frame and the one before it, which has the band, each
averaged over the regions of :data:`aia_94_leak_paper.noise.region` pixels.
It also holds the first run of consecutive frames through the thin filter
after them, accumulated as :func:`aia_94_leak_paper.noise.accumulate` does a
run of the noise archive.

The comparison rests on a single frame, taken during an anomaly. If the AIA
team can take frames through the spare outside an anomaly that are not
followed by a move between the two thick filters, the archive should be
rebuilt from them.
"""

import datetime
import importlib.metadata
import pathlib

import numpy as np
from astropy.io import fits

from .. import _jsoc, noise

__all__ = [
    "average",
    "build",
    "switch",
]

switch: datetime.datetime = datetime.datetime(2026, 3, 7, 23, 20)
"""When the channel changed from its thick focal-plane filter back to its thin one."""

_window = datetime.timedelta(minutes=5)
"""How far from :data:`switch` the frames are searched for."""

_calibration = 1 << 21
"""
The bit of ``QUALITY`` set on the frames of the minutes before
:data:`switch`, through either filter, as on every frame of a routine check
of the thick filter, which the frames through the thick filter may have.
"""


def _consecutive(records: list[dict]) -> bool:
    """Whether a run of frames was taken at the cadence of the channel, 12 s, with none missing."""
    times = [np.datetime64(r["T_OBS"].rstrip("Z")) for r in records]
    span = (times[-1] - times[0]) / np.timedelta64(1, "s")
    return abs(span - 12 * (len(records) - 1)) < 1


def _records() -> dict[str, list[dict]]:
    """
    The last two frames through the thick filter before :data:`switch`, each
    with the ``WAVE_STR`` of the 304 Å frame that followed it, and the first
    usable run of :data:`aia_94_leak_paper.noise.frames` consecutive frames
    through the thin filter after them.
    """
    start = np.datetime64(switch - _window)
    records = _jsoc.records(start, 2 * _window.seconds, 94, ("image", "spikes"))
    telescope = _jsoc.records(start, 2 * _window.seconds, 304)
    before = [r for r in records if r["T_OBS"] < f"{switch:%Y-%m-%dT%H:%M:%S}"]
    thick = [r for r in before if r["WAVE_STR"] == "94_THICK"][-2:]
    if not (
        _consecutive(thick) and noise.usable(thick, "94_THICK", allowed=_calibration)
    ):
        raise ValueError(f"no usable frames through the thick filter before {switch}")
    result = {}
    for name, record, expected in (
        ("banded", thick[0], "304_THICK"),
        ("clean", thick[1], "304_THIN"),
    ):
        following = next(r for r in telescope if r["T_OBS"] > record["T_OBS"])
        if following["WAVE_STR"] != expected:
            raise ValueError(
                f"the frame at {record['T_OBS']} is followed by {following['WAVE_STR']}"
            )
        result[name] = [record | {"FOLLOWED": following["WAVE_STR"]}]
    after = [r for r in records if r["T_OBS"] > thick[-1]["T_OBS"]]
    for i in range(len(after) - noise.frames + 1):
        run = after[i : i + noise.frames]
        if _consecutive(run) and noise.usable(run, "94_THIN"):
            result["thin"] = run
            return result
    raise ValueError(f"no usable run through the thin filter after {switch}")


def average(
    record: dict,
    directory: pathlib.Path,
    flats: pathlib.Path,
) -> fits.HDUList:
    """
    One 94 Å frame in the units of the detector, averaged over the regions of
    :data:`aia_94_leak_paper.noise.region` pixels, leaving out the pixels the
    despiker changed, as the file of the archive holds it.

    Parameters
    ----------
    record
        The frame, as :func:`aia_94_leak_paper._jsoc.records` gives it with
        the ``image`` and ``spikes`` segments, and with the ``WAVE_STR`` of
        the 304 Å frame that followed it as ``FOLLOWED``.
    directory
        Where the level-1 image and its despike record are kept.
    flats
        Where the flat fields are kept.
    """
    stamp = record["T_OBS"][:19].replace(":", "").replace("-", "")
    image = _jsoc.fetch(record["image"], directory / f"aia94_{stamp}.fits")
    spikes = _jsoc.fetch(record["spikes"], directory / f"spikes94_{stamp}.fits")
    with fits.open(image) as hdul:
        header = hdul[1].header
        data = hdul[1].data.astype(np.float64) * _jsoc.flat(header["FLAT_REC"], flats)
    with fits.open(spikes) as hdul:
        spiked = hdul[1].data[0].astype(np.int64)
    data.flat[spiked] = np.nan
    region = noise.region
    n = data.shape[0] // region
    regions = np.nanmean(data.reshape(n, region, n, region), axis=(1, 3))

    primary = fits.PrimaryHDU()
    p = primary.header
    p["T_OBS"] = (record["T_OBS"], "of the frame")
    p["FSN"] = (int(record["FSN"]), "the frame serial number")
    p["EXPTIME"] = (float(record["EXPTIME"]), "exposure, s")
    p["WAVE_STR"] = (record["WAVE_STR"], "the focal-plane filter")
    p["QUALITY"] = (int(record["QUALITY"]), "of the frame")
    p["FOLLOWED"] = (record["FOLLOWED"], "WAVE_STR of the 304 A frame after it")
    p["REGION"] = (region, "pixels along each side of a region")
    p["NSPIKED"] = (len(spiked), "pixels the despiker changed")
    for key in ("CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "R_SUN", "RSUN_OBS"):
        p[key] = (header[key], "of the frame")
    p["FLAT_REC"] = (header["FLAT_REC"], "the flat field the frame is times")
    p["UNITS"] = ("DN", "level 1 times its flat field, per frame")
    p["CREATOR"] = (
        f"aia-94-leak-paper {importlib.metadata.version('aia-94-leak-paper')}",
        "",
    )
    image = fits.ImageHDU(regions.astype(np.float32), name="IMAGE")
    return fits.HDUList([primary, image])


def build(
    cache: pathlib.Path,
    flats: pathlib.Path,
    output: pathlib.Path,
) -> list[pathlib.Path]:
    """
    Write the archive files and return their paths.

    Parameters
    ----------
    cache
        Where the level-1 images and their despike records are kept.
    flats
        Where the flat fields are kept.
    output
        Where to write the archive files.
    """
    output.mkdir(parents=True, exist_ok=True)
    records = _records()
    result = []
    thin = noise.accumulate(records.pop("thin"), cache / "thin", flats)
    thin[0].header["WAVE_STR"] = ("94_THIN", "the focal-plane filter of every frame")
    path = output / "aia94_noise_thin.fits"
    thin.writeto(path, overwrite=True)
    result.append(path)
    for name, (record,) in records.items():
        path = output / f"aia94_frame_{name}.fits"
        average(record, cache / "thick", flats).writeto(path, overwrite=True)
        result.append(path)
    return result
