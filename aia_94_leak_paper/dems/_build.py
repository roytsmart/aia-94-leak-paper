"""
Build the archive of the six coronal channels of AIA from the synoptic
images at JSOC.

For each day of :data:`aia_94_leak_paper.eclipses.dates`, the synoptic
images of the six channels a differential emission measure is inverted from
are fetched from the first slot at or after noon in which all six were
taken and none is flagged, hours away from the eclipse and at the time of
the noise run. The
synoptic images are level-1.5 images, registered to one another and to
solar north, and averaged down 4 x 4 to 1024 x 1024 pixels of 2.4 arcsec;
they are averaged down again 2 x 2 to bins of 4.8 arcsec. With them is kept
the flat field the 94 Å image was divided by, carried onto the same bins by
the master pointing in its header, since a leak that adds to the detector is
divided by the flat field with everything else.
"""

import datetime
import importlib.metadata
import pathlib
import re
import urllib.request

import numpy as np
import scipy.ndimage
from astropy.io import fits

from .. import _jsoc

__all__ = [
    "binning",
    "build",
    "slots",
    "url_synoptic",
    "wavelengths",
]

wavelengths: tuple[int, ...] = (94, 131, 171, 193, 211, 335)
"""The channels a differential emission measure is inverted from, in angstroms."""

binning: int = 2
"""The synoptic pixels along each side of a bin."""

url_synoptic: str = "https://jsoc1.stanford.edu/data/aia/synoptic"
"""Where JSOC keeps the definitive synoptic images, by day and hour."""

_noon = 12
"""The hour of each day from which the first slot with all six channels is taken."""

_samples = 16
"""The level-1.5 pixels each synoptic pixel is the mean of."""

_harmless = 0b100
"""
The bits of ``QUALITY`` an image may have set: only that its ancillary
science data are not available, which changes no pixel, and which is set on
every image of 2026 January 19.
"""


def _listing(url: str) -> str:
    """The text of a directory listing at JSOC."""
    return _jsoc.retry(
        lambda: urllib.request.urlopen(url, timeout=120).read().decode(),
        what=url,
    )


def slots(date: datetime.date) -> list[datetime.datetime]:
    """
    The synoptic slots in the two hours from noon of a day in which all six
    channels were taken, in order.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    result = []
    for hour in (_noon, _noon + 1):
        listing = _listing(f"{url_synoptic}/{date:%Y/%m/%d}/H{hour:02d}00/")
        taken: dict[str, set[int]] = {}
        for hhmm, channel in re.findall(
            rf"AIA{date:%Y%m%d}_(\d{{4}})_(\d{{4}})\.fits", listing
        ):
            taken.setdefault(hhmm, set()).add(int(channel))
        for hhmm in sorted(taken):
            if set(wavelengths) <= taken[hhmm]:
                time = datetime.time(int(hhmm[:2]), int(hhmm[2:]))
                result.append(datetime.datetime.combine(date, time))
    return result


def _fetch(when: datetime.datetime, channel: int, cache: pathlib.Path) -> pathlib.Path:
    """A synoptic image, downloaded unless it is already in the cache."""
    name = f"AIA{when:%Y%m%d_%H%M}_{channel:04d}.fits"
    path = cache / name
    if not path.exists():
        url = f"{url_synoptic}/{when:%Y/%m/%d}/H{when:%H}00/{name}"
        data = _jsoc.retry(
            lambda: urllib.request.urlopen(url, timeout=120).read(), what=url
        )
        cache.mkdir(parents=True, exist_ok=True)
        path.with_suffix(".part").write_bytes(data)
        path.with_suffix(".part").replace(path)
    return path


def _flat(header: fits.Header, flat: np.ndarray, n: int) -> np.ndarray:
    """
    A level-1 flat field carried onto the bins of a synoptic image.

    The flat field is averaged over the detector pixels of a bin and sampled
    at the center of each bin, found by undoing the registration of the
    synoptic image: the master pointing gives the center of the disk on the
    detector (``X0_MP``, ``Y0_MP``), its plate scale (``IMSCL_MP``), and,
    with the roll of the spacecraft, its rotation from solar north.
    """
    scale = header["CDELT1"] * binning
    pixels = scale / header["IMSCL_MP"]
    flat = scipy.ndimage.uniform_filter(flat, size=round(pixels))
    i, j = np.indices((n, n))
    x = (j - (n - 1) / 2) * scale
    y = (i - (n - 1) / 2) * scale
    # the roll is not recorded when the ancillary science data are missing,
    # as on 2026 January 19, and is under a millidegree when it is
    roll = header["SAT_ROT"] or 0
    theta = np.deg2rad(roll + header["INST_ROT"])
    columns = (
        header["X0_MP"] + (np.cos(theta) * x + np.sin(theta) * y) / header["IMSCL_MP"]
    )
    rows = (
        header["Y0_MP"] + (-np.sin(theta) * x + np.cos(theta) * y) / header["IMSCL_MP"]
    )
    return scipy.ndimage.map_coordinates(flat, [rows, columns], order=1, mode="nearest")


def build(
    date: datetime.date,
    cache: pathlib.Path,
    flats: pathlib.Path,
    output: pathlib.Path,
) -> pathlib.Path:
    """
    Write the archive file of one day.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    cache
        Where the synoptic images are kept.
    flats
        Where the flat fields are kept.
    output
        Where to write the archive file.
    """
    for when in slots(date):
        images, headers = [], []
        for channel in wavelengths:
            with fits.open(_fetch(when, channel, cache)) as hdul:
                header = hdul[1].header
                if int(header["WAVELNTH"]) != channel:
                    raise ValueError(f"{when} {channel} A is {header['WAVELNTH']} A")
                images.append(hdul[1].data.astype(np.float64))
                headers.append(header)
        if not any(h["QUALITY"] & ~_harmless for h in headers):
            break
    else:
        raise ValueError(
            f"no slot of {date} after noon is clean in all of {wavelengths}"
        )
    images = np.array(images)
    n = images.shape[-1] // binning
    images = images.reshape(len(wavelengths), n, binning, n, binning).mean(axis=(2, 4))
    first = headers[0]
    flat = _flat(first, _jsoc.flat(first["FLAT_REC"], flats), n)

    primary = fits.PrimaryHDU()
    p = primary.header
    p["DATE_DEM"] = (f"{date:%Y-%m-%d}", "the day")
    p["SLOT"] = (f"{when:%H:%M}", "the synoptic slot")
    for channel, header in zip(wavelengths, headers):
        p[f"T_{channel:04d}"] = (header["T_OBS"], f"T_OBS of {channel} A")
        p[f"X_{channel:04d}"] = (header["EXPTIME"], f"exposure of {channel} A, s")
    p["RSUN_OBS"] = (first["RSUN_OBS"], "radius of the photosphere, arcsec")
    p["CDELT"] = (first["CDELT1"] * binning, "side of a bin, arcsec")
    p["NSAMPLE"] = (_samples * binning**2, "level-1.5 pixels in the mean of a bin")
    p["FLAT_REC"] = (first["FLAT_REC"], "the flat field of 94 A")
    p["UNITS"] = ("DN", "level 1.5, the mean of a bin")
    p["CREATOR"] = (
        f"aia-94-leak-paper {importlib.metadata.version('aia-94-leak-paper')}",
        "",
    )
    hdus = [primary]
    for name, value in (("IMAGES", images), ("FLAT", flat)):
        hdus.append(
            fits.CompImageHDU(
                value.astype(np.float32),
                name=name,
                compression_type="GZIP_2",
                quantize_level=0.0,
            )
        )

    output.mkdir(parents=True, exist_ok=True)
    path = output / f"aia_dems_{date:%Y%m%d}.fits"
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path
