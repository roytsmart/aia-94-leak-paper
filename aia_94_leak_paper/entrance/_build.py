"""
Build the archive of the images AIA takes without a focal-plane filter,
from the level-1 images at JSOC.

The entrance filters of AIA are monitored for leakage of visible light by
images taken with the filter wheel at its open position, without a
focal-plane filter (Boerner et al. 2012). The telescope that holds the
94 Å and 304 Å channels names these images :data:`wave_str`. A set of them
is taken a little after 20:00 UT, every week early in the mission and every
two weeks later; for each month of the mission, the first image of the
month with a short exposure is fetched, returned to the units of the
detector by multiplying it by the flat field it was divided by, divided by
its exposure time, and averaged over square bins. Every pixel is kept, also
where the disk is brighter than the camera is linear to, and a bin is blank
only where a pixel is missing. Every month is written to one FITS file.
"""

import datetime
import importlib.metadata
import json
import pathlib

import numpy as np
from astropy.io import fits

from .. import _jsoc

__all__ = [
    "binning",
    "build",
    "months",
    "wave_str",
]

wave_str: str = "304_OPEN"
"""The name, in ``WAVE_STR``, of the images of the 94 Å and 304 Å telescope without a focal-plane filter."""

binning: int = 16
"""The side of the square bins the images are averaged over, in pixels."""

_series = "aia.lev1"
"""The series that holds every AIA level-1 image, calibration images among them."""

_exposure = (0.003, 0.006)
"""The range of the short exposures, in seconds; the longer ones saturate the disk."""

_evening = (datetime.time(20), 45 * 60)
"""When on each day the sets are searched for, and for how many seconds."""

_first = datetime.date(2011, 12, 1)
"""The first month with these images; the first set was taken on 2011 December 6."""


def months(last: datetime.date) -> list[datetime.date]:
    """
    The first day of every month from :data:`_first` to a day.

    Parameters
    ----------
    last
        A day of the last month.
    """
    result = []
    month = _first
    while month <= last:
        result.append(month)
        month = (month + datetime.timedelta(days=32)).replace(day=1)
    return result


def _records(day: datetime.date) -> list[dict]:
    """The short images without a focal-plane filter on the evening of a day, sorted by time."""
    start, seconds = _evening
    t = datetime.datetime.combine(day, start)
    query = (
        f"{_series}[{t:%Y-%m-%dT%H:%M:%S}Z/{seconds}s]"
        f"[? WAVE_STR='{wave_str}' AND IMG_TYPE='LIGHT' ?]"
    )
    keys, paths = _jsoc.retry(
        lambda: _jsoc.client().query(
            query,
            key="T_OBS,WAVE_STR,EXPTIME,FLAT_REC,FSN,QUALITY",
            seg="image_lev1",
        ),
        what=query,
    )
    result = []
    for (_, k), (_, s) in zip(keys.iterrows(), paths.iterrows()):
        exposure = float(k.EXPTIME)
        if _exposure[0] < exposure < _exposure[1] and "SUM" in str(s.image_lev1):
            result.append(dict(k) | {"image_lev1": s.image_lev1})
    return sorted(result, key=lambda r: r["T_OBS"])


def _month(month: datetime.date) -> None | dict:
    """The first short image without a focal-plane filter in a month, if any."""
    day = month
    while day.month == month.month:
        records = _records(day)
        if records:
            return records[0]
        day += datetime.timedelta(days=1)
    return None


def _reduce(path: pathlib.Path, flats: pathlib.Path) -> tuple[np.ndarray, dict]:
    """
    An image per second in the units of the detector, averaged over bins,
    and the keywords of its header the archive keeps.
    """
    with fits.open(path) as hdul:
        header = hdul[1].header
        image = hdul[1].data.astype(np.float32)
    raw = image * _jsoc.flat(header["FLAT_REC"], flats)
    n = 4096 // binning
    blocks = raw.reshape(n, binning, n, binning)
    good = np.isfinite(blocks)
    with np.errstate(invalid="ignore"):
        mean = np.where(good, blocks, 0).sum(axis=(1, 3)) / good.sum(axis=(1, 3))
    mean = np.where(good.all(axis=(1, 3)), mean, np.nan)
    keys = {
        "T_OBS": header["T_OBS"],
        "EXPTIME": header["EXPTIME"],
        "FSN": header["FSN"],
        "FLAT_REC": header["FLAT_REC"],
        "CRPIX1": header["CRPIX1"],
        "CRPIX2": header["CRPIX2"],
        "R_SUN": header["R_SUN"],
        "DATAMAX": header["DATAMAX"],
        # the pixels above 15000 DN of level 1, which AIA counts as saturated
        "NSATPIX": header["NSATPIX"],
    }
    return (mean / header["EXPTIME"]).astype(np.float32), keys


def build(
    last: datetime.date,
    cache: pathlib.Path,
    flats: pathlib.Path,
    output: pathlib.Path,
) -> pathlib.Path:
    """
    Write the archive file, of every month up to a day, and return its path.

    Parameters
    ----------
    last
        A day of the last month.
    cache
        Where the level-1 images are kept.
    flats
        Where the flat fields are kept.
    output
        Where to write the archive file.
    """
    # the image found for each month is kept with the images, so that a build
    # broken off by JSOC need not search the months it has done again
    found = cache / "months.json"
    index = json.loads(found.read_text()) if found.exists() else {}
    images = []
    rows = []
    for month in months(last):
        key = f"{month:%Y-%m}"
        if key not in index:
            index[key] = _month(month)
            cache.mkdir(parents=True, exist_ok=True)
            found.write_text(json.dumps(index, indent=1, default=str))
        record = index[key]
        if record is None:
            print(f"{month:%Y-%m}: no image", flush=True)
            continue
        stamp = record["T_OBS"][:19].replace(":", "").replace("-", "")
        path = _jsoc.fetch(record["image_lev1"], cache / f"open304_{stamp}.fits")
        image, keys = _reduce(path, flats)
        images.append(image)
        rows.append(keys)
        print(f"{month:%Y-%m}: {keys['T_OBS']}", flush=True)

    primary = fits.PrimaryHDU()
    p = primary.header
    p["WAVE_STR"] = (wave_str, "every image is without a focal-plane filter")
    p["BINNING"] = (binning, "pixels along each side of a bin")
    p["UNITS"] = ("DN/s", "level 1 times its flat field, per second")
    p["CREATOR"] = (
        f"aia-94-leak-paper {importlib.metadata.version('aia-94-leak-paper')}",
        "",
    )
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="T_OBS", format="24A", array=[r["T_OBS"] for r in rows]),
            fits.Column(name="EXPTIME", format="D", array=[r["EXPTIME"] for r in rows]),
            fits.Column(name="FSN", format="K", array=[r["FSN"] for r in rows]),
            fits.Column(
                name="FLAT_REC", format="24A", array=[r["FLAT_REC"] for r in rows]
            ),
            fits.Column(name="CRPIX1", format="D", array=[r["CRPIX1"] for r in rows]),
            fits.Column(name="CRPIX2", format="D", array=[r["CRPIX2"] for r in rows]),
            fits.Column(name="R_SUN", format="D", array=[r["R_SUN"] for r in rows]),
            fits.Column(name="DATAMAX", format="D", array=[r["DATAMAX"] for r in rows]),
            fits.Column(name="NSATPIX", format="K", array=[r["NSATPIX"] for r in rows]),
        ],
        name="FRAMES",
    )
    hdus = [
        primary,
        fits.CompImageHDU(
            np.array(images),
            name="IMAGES",
            compression_type="GZIP_2",
            quantize_level=0.0,
        ),
        table,
    ]
    output.mkdir(parents=True, exist_ok=True)
    path = output / "aia_entrance.fits"
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path
