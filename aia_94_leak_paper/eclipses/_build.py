"""
Build the archive of eclipse ingresses from the AIA level-1 images at JSOC.

For each day in :data:`dates`, the ingress is found from the disk means of
the 171 Å channel, and the 94 Å and 171 Å images of three windows are
fetched: a few frames eight minutes before the ingress, every frame of the
ingress, and a few frames in the umbra, which are the dark reference. The
94 Å images are multiplied by the flat field they were divided by, which
returns them to the units of the detector, and every image is averaged over
square bins. Each day is written to a FITS file of its own.
"""

import datetime
import importlib.metadata
import pathlib
import time
import urllib.request

import numpy as np
from astropy.io import fits

__all__ = [
    "binning",
    "build",
    "channels",
    "ingress",
]

binning: int = 16
"""The side of the square bins the images are averaged over, in pixels."""

channels: tuple[int, ...] = (94, 171)
"""The wavelengths of the channels in the archive, in angstroms."""

_jsoc = "http://jsoc.stanford.edu"

_series = "aia.lev1_euv_12s"


def _retry(function, what: str, attempts: int = 5):
    """Call a function, retrying when JSOC times out or refuses."""
    for attempt in range(attempts):
        try:
            return function()
        # JSOC fails in many ways: timeouts, refused connections, HTTP errors,
        # and the query errors of `drms`, all of which are worth a retry.
        except Exception as e:  # noqa: BLE001
            print(f"{what}: {e!r}, retrying", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{what} failed {attempts} times")


def _client():
    import drms

    return drms.Client()


def ingress(date: datetime.date) -> tuple[np.datetime64, np.datetime64]:
    """
    When the 171 Å disk mean starts to fall on a day, and when it is gone.

    The disk mean is divided by the exposure time, since the automatic
    exposure control shortens the exposures in a flare, which would
    otherwise look like the start of an eclipse.

    Parameters
    ----------
    date
        A day of an eclipse season.
    """
    query = f"{_series}[{date:%Y-%m-%d}T05:30:00Z/180m][171]"
    keys = _retry(
        lambda: _client().query(query, key="T_OBS,DATAMEAN,EXPTIME"),
        what=query,
    )
    keys = keys[keys.EXPTIME.astype(float) > 0.1]
    t = np.array([np.datetime64(s.rstrip("Z")) for s in keys.T_OBS])
    rate = (
        keys.DATAMEAN.astype(float).to_numpy() / keys.EXPTIME.astype(float).to_numpy()
    )
    order = np.argsort(t)
    t, rate = t[order], rate[order]
    before = np.median(rate[:60])
    falling = np.flatnonzero(rate < 0.97 * before)
    if not falling.size:
        raise ValueError(f"no eclipse on {date}")
    start = falling[0]
    gone = start + np.flatnonzero(rate[start:] < 0.01 * before)[0]
    return t[start], t[gone]


def _records(start: np.datetime64, seconds: int, channel: int) -> list[dict]:
    """The images of a channel observed in a window, with the path of each at JSOC."""
    query = (
        f"{_series}[{np.datetime_as_string(start, unit='s')}Z/{seconds}s][{channel}]"
    )
    keys, segments = _retry(
        lambda: _client().query(
            query,
            key="T_OBS,WAVELNTH,EXPTIME,FLAT_REC,FSN,QUALITY",
            seg="image",
        ),
        what=query,
    )
    result = []
    for (_, k), (_, s) in zip(keys.iterrows(), segments.iterrows()):
        if int(k.WAVELNTH) != channel:
            raise ValueError(f"{query} returned {k.WAVELNTH} A")
        if float(k.EXPTIME) < 0.1 or "SUM" not in str(s.image):
            continue
        result.append(dict(k) | {"path": s.image})
    return result


def _fetch(url: str, path: pathlib.Path) -> pathlib.Path:
    """A file, downloaded unless it is already there."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        data = _retry(
            lambda: urllib.request.urlopen(url, timeout=120).read(),
            what=url,
        )
        path.with_suffix(".part").write_bytes(data)
        path.with_suffix(".part").replace(path)
    return path


def _flat(record: str, flats: pathlib.Path) -> np.ndarray:
    """The flat field named by the ``FLAT_REC`` keyword of an image."""
    path = flats / record.replace("aia.flatfield[:#", "flat_").replace("]", ".fits")
    if not path.exists():
        _, segments = _retry(
            lambda: _client().query(record, key="WAVELNTH", seg="flatfield"),
            what=record,
        )
        _fetch(_jsoc + segments.flatfield.iloc[0], path)
    return fits.getdata(path).astype(np.float32)


def build(
    date: datetime.date,
    cache: pathlib.Path,
    flats: pathlib.Path,
    output: pathlib.Path,
) -> pathlib.Path:
    """
    Write the archive file of one day and return its path.

    Parameters
    ----------
    date
        A day of :data:`dates`.
    cache
        Where the level-1 images are kept, one directory per day. A file of
        the right name anywhere in the directory of its day is used rather
        than fetched again.
    flats
        Where the flat fields are kept.
    output
        Where to write the archive file.
    """
    start, gone = ingress(date)
    windows = {
        "PRE": (start - np.timedelta64(8, "m"), 36),
        "INGRESS": (
            start - np.timedelta64(24, "s"),
            int((gone - start) / np.timedelta64(1, "s")) + 114,
        ),
        "UMBRA": (gone + np.timedelta64(8, "m"), 36),
    }
    day = cache / f"{date:%Y%m%d}"
    have = {p.name: p for p in day.rglob("aia*.fits")} if day.exists() else {}

    hdus = []
    header = None
    for channel in channels:
        for group, (begin, seconds) in windows.items():
            records = sorted(
                _records(begin, seconds, channel), key=lambda r: r["T_OBS"]
            )
            images, rows = [], []
            for record in records:
                stamp = record["T_OBS"][:19].replace(":", "").replace("-", "")
                name = f"aia{channel}_{stamp}.fits"
                path = have.get(name) or _fetch(
                    _jsoc + record["path"], day / str(channel) / name
                )
                with fits.open(path) as hdul:
                    h = hdul[1].header
                    image = hdul[1].data.astype(np.float32)
                if channel == 94:
                    image = image * _flat(h["FLAT_REC"], flats)
                if header is None:
                    header = h
                n = 4096 // binning
                images.append(image.reshape(n, binning, n, binning).mean(axis=(1, 3)))
                rows.append(
                    (
                        h["T_OBS"],
                        h["EXPTIME"],
                        h["FSN"],
                        h["QUALITY"],
                        h["CRPIX1"],
                        h["CRPIX2"],
                        h["FLAT_REC"],
                    )
                )
            name = f"{group}{channel}"
            hdus.append(
                fits.CompImageHDU(
                    np.array(images),
                    name=name,
                    compression_type="GZIP_2",
                    quantize_level=0.0,
                )
            )
            columns = list(zip(*rows))
            hdus.append(
                fits.BinTableHDU.from_columns(
                    [
                        fits.Column("T_OBS", "26A", array=columns[0]),
                        fits.Column("EXPTIME", "E", unit="s", array=columns[1]),
                        fits.Column("FSN", "J", array=columns[2]),
                        fits.Column("QUALITY", "J", array=columns[3]),
                        fits.Column("CRPIX1", "E", array=columns[4]),
                        fits.Column("CRPIX2", "E", array=columns[5]),
                        fits.Column("FLAT_REC", "24A", array=columns[6]),
                    ],
                    name=f"{name}_FRAMES",
                )
            )

    primary = fits.PrimaryHDU()
    p = primary.header
    p["DATE_ECL"] = (f"{date:%Y-%m-%d}", "the day of the eclipse")
    p["T_START"] = (
        np.datetime_as_string(start, unit="ms"),
        "171 disk mean begins to fall",
    )
    p["T_GONE"] = (
        np.datetime_as_string(gone, unit="ms"),
        "171 disk mean under 1% of before",
    )
    p["BINNING"] = (binning, "pixels averaged along each side of a bin")
    for key in (
        "CRPIX1",
        "CRPIX2",
        "CDELT1",
        "CDELT2",
        "R_SUN",
        "RSUN_OBS",
        "DSUN_OBS",
    ):
        p[key] = (header[key], "of the first 94 frame before the eclipse")
    p["UNITS94"] = ("DN", "94: level 1 times its flat field")
    p["UNITS171"] = ("DN", "171: level 1")
    p["CREATOR"] = (
        f"aia-94-leak-paper {importlib.metadata.version('aia-94-leak-paper')}",
        "",
    )

    output.mkdir(parents=True, exist_ok=True)
    path = output / f"aia94_eclipse_{date:%Y%m%d}.fits"
    fits.HDUList([primary] + hdus).writeto(path, overwrite=True)
    return path
