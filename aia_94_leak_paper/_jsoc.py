"""
Fetch AIA level-1 images and flat fields from JSOC, for the builders of the
archives this article reads.
"""

import pathlib
import time
import urllib.request

import numpy as np
from astropy.io import fits

__all__ = [
    "fetch",
    "flat",
    "records",
    "retry",
    "series",
    "url",
]

url: str = "http://jsoc.stanford.edu"
"""Where JSOC serves the files of a record, by their path in SUMS."""

series: str = "aia.lev1_euv_12s"
"""The series of AIA level-1 EUV images."""


def retry(function, what: str, attempts: int = 5):
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


def client():
    """A client for the JSOC database, from the optional dependency `drms`."""
    import drms

    return drms.Client()


def records(
    start: np.datetime64,
    seconds: int,
    channel: int,
    segments: tuple[str, ...] = ("image",),
) -> list[dict]:
    """
    The images of a channel observed in a window, sorted by time, with the
    path in SUMS of each of their segments.

    Parameters
    ----------
    start
        The start of the window.
    seconds
        The length of the window.
    channel
        The wavelength of the channel, in angstroms.
    segments
        The segments of each record to give the paths of.
    """
    query = f"{series}[{np.datetime_as_string(start, unit='s')}Z/{seconds}s][{channel}]"
    keys, paths = retry(
        lambda: client().query(
            query,
            key="T_OBS,WAVELNTH,EXPTIME,FLAT_REC,FSN,QUALITY",
            seg=",".join(segments),
        ),
        what=query,
    )
    result = []
    for (_, k), (_, s) in zip(keys.iterrows(), paths.iterrows()):
        if int(k.WAVELNTH) != channel:
            raise ValueError(f"{query} returned {k.WAVELNTH} A")
        if float(k.EXPTIME) < 0.1 or "SUM" not in str(s[segments[0]]):
            continue
        result.append(dict(k) | {segment: s[segment] for segment in segments})
    return sorted(result, key=lambda r: r["T_OBS"])


def fetch(source: str, path: pathlib.Path) -> pathlib.Path:
    """A file at a path in SUMS, downloaded unless it is already there."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        data = retry(
            lambda: urllib.request.urlopen(url + source, timeout=120).read(),
            what=source,
        )
        path.with_suffix(".part").write_bytes(data)
        path.with_suffix(".part").replace(path)
    return path


def flat(record: str, flats: pathlib.Path) -> np.ndarray:
    """The flat field named by the ``FLAT_REC`` keyword of an image."""
    path = flats / record.replace("aia.flatfield[:#", "flat_").replace("]", ".fits")
    if not path.exists():
        _, paths = retry(
            lambda: client().query(record, key="WAVELNTH", seg="flatfield"),
            what=record,
        )
        fetch(paths.flatfield.iloc[0], path)
    return fits.getdata(path).astype(np.float32)
