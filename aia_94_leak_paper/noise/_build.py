"""
Build the archive of the noise of AIA 94 Å from the level-1 images at JSOC.

For each day of :data:`aia_94_leak_paper.eclipses.dates`, a run of
consecutive 94 Å frames is fetched at noon, hours away from the eclipse,
and returned to the units of the detector by multiplying each by the flat
field it was divided by. For every three consecutive frames the residual

.. math::

    d = I_k - \\tfrac{1}{2} (I_{k-1} + I_{k+1})

cancels the solar signal and anything that changes linearly with time,
leaving noise with 3/2 the variance of one frame. Each pixel's residual is
sorted into a bin of signal by the mean of the *other* frames of the run,
whose noise is independent of the residual, and the signal itself is
estimated from the three frames. The residuals are accumulated, over the
pixels of square regions and over the triplets of the run, into the count,
the sums of the signal and of the residual, and the sum of the squared
residual of every region and bin of signal. Pixels the despiker changed in
any frame are left out, and so are residuals beyond eight robust standard
deviations of all the residuals of their bin.
"""

import datetime
import importlib.metadata
import pathlib

import numpy as np
from astropy.io import fits

from .. import _jsoc

__all__ = [
    "accumulate",
    "build",
    "edges",
    "frames",
    "region",
    "usable",
]

frames: int = 15
"""The number of consecutive frames in a run."""

region: int = 64
"""The side of the square regions the statistics are accumulated over, in pixels."""

edges: np.ndarray = np.concatenate([[-1e3], np.geomspace(0.5, 64, 22), [1e5]])
"""The edges of the bins of signal, in DN of the detector per frame."""

_noon = datetime.time(12)
"""When on each day the run begins, hours away from the eclipse."""

_clip = 8
"""How many robust standard deviations of its bin a residual may lie from zero."""

_harmless = 0b100
"""
The bits of ``QUALITY`` a frame of a run may have set: only that its
ancillary science data are not available, which changes no pixel, and which
is set on every frame of 2026 January 19.
"""


def usable(records: list[dict], filter: str = "94_THIN", allowed: int = 0) -> bool:
    """
    Whether a run of 94 Å frames has equal exposures, all through one
    focal-plane filter, and no flags in ``QUALITY`` but :data:`_harmless`.

    Parameters
    ----------
    records
        The frames of the run, as :func:`aia_94_leak_paper._jsoc.records`
        gives them.
    filter
        The focal-plane filter every frame must have been taken through, as
        ``WAVE_STR`` names it.
    allowed
        Further bits of ``QUALITY`` the frames may have set.
    """
    exposures = {round(float(r["EXPTIME"]), 3) for r in records}
    filters = {r["WAVE_STR"] for r in records}
    quality = {int(r["QUALITY"]) & ~(_harmless | allowed) for r in records}
    return len(exposures) == 1 and filters == {filter} and quality == {0}


def _run(date: datetime.date) -> list[dict]:
    """The first :func:`usable` run of :data:`frames` consecutive 94 Å frames from noon on a day."""
    start = np.datetime64(datetime.datetime.combine(date, _noon))
    for attempt in range(20):
        records = _jsoc.records(start, 12 * frames + 6, 94, ("image", "spikes"))
        if len(records) >= frames:
            records = records[:frames]
            if usable(records):
                return records
        start = start + np.timedelta64(5, "m")
    raise ValueError(f"no clean run of 94 frames near noon on {date}")


def accumulate(
    records: list[dict],
    directory: pathlib.Path,
    flats: pathlib.Path,
) -> fits.HDUList:
    """
    The noise of a run of consecutive 94 Å frames, accumulated over regions
    and bins of signal, as the file of the archive holds it.

    Parameters
    ----------
    records
        The frames of the run, as :func:`aia_94_leak_paper._jsoc.records`
        gives them with the ``image`` and ``spikes`` segments.
    directory
        Where the level-1 images of the run and their despike records are kept.
    flats
        Where the flat fields are kept.
    """
    n_frames = len(records)
    images = []
    spiked = np.zeros((4096, 4096), dtype=bool)
    for i, record in enumerate(records):
        stamp = record["T_OBS"][:19].replace(":", "").replace("-", "")
        image = _jsoc.fetch(record["image"], directory / f"aia94_{stamp}.fits")
        spikes = _jsoc.fetch(record["spikes"], directory / f"spikes94_{stamp}.fits")
        with fits.open(image) as hdul:
            header = hdul[1].header
            images.append(
                hdul[1].data.astype(np.float32) * _jsoc.flat(header["FLAT_REC"], flats)
            )
        with fits.open(spikes) as hdul:
            spiked.flat[hdul[1].data[0].astype(np.int64)] = True
        if i == 0:
            first = header
    images = np.array(images)

    n = 4096 // region
    y, x = np.indices((4096, 4096))
    cell = ((y // region) * n + (x // region)).ravel()
    total = images.sum(axis=0, dtype=np.float64)
    nbins = len(edges) - 1
    shape = (nbins, n * n)
    count = np.zeros(shape)
    signal = np.zeros(shape)
    residual = np.zeros(shape)
    squared = np.zeros(shape)
    keep = ~spiked.ravel()

    triplets = range(1, n_frames - 1)
    samples = []
    for k in triplets:
        d = images[k] - (images[k - 1] + images[k + 1]) / 2
        s = (images[k - 1] + images[k] + images[k + 1]) / 3
        other = (total - images[k - 1] - images[k] - images[k + 1]) / (n_frames - 3)
        b = np.clip(np.digitize(other.ravel(), edges) - 1, 0, nbins - 1)
        samples.append((d.ravel(), s.ravel(), b))

    # the threshold of each bin, from the residuals of every region in it
    threshold = np.zeros(nbins)
    for i in range(nbins):
        d = np.concatenate([d[(b == i) & keep] for d, _, b in samples])
        if d.size:
            mad = np.median(np.abs(d - np.median(d)))
            threshold[i] = _clip * 1.4826 * mad

    for d, s, b in samples:
        use = keep & (np.abs(d) <= threshold[b])
        index = b[use] * n * n + cell[use]
        size = nbins * n * n
        count += np.bincount(index, minlength=size).reshape(shape)
        signal += np.bincount(index, weights=s[use], minlength=size).reshape(shape)
        residual += np.bincount(index, weights=d[use], minlength=size).reshape(shape)
        squared += np.bincount(index, weights=d[use] ** 2, minlength=size).reshape(
            shape
        )

    primary = fits.PrimaryHDU()
    p = primary.header
    p["DATE_RUN"] = (records[0]["T_OBS"][:10], "the day of the run")
    p["T_FIRST"] = (records[0]["T_OBS"], "T_OBS of the first frame")
    p["T_LAST"] = (records[-1]["T_OBS"], "T_OBS of the last frame")
    p["NFRAMES"] = (n_frames, "consecutive frames in the run")
    p["EXPTIME"] = (float(records[0]["EXPTIME"]), "exposure of every frame, s")
    p["REGION"] = (region, "pixels along each side of a region")
    p["CLIP"] = (_clip, "robust standard deviations kept")
    p["NSPIKED"] = (int(spiked.sum()), "pixels the despiker changed in any frame")
    for key in ("CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "R_SUN", "RSUN_OBS"):
        p[key] = (first[key], "of the first frame")
    p["FLAT_REC"] = (first["FLAT_REC"], "the flat field each frame is times")
    p["UNITS"] = ("DN", "level 1 times its flat field, per frame")
    p["CREATOR"] = (
        f"aia-94-leak-paper {importlib.metadata.version('aia-94-leak-paper')}",
        "",
    )
    hdus = [primary, fits.ImageHDU(edges, name="EDGES")]
    for name, value in (
        ("COUNT", count),
        ("SIGNAL", signal),
        ("RESIDUAL", residual),
        ("SQUARED", squared),
    ):
        hdus.append(
            fits.ImageHDU(value.reshape(nbins, n, n).astype(np.float32), name=name)
        )
    return fits.HDUList(hdus)


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
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    cache
        Where the level-1 images and their despike records are kept, one
        directory per day.
    flats
        Where the flat fields are kept.
    output
        Where to write the archive file.
    """
    hdul = accumulate(_run(date), cache / f"{date:%Y%m%d}", flats)
    output.mkdir(parents=True, exist_ok=True)
    path = output / f"aia94_noise_{date:%Y%m%d}.fits"
    hdul.writeto(path, overwrite=True)
    return path
