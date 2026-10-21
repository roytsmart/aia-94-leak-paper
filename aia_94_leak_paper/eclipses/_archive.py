"""
Fetch the archive of eclipse ingresses and read a day from it.

The archive is built by :func:`aia_94_leak_paper.eclipses.build` and
published as an asset of a release of this repository, so that building the
article needs neither JSOC nor the level-1 images.
"""

import dataclasses
import datetime
import hashlib
import os
import pathlib
import tarfile
import urllib.request

import astropy.time
import astropy.units as u
import named_arrays as na
import numpy as np
from astropy.io import fits

__all__ = [
    "Eclipse",
    "archive",
    "directory_default",
    "load",
    "sha256",
    "url",
    "version",
]

version: str = "eclipses-v1"
"""The release of this repository the archive is an asset of."""

url: str = (
    "https://github.com/roytsmart/aia-94-leak-paper/releases/download/"
    f"{version}/aia-94-eclipses.tar"
)
"""Where the archive is published."""

sha256: str = ""
"""The SHA-256 digest of the archive, checked after every download."""

directory_default = pathlib.Path(
    os.environ.get(
        "AIA_94_LEAK_CACHE",
        pathlib.Path.home() / ".cache" / "aia_94_leak_paper",
    )
)
"""
Where the archive is kept, ``~/.cache/aia_94_leak_paper`` unless the
``AIA_94_LEAK_CACHE`` environment variable says otherwise.
"""


def archive(directory: None | pathlib.Path = None) -> pathlib.Path:
    """
    The directory of the extracted archive, downloaded and checked if needed.

    Parameters
    ----------
    directory
        Where the archive is kept, :data:`directory_default` if :obj:`None`.
    """
    if directory is None:
        directory = directory_default
    result = pathlib.Path(directory) / version
    if result.exists():
        return result
    result.parent.mkdir(parents=True, exist_ok=True)
    tar = result.parent / f"{version}.tar"
    if not tar.exists():
        part = tar.with_suffix(".part")
        with (
            urllib.request.urlopen(url, timeout=600) as response,
            open(part, "wb") as f,
        ):
            while chunk := response.read(1 << 20):
                f.write(chunk)
        part.replace(tar)
    digest = hashlib.sha256(tar.read_bytes()).hexdigest()
    if digest != sha256:
        tar.unlink()
        raise ValueError(f"{url} has SHA-256 {digest}, not {sha256}")
    with tarfile.open(tar) as t:
        t.extractall(result.with_suffix(".part"), filter="data")
    result.with_suffix(".part").replace(result)
    return result


@dataclasses.dataclass
class Eclipse:
    """
    The ingress of one eclipse, averaged over square bins.

    Every image has the axes ``time``, ``y``, and ``x``. The 94 Å images are
    in the units of the detector (level 1 times its flat field) and the
    171 Å images are level 1; neither is divided by the exposure time.
    """

    date: datetime.date
    """The day of the eclipse."""

    images: dict[str, na.ScalarArray]
    """
    The images of each window and channel, keyed ``PRE94``, ``INGRESS94``,
    ``UMBRA94``, ``PRE171``, ``INGRESS171``, and ``UMBRA171``.
    """

    times: dict[str, na.ScalarArray]
    """The time each image was observed (``T_OBS``), keyed as :attr:`images`."""

    exposures: dict[str, na.ScalarArray]
    """The exposure time of each image, keyed as :attr:`images`."""

    start: astropy.time.Time
    """When the 171 Å disk mean began to fall."""

    gone: astropy.time.Time
    """When the 171 Å disk mean fell below a hundredth of what it was before."""

    binning: int
    """The side of the square bins, in pixels."""

    center: na.Cartesian2dVectorArray
    """The center of the solar disk, in bins from the corner of the first."""

    radius: float
    """The radius of the photosphere, in bins."""

    @property
    def position(self) -> na.Cartesian2dVectorArray:
        """The center of each bin, in solar radii from the center of the disk."""
        n = self.images["PRE94"].shape["x"]
        index = na.Cartesian2dVectorArray(
            x=na.arange(0, n, axis="x"),
            y=na.arange(0, n, axis="y"),
        )
        return (index - self.center) / self.radius


def load(date: datetime.date, directory: None | pathlib.Path = None) -> Eclipse:
    """
    Read one day from the archive.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    directory
        Where the archive is kept, :data:`directory_default` if :obj:`None`.
    """
    path = archive(directory) / f"aia94_eclipse_{date:%Y%m%d}.fits"
    images, times, exposures = {}, {}, {}
    with fits.open(path) as hdul:
        p = hdul[0].header
        binning = p["BINNING"]
        for hdu in hdul[1:]:
            if hdu.name.endswith("_FRAMES"):
                key = hdu.name.removesuffix("_FRAMES")
                t = [s.strip().rstrip("Z") for s in hdu.data["T_OBS"]]
                times[key] = na.ScalarArray(
                    astropy.time.Time(t, scale="utc"), axes="time"
                )
                exposures[key] = na.ScalarArray(
                    np.array(hdu.data["EXPTIME"], dtype=float) * u.s, axes="time"
                )
            else:
                images[hdu.name] = na.ScalarArray(
                    np.array(hdu.data, dtype=float) * u.DN, axes=("time", "y", "x")
                )
    return Eclipse(
        date=date,
        images=images,
        times=times,
        exposures=exposures,
        start=astropy.time.Time(p["T_START"], scale="utc"),
        gone=astropy.time.Time(p["T_GONE"], scale="utc"),
        binning=binning,
        center=na.Cartesian2dVectorArray(
            x=((p["CRPIX1"] - 1) - (binning - 1) / 2) / binning,
            y=((p["CRPIX2"] - 1) - (binning - 1) / 2) / binning,
        ),
        radius=p["R_SUN"] / binning,
    )
