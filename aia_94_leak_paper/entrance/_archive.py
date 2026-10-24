"""
Fetch the archive of the images AIA takes without a focal-plane filter and
read it.

The archive is built by :func:`aia_94_leak_paper.entrance.build` and
published as an asset of a release of this repository.
"""

import dataclasses
import hashlib
import pathlib
import tarfile
import urllib.request

import astropy.time
import astropy.units as u
import named_arrays as na
import numpy as np
from astropy.io import fits

from ..eclipses._archive import directory_default
from ._build import binning

__all__ = [
    "Monitor",
    "archive",
    "load",
    "sha256",
    "url",
    "version",
]

version: str = "entrance-v2"
"""The release of this repository the archive is an asset of."""

url: str = (
    "https://github.com/roytsmart/aia-94-leak-paper/releases/download/"
    f"{version}/aia-entrance.tar"
)
"""Where the archive is published."""

sha256: str = "4654597f5dfc821b4bb1afd810572de53c1a346628d3a7f9c3e75688f0520fa8"
"""The SHA-256 digest of the archive, checked after every download."""


def archive(directory: None | pathlib.Path = None) -> pathlib.Path:
    """
    The directory of the extracted archive, downloaded and checked if needed.

    Parameters
    ----------
    directory
        Where the archive is kept,
        :data:`aia_94_leak_paper.eclipses.directory_default` if :obj:`None`.
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
class Monitor:
    """
    One image a month of the 94 Å and 304 Å telescope without a focal-plane
    filter, in the units of the detector per second, averaged over bins.

    Every array has the axis ``time``; the images also have ``y`` and ``x``.
    """

    time: astropy.time.Time
    """When each image was taken."""

    exposure: u.Quantity
    """The exposure time of each image."""

    images: na.ScalarArray
    """The mean of each bin per second, blank where a pixel of the bin is missing."""

    center: na.Cartesian2dVectorArray
    """The center of the solar disk in each image, in bins from the corner of the first."""

    radius: na.ScalarArray
    """The radius of the photosphere in each image, in bins."""

    @property
    def position(self) -> na.Cartesian2dVectorArray:
        """The center of each bin, in solar radii from the center of the disk."""
        n = self.images.shape["x"]
        index = na.Cartesian2dVectorArray(
            x=na.arange(0, n, axis="x"),
            y=na.arange(0, n, axis="y"),
        )
        return (index - self.center) / self.radius

    @property
    def r(self) -> na.ScalarArray:
        """The distance of each bin from the center of the disk, in solar radii."""
        position = self.position
        return np.sqrt(np.square(position.x) + np.square(position.y))


def load(directory: None | pathlib.Path = None) -> Monitor:
    """
    Read the archive.

    Parameters
    ----------
    directory
        Where the archive is kept,
        :data:`aia_94_leak_paper.eclipses.directory_default` if :obj:`None`.
    """
    path = archive(directory) / "aia_entrance.fits"
    with fits.open(path) as hdul:
        images = np.array(hdul["IMAGES"].data, dtype=float)
        frames = hdul["FRAMES"].data
        t = [s.strip().rstrip("Z") for s in frames["T_OBS"]]
        exposure = np.array(frames["EXPTIME"], dtype=float)
        crpix1 = np.array(frames["CRPIX1"], dtype=float)
        crpix2 = np.array(frames["CRPIX2"], dtype=float)
        r_sun = np.array(frames["R_SUN"], dtype=float)
    return Monitor(
        time=astropy.time.Time(t, scale="utc"),
        exposure=na.ScalarArray(exposure * u.s, axes="time"),
        images=na.ScalarArray(images * u.DN / u.s, axes=("time", "y", "x")),
        center=na.Cartesian2dVectorArray(
            x=na.ScalarArray(((crpix1 - 1) - (binning - 1) / 2) / binning, "time"),
            y=na.ScalarArray(((crpix2 - 1) - (binning - 1) / 2) / binning, "time"),
        ),
        radius=na.ScalarArray(r_sun / binning, axes="time"),
    )
