"""
Fetch the archive of the six coronal channels of AIA and read a day from it.

The archive is built by :func:`aia_94_leak_paper.dems.build` and published
as an asset of a release of this repository.
"""

import dataclasses
import datetime
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
from ._build import wavelengths

__all__ = [
    "Synoptic",
    "archive",
    "load",
    "sha256",
    "url",
    "version",
]

version: str = "dems-v1"
"""The release of this repository the archive is an asset of."""

url: str = (
    "https://github.com/roytsmart/aia-94-leak-paper/releases/download/"
    f"{version}/aia-dems.tar"
)
"""Where the archive is published."""

sha256: str = "e572f09f9dba5e4c0d0a81f06ad0c50d932328fa0068c368a69f98d72f843703"
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
class Synoptic:
    """
    The six coronal channels on one day, registered to one another and to
    solar north, with the disk at the center.

    The images have the axes ``wavelength``, ``y``, and ``x``, and are in
    the units of level 1.5 per bin.
    """

    date: datetime.date
    """The day."""

    time: astropy.time.Time
    """When the 94 Å image was observed."""

    wavelength: na.ScalarArray
    """The channel of each image, along ``wavelength``."""

    exposure: na.ScalarArray
    """The exposure time of each image, along ``wavelength``."""

    images: na.ScalarArray
    """The mean of each bin."""

    flat: na.ScalarArray
    """
    The flat field the 94 Å image was divided by, averaged over the detector
    pixels of each bin, with the axes ``y`` and ``x``.
    """

    samples: int
    """The level-1.5 pixels in the mean of a bin."""

    scale: u.Quantity
    """The side of a bin."""

    radius: float
    """The radius of the photosphere, in bins."""

    @property
    def position(self) -> na.Cartesian2dVectorArray:
        """The center of each bin, in solar radii from the center of the disk."""
        n = self.images.shape["x"]
        index = na.Cartesian2dVectorArray(
            x=na.arange(0, n, axis="x"),
            y=na.arange(0, n, axis="y"),
        )
        return (index - (n - 1) / 2) / self.radius

    @property
    def r(self) -> na.ScalarArray:
        """The distance of each bin from the center of the disk, in solar radii."""
        position = self.position
        return np.sqrt(np.square(position.x) + np.square(position.y))

    @property
    def rate(self) -> na.ScalarArray:
        """The images per second."""
        return self.images / self.exposure


def load(date: datetime.date, directory: None | pathlib.Path = None) -> Synoptic:
    """
    Read one day from the archive.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    directory
        Where the archive is kept,
        :data:`aia_94_leak_paper.eclipses.directory_default` if :obj:`None`.
    """
    path = archive(directory) / f"aia_dems_{date:%Y%m%d}.fits"
    with fits.open(path) as hdul:
        p = hdul[0].header
        images = np.array(hdul["IMAGES"].data, dtype=float)
        flat = np.array(hdul["FLAT"].data, dtype=float)
    return Synoptic(
        date=date,
        time=astropy.time.Time(p["T_0094"].rstrip("Z"), scale="utc"),
        wavelength=na.ScalarArray(np.array(wavelengths) * u.AA, axes="wavelength"),
        exposure=na.ScalarArray(
            np.array([p[f"X_{w:04d}"] for w in wavelengths]) * u.s, axes="wavelength"
        ),
        images=na.ScalarArray(images * u.DN, axes=("wavelength", "y", "x")),
        flat=na.ScalarArray(flat, axes=("y", "x")),
        samples=p["NSAMPLE"],
        scale=p["CDELT"] * u.arcsec,
        radius=p["RSUN_OBS"] / p["CDELT"],
    )
