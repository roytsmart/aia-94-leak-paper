"""
Fetch the archive of the noise of AIA 94 Å and read a day from it.

The archive is built by :func:`aia_94_leak_paper.noise.build` and published
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

__all__ = [
    "Run",
    "archive",
    "load",
    "read",
    "sha256",
    "url",
    "version",
]

version: str = "noise-v1"
"""The release of this repository the archive is an asset of."""

url: str = (
    "https://github.com/roytsmart/aia-94-leak-paper/releases/download/"
    f"{version}/aia-94-noise.tar"
)
"""Where the archive is published."""

sha256: str = "b92ed3ea7a145e6f7e8945dbf9c6b530439c1def546c7f008b2c8e9dcfaef539"
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
class Run:
    """
    The noise of a run of consecutive 94 Å frames, accumulated over square
    regions and bins of signal.

    Every array has the axes ``signal``, ``y``, and ``x``, and is in the
    units of the detector per frame.
    """

    date: datetime.date
    """The day of the run."""

    time: astropy.time.Time
    """When the first frame was observed."""

    exposure: u.Quantity
    """The exposure time of every frame."""

    edges: na.ScalarArray
    """The edges of the bins of signal, along ``signal_vertex``."""

    count: na.ScalarArray
    """How many residuals each region and bin holds."""

    signal: na.ScalarArray
    """The mean signal of each region and bin, from the three frames of each residual."""

    mean: na.ScalarArray
    """The mean residual of each region and bin."""

    variance: na.ScalarArray
    """
    The variance of a single frame in each region and bin: the variance of
    the residual divided by 3/2.
    """

    region: int
    """The side of a region, in pixels."""

    center: na.Cartesian2dVectorArray
    """The center of the solar disk, in regions from the corner of the first."""

    radius: float
    """The radius of the photosphere, in regions."""

    @property
    def position(self) -> na.Cartesian2dVectorArray:
        """The center of each region, in solar radii from the center of the disk."""
        n = self.count.shape["x"]
        index = na.Cartesian2dVectorArray(
            x=na.arange(0, n, axis="x"),
            y=na.arange(0, n, axis="y"),
        )
        return (index - self.center) / self.radius

    @property
    def r(self) -> na.ScalarArray:
        """The distance of each region from the center of the disk, in solar radii."""
        position = self.position
        return np.sqrt(np.square(position.x) + np.square(position.y))

    @property
    def image(self) -> na.ScalarArray:
        """The mean signal of each region, over all its bins of signal."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.nan_to_num(self.signal * self.count).sum(
                "signal"
            ) / self.count.sum("signal")


def read(path: pathlib.Path) -> Run:
    """
    Read a run from a file written by :func:`aia_94_leak_paper.noise.accumulate`.

    Parameters
    ----------
    path
        The file of the run.
    """
    with fits.open(path) as hdul:
        p = hdul[0].header
        edges = np.array(hdul["EDGES"].data, dtype=float)
        count, signal, residual, squared = (
            np.array(hdul[name].data, dtype=float)
            for name in ("COUNT", "SIGNAL", "RESIDUAL", "SQUARED")
        )
    axes = ("signal", "y", "x")
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = residual / count
        variance = (squared / count - np.square(mean)) / 1.5
        signal = signal / count
    region = p["REGION"]
    return Run(
        date=datetime.date.fromisoformat(p["DATE_RUN"]),
        time=astropy.time.Time(p["T_FIRST"].rstrip("Z"), scale="utc"),
        exposure=p["EXPTIME"] * u.s,
        edges=na.ScalarArray(edges * u.DN, axes="signal_vertex"),
        count=na.ScalarArray(count, axes=axes),
        signal=na.ScalarArray(signal * u.DN, axes=axes),
        mean=na.ScalarArray(mean * u.DN, axes=axes),
        variance=na.ScalarArray(variance * u.DN**2, axes=axes),
        region=region,
        center=na.Cartesian2dVectorArray(
            x=((p["CRPIX1"] - 1) - (region - 1) / 2) / region,
            y=((p["CRPIX2"] - 1) - (region - 1) / 2) / region,
        ),
        radius=p["R_SUN"] / region,
    )


def load(date: datetime.date, directory: None | pathlib.Path = None) -> Run:
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
    return read(archive(directory) / f"aia94_noise_{date:%Y%m%d}.fits")
