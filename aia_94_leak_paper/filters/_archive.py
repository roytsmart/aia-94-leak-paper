"""
Fetch the archive of AIA 94 Å through each of its focal-plane filters and
read the run through the thin filter, or a frame through the thick one, from
it.

The archive is built by :func:`aia_94_leak_paper.filters.build` and
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

import aia_94_leak_paper

from ..eclipses._archive import directory_default

__all__ = [
    "Frame",
    "archive",
    "frame",
    "frames",
    "run",
    "sha256",
    "url",
    "version",
]

version: str = "filters-v2"
"""The release of this repository the archive is an asset of."""

url: str = (
    "https://github.com/roytsmart/aia-94-leak-paper/releases/download/"
    f"{version}/aia-94-filters.tar"
)
"""Where the archive is published."""

sha256: str = "28e29b0ec369e643ba055eea054f93acea2bea30e1bfe3263eb105a69bee8ba0"
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


frames: tuple[str, ...] = ("clean", "banded")
"""
The frames through the thick filter in the archive: the last before the
channel changed back to the thin filter, followed by a frame of 304 Å
through its thin filter, and the one before it, followed by a frame through
the thick filter of 304 Å, which has the band.
"""


@dataclasses.dataclass
class Frame:
    """One 94 Å frame, averaged over square regions, in the units of the detector."""

    time: astropy.time.Time
    """When the frame was observed."""

    exposure: u.Quantity
    """The exposure time of the frame."""

    followed: str
    """The ``WAVE_STR`` of the 304 Å frame after it, which the filter wheel moved to as it was read out."""

    image: na.ScalarArray
    """The mean signal of each region, with the axes ``y`` and ``x``."""


def run(directory: None | pathlib.Path = None) -> "aia_94_leak_paper.noise.Run":
    """
    Read the run through the thin filter from the archive.

    Parameters
    ----------
    directory
        Where the archive is kept,
        :data:`aia_94_leak_paper.eclipses.directory_default` if :obj:`None`.
    """
    path = archive(directory) / "aia94_noise_thin.fits"
    return aia_94_leak_paper.noise.read(path)


def frame(name: str, directory: None | pathlib.Path = None) -> Frame:
    """
    Read a frame through the thick filter from the archive.

    Parameters
    ----------
    name
        One of :data:`frames`.
    directory
        Where the archive is kept,
        :data:`aia_94_leak_paper.eclipses.directory_default` if :obj:`None`.
    """
    with fits.open(archive(directory) / f"aia94_frame_{name}.fits") as hdul:
        p = hdul[0].header
        image = np.array(hdul["IMAGE"].data, dtype=float)
    return Frame(
        time=astropy.time.Time(p["T_OBS"].rstrip("Z"), scale="utc"),
        exposure=p["EXPTIME"] * u.s,
        followed=p["FOLLOWED"],
        image=na.ScalarArray(image * u.DN, axes=("y", "x")),
    )
