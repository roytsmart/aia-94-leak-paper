"""
Fetch the archive of the noise of AIA 94 Å through each of its focal-plane
filters and read a filter from it.

The archive is built by :func:`aia_94_leak_paper.filters.build` and
published as an asset of a release of this repository.
"""

import hashlib
import pathlib
import tarfile
import urllib.request

import aia_94_leak_paper

from ..eclipses._archive import directory_default

__all__ = [
    "archive",
    "load",
    "sha256",
    "url",
    "version",
]

version: str = "filters-v1"
"""The release of this repository the archive is an asset of."""

url: str = (
    "https://github.com/roytsmart/aia-94-leak-paper/releases/download/"
    f"{version}/aia-94-filters.tar"
)
"""Where the archive is published."""

sha256: str = "e9c8bc9245d2a9aeba80fd34544ded480d863d62a327c9b2ec5d37ea0c4a03ab"
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


def load(
    filter: str,
    directory: None | pathlib.Path = None,
) -> "aia_94_leak_paper.noise.Run":
    """
    Read the run through one focal-plane filter from the archive.

    Parameters
    ----------
    filter
        One of :data:`aia_94_leak_paper.filters.filters`.
    directory
        Where the archive is kept,
        :data:`aia_94_leak_paper.eclipses.directory_default` if :obj:`None`.
    """
    path = archive(directory) / f"aia94_noise_{filter}.fits"
    return aia_94_leak_paper.noise.read(path)
