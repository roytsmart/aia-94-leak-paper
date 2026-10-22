"""
Build the archive of the noise of AIA 94 Å through each of its focal-plane
filters from the level-1 images at JSOC.

The channel has two focal-plane filters, the thin one it normally observes
through and a thick one. On 2026 March 7 it observed through the thick one
until near :data:`switch`, when it changed back to the thin one. The last
run of consecutive frames through the thick filter before the change and
the first through the thin filter after it see nearly the same Sun, and
each is accumulated as :func:`aia_94_leak_paper.noise.accumulate` does a
run of the noise archive.
"""

import datetime
import pathlib

import numpy as np

from .. import _jsoc, noise

__all__ = [
    "build",
    "filters",
    "switch",
]

filters: tuple[str, ...] = ("thin", "thick")
"""The focal-plane filters of the channel, each named as in ``WAVE_STR``."""

switch: datetime.datetime = datetime.datetime(2026, 3, 7, 23, 20)
"""
When the channel changed from its thick focal-plane filter back to its thin
one; the frames of the few minutes before are flagged in ``QUALITY``.
"""

_window = datetime.timedelta(minutes=15)
"""How far from :data:`switch` the runs are searched for."""


def _consecutive(records: list[dict]) -> bool:
    """Whether a run of frames was taken at the cadence of the channel, 12 s, with none missing."""
    times = [np.datetime64(r["T_OBS"].rstrip("Z")) for r in records]
    span = (times[-1] - times[0]) / np.timedelta64(1, "s")
    return abs(span - 12 * (len(records) - 1)) < 1


def _runs() -> dict[str, list[dict]]:
    """
    The last usable run of :data:`aia_94_leak_paper.noise.frames` consecutive
    frames through the thick filter before :data:`switch`, and the first
    through the thin filter after it.
    """
    records = _jsoc.records(
        np.datetime64(switch - _window), 2 * _window.seconds, 94, ("image", "spikes")
    )
    cut = np.datetime64(switch)
    before = [r for r in records if np.datetime64(r["T_OBS"].rstrip("Z")) < cut]
    after = [r for r in records if np.datetime64(r["T_OBS"].rstrip("Z")) >= cut]
    n = noise.frames
    result = {}
    for name, candidates, order in (
        ("thick", before, range(len(before) - n, -1, -1)),
        ("thin", after, range(len(after) - n + 1)),
    ):
        for i in order:
            run = candidates[i : i + n]
            if _consecutive(run) and noise.usable(run, f"94_{name.upper()}"):
                result[name] = run
                break
        else:
            raise ValueError(f"no usable run through the {name} filter near {switch}")
    return result


def build(
    cache: pathlib.Path,
    flats: pathlib.Path,
    output: pathlib.Path,
) -> list[pathlib.Path]:
    """
    Write the archive file of each filter and return their paths.

    Parameters
    ----------
    cache
        Where the level-1 images and their despike records are kept, one
        directory per filter.
    flats
        Where the flat fields are kept.
    output
        Where to write the archive files.
    """
    output.mkdir(parents=True, exist_ok=True)
    result = []
    for name, records in _runs().items():
        hdul = noise.accumulate(records, cache / name, flats)
        hdul[0].header["WAVE_STR"] = (
            records[0]["WAVE_STR"],
            "the focal-plane filter of every frame",
        )
        path = output / f"aia94_noise_{name}.fits"
        hdul.writeto(path, overwrite=True)
        result.append(path)
    return result
