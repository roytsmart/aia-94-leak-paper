"""
The leak through each of the focal-plane filters of AIA 94 Å: the thick
filter attenuates it about as much as it does the EUV.
"""

from ._archive import archive, load, sha256, url, version
from ._build import build, filters, switch
from ._model import (
    Attenuation,
    attenuation,
    band,
    euv_table,
    metal,
    metal_visible,
    tables,
    visible,
    wavelength_euv,
    within,
    zirconium,
)

__all__ = [
    "Attenuation",
    "archive",
    "attenuation",
    "band",
    "build",
    "euv_table",
    "filters",
    "load",
    "metal",
    "metal_visible",
    "sha256",
    "switch",
    "tables",
    "url",
    "version",
    "visible",
    "wavelength_euv",
    "within",
    "zirconium",
]
