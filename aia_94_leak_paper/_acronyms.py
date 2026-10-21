import aastex

__all__ = [
    "acronyms",
]


def acronyms() -> list[aastex.Acronym]:
    """
    A list of acronyms that are used in the body of the article.
    """
    return [
        aastex.Acronym("AIA", "the Atmospheric Imaging Assembly"),
        aastex.Acronym("SDO", "the Solar Dynamics Observatory"),
        aastex.Acronym("SUVI", "the Solar Ultraviolet Imager"),
        aastex.Acronym(
            "GOES", "the Geostationary Operational Environmental Satellite", short=True
        ),
        aastex.Acronym("EVE", "the EUV Variability Experiment"),
        aastex.Acronym("XRT", "the X-Ray Telescope"),
        aastex.Acronym("MaGIXS", "the Marshall Grazing Incidence X-ray Spectrometer"),
        aastex.Acronym("EUV", "extreme ultraviolet"),
        aastex.Acronym("DEM", "differential emission measure", plural=True),
        aastex.Acronym("CCD", "charge-coupled device", plural=True),
        aastex.Acronym("PSF", "point spread function"),
        aastex.Acronym("FOV", "field of view"),
    ]
