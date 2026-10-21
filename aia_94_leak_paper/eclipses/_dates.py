import datetime

__all__ = [
    "dates",
]

dates: tuple[datetime.date, ...] = tuple(
    datetime.date.fromisoformat(date)
    for date in (
        # the first eclipse season after AIA began to observe
        "2010-09-26",
        # the first season of each year, while the leak was faint
        "2011-03-19",
        "2012-03-15",
        "2013-03-11",
        "2014-03-08",
        "2015-03-04",
        # both seasons of each year since the leak began to grow
        "2016-02-29",
        "2016-08-29",
        "2017-02-24",
        "2017-08-25",
        "2018-02-20",
        "2018-08-21",
        "2019-02-17",
        "2019-08-18",
        "2020-02-13",
        "2020-08-13",
        "2021-02-08",
        "2021-08-09",
        "2022-02-05",
        "2022-08-05",
        "2023-02-01",
        "2023-08-01",
        "2024-01-28",
        "2024-07-27",
        "2025-01-23",
        "2025-07-23",
        "2026-01-19",
        "2026-07-15",
    )
)
"""
The days whose eclipse ingress this article measures the leak in.

SDO is eclipsed by the Earth once a day for three to four weeks twice a
year, and the seasons drift about four days earlier each year. Each day
here is near the middle of its season, when the eclipse is longest and the
edge of the atmosphere crosses the solar disk quickly and squarely.
"""
