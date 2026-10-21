"""
Rebuild the archive of the six coronal channels of AIA from JSOC.

``python -m aia_94_leak_paper.dems build`` writes one FITS file per day of
:data:`aia_94_leak_paper.eclipses.dates` (or of each ``--date`` given), and
``python -m aia_94_leak_paper.dems pack`` bundles them into the tar file
published as a release asset, printing its SHA-256 digest for
:data:`aia_94_leak_paper.dems.sha256`. Building needs ``drms`` for the flat
fields (``pip install -e .[data]``).
"""

import argparse
import datetime
import hashlib
import pathlib
import tarfile

import aia_94_leak_paper


def main() -> None:
    dems = aia_94_leak_paper.dems
    dates = aia_94_leak_paper.eclipses.dates
    parser = argparse.ArgumentParser(prog="python -m aia_94_leak_paper.dems")
    commands = parser.add_subparsers(dest="command", required=True)

    build = commands.add_parser("build", help="write the archive files from JSOC")
    build.add_argument(
        "--date",
        action="append",
        type=datetime.date.fromisoformat,
        help="a day to build; every day if not given",
    )
    build.add_argument(
        "--cache",
        type=pathlib.Path,
        required=True,
        help="where the synoptic images are kept",
    )
    build.add_argument(
        "--flats",
        type=pathlib.Path,
        required=True,
        help="where the flat fields are kept",
    )
    build.add_argument(
        "--output",
        type=pathlib.Path,
        required=True,
        help="where to write the archive files",
    )

    pack = commands.add_parser(
        "pack", help="bundle the archive files into the release asset"
    )
    pack.add_argument(
        "--output", type=pathlib.Path, required=True, help="where the archive files are"
    )

    args = parser.parse_args()
    if args.command == "build":
        for date in args.date or dates:
            path = dems.build(
                date, cache=args.cache, flats=args.flats, output=args.output
            )
            print(
                f"{date}: {path.name}, {path.stat().st_size / 2**20:.1f} MiB",
                flush=True,
            )
    elif args.command == "pack":
        tar = args.output / "aia-dems.tar"
        with tarfile.open(tar, "w") as t:
            for date in dates:
                name = f"aia_dems_{date:%Y%m%d}.fits"
                t.add(args.output / name, arcname=name)
        digest = hashlib.sha256(tar.read_bytes()).hexdigest()
        print(f"{tar}: {tar.stat().st_size / 2**20:.1f} MiB, SHA-256 {digest}")


if __name__ == "__main__":
    main()
