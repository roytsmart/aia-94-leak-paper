"""
Rebuild the archive of the images AIA takes without a focal-plane filter
from JSOC.

``python -m aia_94_leak_paper.entrance build`` writes the FITS file of every
month up to ``--last``, and ``python -m aia_94_leak_paper.entrance pack``
bundles it into the tar file published as a release asset, printing its
SHA-256 digest for :data:`aia_94_leak_paper.entrance.sha256`. Building needs
``drms`` (``pip install -e .[data]``).
"""

import argparse
import datetime
import hashlib
import pathlib
import tarfile

import aia_94_leak_paper


def main() -> None:
    entrance = aia_94_leak_paper.entrance
    parser = argparse.ArgumentParser(prog="python -m aia_94_leak_paper.entrance")
    commands = parser.add_subparsers(dest="command", required=True)

    build = commands.add_parser("build", help="write the archive file from JSOC")
    build.add_argument(
        "--last",
        type=datetime.date.fromisoformat,
        required=True,
        help="a day of the last month to include",
    )
    build.add_argument(
        "--cache",
        type=pathlib.Path,
        required=True,
        help="where the level-1 images are kept",
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
        help="where to write the archive file",
    )

    pack = commands.add_parser(
        "pack", help="bundle the archive file into the release asset"
    )
    pack.add_argument(
        "--output", type=pathlib.Path, required=True, help="where the archive file is"
    )

    args = parser.parse_args()
    if args.command == "build":
        path = entrance.build(
            args.last, cache=args.cache, flats=args.flats, output=args.output
        )
        print(f"{path.name}, {path.stat().st_size / 2**20:.1f} MiB", flush=True)
    elif args.command == "pack":
        tar = args.output / "aia-entrance.tar"
        with tarfile.open(tar, "w") as t:
            t.add(args.output / "aia_entrance.fits", arcname="aia_entrance.fits")
        digest = hashlib.sha256(tar.read_bytes()).hexdigest()
        print(f"{tar}: {tar.stat().st_size / 2**20:.1f} MiB, SHA-256 {digest}")


if __name__ == "__main__":
    main()
