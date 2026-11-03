"""
Rebuild the archive of AIA 94 Å through each of its focal-plane filters
from JSOC.

``python -m aia_94_leak_paper.filters build`` writes one FITS file for the
run through the thin filter and one for each frame through the thick filter,
:data:`aia_94_leak_paper.filters.frames`, and
``python -m aia_94_leak_paper.filters pack`` bundles them into the tar file
published as a release asset, printing its SHA-256 digest for
:data:`aia_94_leak_paper.filters.sha256`. Building needs ``drms``
(``pip install -e .[data]``).
"""

import argparse
import hashlib
import pathlib
import tarfile

import aia_94_leak_paper


def main() -> None:
    filters = aia_94_leak_paper.filters
    parser = argparse.ArgumentParser(prog="python -m aia_94_leak_paper.filters")
    commands = parser.add_subparsers(dest="command", required=True)

    build = commands.add_parser("build", help="write the archive files from JSOC")
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
        for path in filters.build(
            cache=args.cache, flats=args.flats, output=args.output
        ):
            print(f"{path.name}, {path.stat().st_size / 2**20:.1f} MiB", flush=True)
    elif args.command == "pack":
        tar = args.output / "aia-94-filters.tar"
        with tarfile.open(tar, "w") as t:
            names = ["aia94_noise_thin"] + [f"aia94_frame_{n}" for n in filters.frames]
            for name in names:
                t.add(args.output / f"{name}.fits", arcname=f"{name}.fits")
        digest = hashlib.sha256(tar.read_bytes()).hexdigest()
        print(f"{tar}: {tar.stat().st_size / 2**20:.1f} MiB, SHA-256 {digest}")


if __name__ == "__main__":
    main()
