#!/usr/bin/env python3
"""
Write a plain-text file listing the foreground species names to be used as
repeated `--branches <name>` args in HyPhy MoleRate invocations.
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Write MoleRate foreground branches from a species list"
    )
    parser.add_argument(
        "--species-file",
        required=True,
        help="Plain-text species list, one species name per line",
    )
    parser.add_argument("--output", required=True, help="Output plain-text file")
    args = parser.parse_args()

    fg_species = []
    with open(args.species_file) as fh:
        for line in fh:
            species = line.strip()
            if species and not species.startswith("#"):
                fg_species.append(species)

    if not fg_species:
        sys.exit(
            f"ERROR extract_fg_branches: no foreground species found in '{args.species_file}'"
        )

    with open(args.output, "w") as fh:
        for sp in fg_species:
            fh.write(sp + "\n")

    sys.stderr.write(
        f"INFO extract_fg_branches: wrote {len(fg_species)} foreground species to '{args.output}'\n"
    )


if __name__ == "__main__":
    main()
