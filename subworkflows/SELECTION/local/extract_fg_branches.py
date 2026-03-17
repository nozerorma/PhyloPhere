#!/usr/bin/env python3
"""
extract_fg_branches.py
──────────────────────
Read a caastools traitfile and write a plain-text file listing the foreground
species names (one per line), to be used as repeated `--branches <name>` args
in HyPhy MoleRate invocations.

Foreground = contrast_group 1 for TOP direction, contrast_group 0 for BOTTOM.

Usage
-----
    python extract_fg_branches.py \\
        --traitfile traitfile.tab \\
        --output    fg_branches.txt \\
        [--flip]

Output
------
One species name per line, matching the taxon labels in the alignment / tree.
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Extract foreground branch names from a caastools traitfile"
    )
    parser.add_argument("--traitfile", required=True, help="caastools traitfile (3-col TSV)")
    parser.add_argument("--output",    required=True, help="Output plain-text file")
    parser.add_argument("--flip",      action="store_true",
                        help="Swap FG/BG: use contrast_group==0 as FG (BOTTOM direction)")
    args = parser.parse_args()

    fg_species = []
    with open(args.traitfile) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            species = parts[0]
            try:
                group = int(parts[1])
            except ValueError:
                continue
            is_fg = (group == 0) if args.flip else (group == 1)
            if is_fg:
                fg_species.append(species)

    if not fg_species:
        sys.exit(
            f"ERROR extract_fg_branches: no foreground species found in '{args.traitfile}' "
            f"(flip={args.flip})"
        )

    with open(args.output, "w") as fh:
        for sp in fg_species:
            fh.write(sp + "\n")

    sys.stderr.write(
        f"INFO extract_fg_branches: wrote {len(fg_species)} foreground species to '{args.output}'\n"
    )


if __name__ == "__main__":
    main()
