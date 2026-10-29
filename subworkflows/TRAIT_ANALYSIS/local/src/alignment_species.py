#!/usr/bin/env python3
# alignment_species.py — Collect species record IDs using the pipeline alignment format.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/

"""Read the same top-level alignment files as CT, including file symlinks.

FASTA is streamed through SeqIO; other formats use AlignIO like CT. Descriptions
are excluded from species IDs. The output can be reused with --ali_sp_names.
"""

import argparse
from pathlib import Path

from Bio import AlignIO, SeqIO

from tree_cleanup import normalize_label


def alignment_species(directory, fmt):
    """Return sorted species IDs, rejecting unreadable files and name collisions."""
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"alignment directory not found: {directory}")
    names = {}
    count = 0
    for path in sorted(directory.iterdir()):
        # Match CT's top-level listing and sidecar exclusions.
        if path.suffix in {".txt", ".tsv", ".csv", ".log", ".map"} or path.is_dir():
            continue
        if not path.is_file():
            raise ValueError(f"alignment file is unreadable or not a regular file: {path}")
        count += 1
        try:
            with path.open() as fh:
                records = SeqIO.parse(fh, "fasta") if fmt == "fasta" else AlignIO.read(fh, fmt)
                found = False
                for record in records:
                    found = True
                    name = normalize_label(record.id)
                    if not name:
                        raise ValueError("empty species record ID")
                    if name in names and names[name] != record.id:
                        raise ValueError(f"record IDs {names[name]!r} and {record.id!r} "
                                         f"normalize to {name!r}")
                    names[name] = record.id
                if not found:
                    raise ValueError(f"no species records found for format {fmt!r}")
        except (ValueError, OSError) as exc:
            raise ValueError(f"{path}: {exc}") from exc
    if not count:
        raise ValueError(f"{directory}: no alignment files found")
    # Preserve the identifiers for downstream normalization and collision checking.
    return [names[name] for name in sorted(names)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--format", default="fasta")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        names = alignment_species(args.alignment, args.format)
        Path(args.output).write_text("\n".join(names) + "\n")
    except (ValueError, OSError) as exc:
        parser.exit(1, f"[NAME_CURATION] ERROR: {exc}\n")
    print(f"[NAME_CURATION] Found {len(names)} unique species record IDs")


if __name__ == "__main__":
    main()
