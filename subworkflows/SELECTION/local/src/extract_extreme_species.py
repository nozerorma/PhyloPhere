#!/usr/bin/env python3
"""
Extract disjoint top/bottom pools from candidate species TSV or legacy trait_stats CSV.

Expected input columns:
  - species
  - global_label

global_label mapping:
  - high_extreme -> top
  - low_extreme  -> bottom
"""

import argparse
import csv
import sys


def load_species_sets(path: str):
    """Read a disjoint global pool, not a paired CT design; pair is ignored."""
    top_species, bottom_species, seen = set(), set(), {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        lines = [line for line in fh if line.strip()]
    if not lines:
        raise ValueError(f"Empty species input: {path}")
    is_tsv = "\t" in lines[0]
    if is_tsv:
        rows = list(csv.reader(lines, delimiter="\t"))
        if (len(rows[0]) in {2, 3} and rows[0][0] in {"species", "taxa", "tip"}
                and rows[0][1] in {"contrast_group", "group"}):
            rows = rows[1:]
        parsed = []
        for row in rows:
            if len(row) not in {2, 3} or row[1] not in {"0", "1"}:
                raise ValueError(f"Invalid species/group row in {path}: {row!r}")
            parsed.append((row[0], row[1]))
    else:
        reader = csv.DictReader(lines)
        fields = [f.strip() for f in (reader.fieldnames or [])]
        if len(fields) != len(set(fields)) or not {"species", "global_label"} <= set(fields):
            raise ValueError(f"Expected species/global_label CSV or species/contrast_group TSV: {path}")
        reader.fieldnames = fields
        parsed = []
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Invalid CSV row in {path}: {row!r}")
            label = row["global_label"].strip().lower()
            # Non-extreme rows are valid in the legacy trait_stats CSV.
            if label not in {"high_extreme", "low_extreme", "intermediate", "normal", "non_extreme", "", "na"}:
                raise ValueError(f"Invalid global_label in {path}: {label!r}")
            parsed.append((row["species"], label))
    for species, group in parsed:
        if (not species or species != species.strip() or "," in species
                or any(c.isspace() for c in species)):
            raise ValueError(f"Invalid species identifier in {path}: {species!r}")
        if species in seen:
            raise ValueError(f"Duplicate or conflicting species in {path}: {species!r}")
        seen[species] = group
        if group in {"1", "high_extreme"}:
            top_species.add(species)
        elif group in {"0", "low_extreme"}:
            bottom_species.add(species)
    if not top_species or not bottom_species:
        raise ValueError(f"Both top and bottom species are required: {path}")
    return top_species, bottom_species


def write_species(path: str, species: set, label: str) -> None:
    with open(path, "w") as fh:
        for sp in sorted(species):
            fh.write(sp + "\n")
    sys.stderr.write(
        f"INFO extract_extreme_species: wrote {len(species)} {label} species to '{path}'\n"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Extract disjoint top/bottom species pools"
    )
    parser.add_argument("--stats-file", required=True, help="Candidate species TSV (pair ignored) or legacy trait_stats CSV")
    parser.add_argument("--out-top", required=True, help="Output top species list")
    parser.add_argument(
        "--out-bottom", required=True, help="Output bottom species list"
    )
    args = parser.parse_args()

    top_species, bottom_species = load_species_sets(args.stats_file)
    write_species(args.out_top, top_species, "top")
    write_species(args.out_bottom, bottom_species, "bottom")


if __name__ == "__main__":
    main()
