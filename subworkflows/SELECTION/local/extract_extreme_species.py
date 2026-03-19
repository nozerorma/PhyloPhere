#!/usr/bin/env python3
"""
Extract top/bottom extreme species lists from trait_stats.csv.

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
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = [f.strip() for f in (reader.fieldnames or [])]
        reader.fieldnames = fieldnames

        missing = {"species", "global_label"} - set(fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns in '{path}': {', '.join(sorted(missing))}"
            )

        top_species = set()
        bottom_species = set()
        seen_labels = {}

        for row in reader:
            species = (row.get("species") or "").strip()
            global_label = (row.get("global_label") or "").strip().lower()
            if not species or not global_label:
                continue
            if global_label not in {"high_extreme", "low_extreme"}:
                continue
            if species in seen_labels and seen_labels[species] != global_label:
                raise ValueError(
                    f"Species '{species}' has conflicting global_label values "
                    f"('{seen_labels[species]}' and '{global_label}') in '{path}'"
                )
            seen_labels[species] = global_label
            if global_label == "high_extreme":
                top_species.add(species)
            elif global_label == "low_extreme":
                bottom_species.add(species)

    if not top_species:
        raise ValueError(
            f"No top species found in '{path}' (expected global_label == high_extreme)"
        )
    if not bottom_species:
        raise ValueError(
            f"No bottom species found in '{path}' (expected global_label == low_extreme)"
        )

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
        description="Extract top/bottom extreme species from trait_stats.csv"
    )
    parser.add_argument("--stats-file", required=True, help="Input trait_stats.csv file")
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
