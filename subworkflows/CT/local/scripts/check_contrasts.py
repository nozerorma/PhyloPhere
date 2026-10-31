#!/usr/bin/env python3
"""Validate and gate the generated contrast design; retain only active metadata."""
import csv
import json
import math
from pathlib import Path
import re
import shutil
import sys


def positive_integer(value, name):
    if not re.fullmatch(r"[0-9]+", str(value)) or not 1 <= int(value) <= 2147483647:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def read_design(path):
    pairs, species = {}, set()
    with Path(path).open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            if len(row) != 3 or row[1] not in {"0", "1"} or not row[0]:
                raise ValueError(f"Invalid species/group/pair row in {path}: {row}")
            pair = positive_integer(row[2], "pair")
            if row[0] in species or row[1] in pairs.setdefault(pair, {}):
                raise ValueError(f"Duplicate species or pair side in {path}: {row}")
            species.add(row[0])
            pairs[pair][row[1]] = row[0]
    if any(set(sides) != {"0", "1"} for sides in pairs.values()):
        raise ValueError(f"Unmatched foreground/background pair in {path}")
    if set(pairs) != set(range(1, len(pairs) + 1)):
        raise ValueError(f"Pair IDs must be contiguous from 1 in {path}")
    return pairs


def filter_metadata(source, destination, designs, *, pairs):
    if not source.exists():
        return  # Legacy unweighted designs remain supported.
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"hypothesis_id"}
        if pairs:
            required |= {"pair", "species1", "species2", "pss_score"}
        if not required <= set(fields):
            raise ValueError(f"Missing metadata columns in {source}")
        rows = [row for row in reader if row["hypothesis_id"] in designs]
    seen = set()
    for row in rows:
        hyp = row["hypothesis_id"]
        if pairs:
            pair = positive_integer(row["pair"], "metadata pair")
            key = (hyp, pair)
            design = designs[hyp].get(pair)
            if design is None or (row["species1"], row["species2"]) != (design["1"], design["0"]):
                raise ValueError(f"Metadata species/pair mismatch: {key}")
            weight = float(row["pss_score"])
            if not math.isfinite(weight) or weight < 0:
                raise ValueError(f"Invalid PSS weight: {key}")
        else:
            key = hyp
            if "num_pairs" in row and positive_integer(row["num_pairs"], "num_pairs") != len(designs[hyp]):
                raise ValueError(f"Metadata pair count mismatch: {hyp}")
        if key in seen:
            raise ValueError(f"Duplicate metadata entry: {key}")
        seen.add(key)
    expected = {(h, p) for h, design in designs.items() for p in design} if pairs else set(designs)
    if seen != expected:
        raise ValueError(f"Metadata does not cover the active design: {source}")
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def gate(arguments):
    minimum = positive_integer(arguments["min_contrasts"], "min_contrasts")
    canonical = read_design(arguments["traitfile"])
    count = len(canonical)
    skipped = count < minimum
    if skipped:
        Path("low_contrasts.skip").write_text(
            f"trait={arguments['traitname']}\tn_contrasts={count}\tmin_required={minimum}\n"
        )
    else:
        source = Path(arguments["trait_dir"])
        files = sorted(source.glob("traitfile_H*.tab")) if arguments["multi_hypothesis"] else []
        if not files:
            files = [source / "traitfile_H1.tab"] if (source / "traitfile_H1.tab").is_file() else []
        designs = {f.stem.removeprefix("traitfile_"): read_design(f) for f in files}
        if not files:
            designs = {"H1": canonical}
        if designs.get("H1") != canonical or any(len(d) != count for d in designs.values()):
            raise ValueError("Active hypotheses must include canonical H1 and have its pair count")
        if any(not re.fullmatch(r"H[1-9][0-9]*", hyp) for hyp in designs):
            raise ValueError("Invalid hypothesis filename")
        output = Path("traitfiles_ok_dir")
        output.mkdir()
        for f in files:
            shutil.copyfile(f, output / f.name)
        if not files:
            shutil.copyfile(arguments["traitfile"], output / "traitfile_H1.tab")
        if (source / "observed_selection.rds").is_file():
            shutil.copyfile(source / "observed_selection.rds", output / "observed_selection.rds")
        for name, pair_table in [("contrast_hypotheses_pairs.tsv", True),
                                 ("contrast_hypotheses_summary.tsv", False)]:
            filter_metadata(source / name, output / name, designs, pairs=pair_table)
        shutil.copyfile(arguments["traitfile"], "traitfile_ok.tab")
        shutil.copyfile(arguments["permulation_traitfile"], "permulation_traitfile_ok.tab")
    Path("gate_status.json").write_text(json.dumps({
        "skipped": skipped, "trait": arguments["traitname"],
        "n_contrasts": count, "min_required": minimum,
    }) + "\n")
    print(f"CHECK_MIN_CONTRASTS: {count} pairs; minimum {minimum}; {'skipped' if skipped else 'passing'}")


if __name__ == "__main__":
    gate(json.loads(Path(sys.argv[1]).read_text()))
