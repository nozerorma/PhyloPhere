#!/usr/bin/env python3
"""
tree_cleanup.py  —  Curate a species tree to match alignment species names.

Uses NCBI taxonomy IDs as a shared key to translate tree tip labels that
differ from alignment header names (e.g. taxonomic synonyms / genus renames).
Tips with no match in the alignment species set are pruned.

Expected taxid file format (TSV or CSV):
    tax_id  species  [other columns ignored]
  where both tree tip names AND alignment species names appear under 'species',
  each paired with their stable NCBI tax_id.

Inputs
------
  --tree          Newick species tree
  --ali-sp-names  Flat file: one alignment species name per line (\\n-separated)
  --tax-id        TSV/CSV with columns 'tax_id' and 'species'
  --output        Output newick (curated tree)
  --report        Output TSV: original_name, curated_name, fate (kept/renamed/pruned)
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import dendropy


def normalize_label(label: str) -> str:
    label = (label or "").strip()
    label = label.strip("'\"")
    label = re.sub(r"\s+", "_", label)
    return label


def load_ali_sp_names(path: str) -> set:
    names = set()
    with open(path) as fh:
        for line in fh:
            name = normalize_label(line)
            if name:
                names.add(name)
    return names


def load_tax_id_map(path: str):
    """Return (name_to_taxid, taxid_to_names) from a TSV/CSV taxid file."""
    sep = "\t" if path.endswith(".tsv") else ","
    name_to_taxid: dict[str, str] = {}
    taxid_to_names: dict[str, list[str]] = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        for row in reader:
            row = {k: v.strip() for k, v in row.items() if k}
            taxid   = row.get("tax_id", "").strip()
            species = normalize_label(row.get("species", ""))
            if not taxid or not species:
                continue
            name_to_taxid[species] = taxid
            taxid_to_names.setdefault(taxid, []).append(species)

    return name_to_taxid, taxid_to_names


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree",         required=True, help="Input newick tree")
    parser.add_argument("--ali-sp-names", required=True, help="Alignment species names (flat list)")
    parser.add_argument("--tax-id",       required=True, help="Taxid TSV/CSV mapping file")
    parser.add_argument("--output",       required=True, help="Output curated newick tree")
    parser.add_argument("--report",       required=True, help="Output TSV report")
    args = parser.parse_args()

    ali_sp = load_ali_sp_names(args.ali_sp_names)
    name_to_taxid, taxid_to_names = load_tax_id_map(args.tax_id)

    # Build taxid → alignment species name (for ali species that have a taxid entry)
    taxid_to_ali: dict[str, str] = {}
    for sp in ali_sp:
        taxid = name_to_taxid.get(sp)
        if taxid and taxid not in taxid_to_ali:
            taxid_to_ali[taxid] = sp

    tree = dendropy.Tree.get(path=args.tree, schema="newick",
                             preserve_underscores=True)
    sample_labels = [taxon.label for taxon in list(tree.taxon_namespace)[:5]]
    print(f"[tree_cleanup] Input tips={len(tree.taxon_namespace)} sample={sample_labels}", flush=True)

    kept = renamed = 0
    taxa_to_prune = []
    report_rows: list[tuple[str, str, str]] = []

    for taxon in list(tree.taxon_namespace):
        tip = normalize_label(taxon.label)

        if tip in ali_sp:
            report_rows.append((tip, tip, "kept"))
            kept += 1
        else:
            taxid   = name_to_taxid.get(tip)
            ali_name = taxid_to_ali.get(taxid) if taxid else None
            if ali_name:
                report_rows.append((tip, ali_name, "renamed"))
                taxon.label = normalize_label(ali_name)
                renamed += 1
            else:
                report_rows.append((tip, "", "pruned"))
                taxa_to_prune.append(taxon)

    tree.prune_taxa(taxa_to_prune)
    tree.purge_taxon_namespace()

    tree.write(
        path=args.output,
        schema="newick",
        preserve_spaces=False,
        unquoted_underscores=True,
    )
    print(f"[tree_cleanup] Wrote curated tree to {args.output}", flush=True)

    with open(args.report, "w", newline="") as fh:
        fh.write("original_name\tcurated_name\tfate\n")
        for row in report_rows:
            fh.write("\t".join(row) + "\n")

    pruned = len(taxa_to_prune)
    print(f"[tree_cleanup] kept={kept}  renamed={renamed}  pruned={pruned}", flush=True)
    print(f"[tree_cleanup] Output tree: {args.output}  ({kept + renamed} tips retained)", flush=True)

    if kept + renamed == 0:
        print("[tree_cleanup] ERROR: no tree tips matched alignment species — "
              "check that tax_id file covers both tree and alignment naming conventions.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
