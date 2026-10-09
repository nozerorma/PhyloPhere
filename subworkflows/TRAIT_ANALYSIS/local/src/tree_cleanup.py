#!/usr/bin/env python3
# tree_cleanup.py — Curate a species tree so its tip labels match the alignment species names.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/

"""
TreeCleanup: Renames or prunes the tips of a species tree so that every
remaining tip carries the species name used in the alignment FASTA headers.

NCBI tax_ids are the shared key: a tip whose name is not among the alignment
species is renamed to the alignment species that has the same tax_id (taxonomic
synonyms, genus renames), and pruned when there is none. Without a tax_id map,
tips are matched by exact name only. Labels are normalized (quotes stripped,
whitespace to underscores) before comparison. The script exits with an error
when no tip survives.

Called by:  TREE_CLEANUP Nextflow process (ta_name_curation.nf → tree_cleanup.py)
Inputs:     --tree          Newick species tree
            --ali-sp-names  Text file with one alignment species name per line. The sentinel
                            NO_FILE means no alignment names: every tree tip is then a
                            canonical species and only the trait table is curated
            --tax-id        TSV/CSV with columns [tax_id, species] (other columns
                            ignored); tree-side and alignment-side names both
                            appear under `species`. The sentinel NO_FILE, or a
                            missing file, means no map
            --output        Path for the curated Newick tree
            --report        Path for the per-tip report TSV
            --traits        Optional trait table (CSV/TSV). Its species are matched to
                            the curated tree tips: by exact name, then through tax_id
                            (the table's own `tax_id` column when present, else the map)
                            Without it only the tree is curated
            --sp-col        Species column of --traits (default: species)
            --traits-out    Path for the curated trait table (needed with --traits; species
                            renamed to the tip names, unmatched rows removed)
            --species-table Path for the table of canonical species and their tax_ids
            --species-report Path for the plain-text curation report
            --taxid-map     Path for the curated tax_id map: one row per canonical species with
                            its resolved tax_id and its family, in the 5-column layout of the
                            taxonomy file [tax_id, species, family, rank, name_class]
Outputs:    Curated Newick tree; TSV [original_name, curated_name, fate] with
            fate in {kept, renamed, pruned} (curated_name is empty for pruned).
            With --species-table, --taxid-map and --species-report: the species table
            [species, tax_id, tax_id_resolved, note], the curated tax_id map and the report
            [source, original_name, curated_name, status, reason], status in
            {maintained, changed, removed}, covering tree tips, trait species and
            synthetic tax_ids

Canonical species are the tips of the curated tree. When several of them share a
tax_id, the alphabetically first keeps it and each other one receives
tax_id + i (probing past every tax_id of the map). Removals are reported, never fatal.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import dendropy


def normalize_label(label: str) -> str:
    """Strip whitespace and quotes and turn inner whitespace into underscores."""
    label = (label or "").strip()
    label = label.strip("'\"")
    label = re.sub(r"\s+", "_", label)
    return label


def load_ali_sp_names(path: str) -> set:
    """Return the set of normalized species names listed in `path`, one per line."""
    names = set()
    with open(path) as fh:
        for line in fh:
            name = normalize_label(line)
            if name:
                names.add(name)
    return names


def load_tax_id_map(path: str):
    """Return (name_to_taxid, taxid_to_names, name_to_family) from a TSV/CSV taxid file.

    "NO_FILE" is the Nextflow sentinel for an absent optional file (the workflows
    pass it when params.tax_id is unset). The sentinel, an empty path and a
    missing file all give empty maps, so tips are then matched by exact name only,
    with no translation between naming conventions.
    """
    if not path or path == "NO_FILE" or not Path(path).exists():
        return {}, {}, {}

    sep = "\t" if path.endswith(".tsv") else ","
    name_to_taxid: dict[str, str] = {}
    taxid_to_names: dict[str, list[str]] = {}
    name_to_family: dict[str, str] = {}

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
            name_to_family[species] = row.get("family", "")

    return name_to_taxid, taxid_to_names, name_to_family


def read_table(path: str):
    """Return (header, rows, delimiter) of a CSV/TSV trait table; rows are lists of strings."""
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader)
        rows = [r for r in reader if r]
    return header, rows, sep


def resolve_tip_taxids(tips, name_to_taxid, reserved):
    """Resolve the tax_id of every canonical species.

    Returns {species: (real_tax_id, resolved_tax_id, note)}. Species that share a real
    tax_id compete alphabetically: the first keeps it, the i-th other one gets
    real + i, probing forward past every id in `reserved`. A species without a
    tax_id in the map keeps an empty id.
    """
    by_id: dict[str, list[str]] = {}
    for sp in tips:
        tid = name_to_taxid.get(sp, "")
        by_id.setdefault(tid, []).append(sp)

    out: dict[str, tuple[str, str, str]] = {}
    used = set(reserved)
    for tid in sorted(by_id, key=lambda t: (not t.isdigit(), int(t) if t.isdigit() else 0, t)):
        group = sorted(by_id[tid])
        if not tid:
            for sp in group:
                out[sp] = ("", "", "no tax_id in the map")
            continue
        out[group[0]] = (tid, tid, "")
        for i, sp in enumerate(group[1:], start=1):
            new = int(tid) + i
            while str(new) in used:
                new += 1
            used.add(str(new))
            out[sp] = (tid, str(new),
                       f"synthetic tax_id: real tax_id {tid} is shared with {group[0]}")
    return out


def curate_traits(header, rows, sp_col, tips, name_to_taxid, resolved):
    """Match trait rows to the curated tree tips.

    Returns (kept_rows, decisions). Each decision is
    (original_name, curated_name, status, reason). A row is maintained when its name is a
    tip, changed when its tax_id points to a tip, and removed otherwise or when another
    row already stands for the same tip (the row whose name is the tip wins, else the
    alphabetically first original name).
    """
    sp_i = header.index(sp_col)
    tid_i = header.index("tax_id") if "tax_id" in header else None
    tip_by_id: dict[str, list[str]] = {}
    for sp in sorted(tips):
        tid = name_to_taxid.get(sp, "")
        if tid:
            tip_by_id.setdefault(tid, []).append(sp)

    cand = []  # (row index, original name, target tip or None, status, reason)
    for k, row in enumerate(rows):
        name = normalize_label(row[sp_i])
        if name in tips:
            cand.append((k, name, name, "maintained", "name matches a tree tip"))
            continue
        tid = (row[tid_i].strip() if tid_i is not None else "") or name_to_taxid.get(name, "")
        if not tid:
            cand.append((k, name, None, "removed",
                         "no overlap with the tree: the name is not a tip and has no tax_id"))
        elif tid in tip_by_id:
            target = tip_by_id[tid][0]
            cand.append((k, name, target, "changed",
                         f"renamed to {target} through the shared tax_id {tid}"))
        else:
            cand.append((k, name, None, "removed",
                         f"no overlap with the tree: no tip has the tax_id {tid}"))

    # One row per tip: the row named like the tip wins, then the alphabetically first name.
    winner: dict[str, int] = {}
    for k, name, target, status, _ in sorted(cand, key=lambda c: (c[2] or "", c[1] != c[2], c[1], c[0])):
        if target is not None and target not in winner:
            winner[target] = k

    kept, decisions = [], []
    for k, name, target, status, reason in cand:
        if target is not None and winner[target] != k:
            decisions.append((name, "", "removed",
                              f"duplicate of the row kept for {target} (same canonical species)"))
            continue
        decisions.append((name, target or "", status, reason))
        if target is not None:
            row = list(rows[k])
            row[sp_i] = target
            if tid_i is not None and resolved.get(target, ("", "", ""))[1]:
                row[tid_i] = resolved[target][1]
            kept.append(row)
    return kept, decisions


def write_species_outputs(args, header, kept, decisions, tree_decisions, resolved, sep,
                          name_to_family=None):
    """Write the curated trait table (when there is one), the species table, the curated
    tax_id map and the text report."""
    name_to_family = name_to_family or {}
    if header is not None:
        with open(args.traits_out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter=sep)
            w.writerow(header)
            w.writerows(kept)

    with open(args.species_table, "w", newline="") as fh:
        fh.write("species\ttax_id\ttax_id_resolved\tnote\n")
        for sp in sorted(resolved):
            fh.write("\t".join((sp, *resolved[sp])) + "\n")

    # Same layout as the --tax-id taxonomy file, so every consumer of a tax_id map (the name to
    # tax_id readers and the clade variability, which reads the family) takes it as is: one row
    # per canonical species, with the resolved (unique) tax_id.
    with open(args.taxid_map, "w", newline="") as fh:
        fh.write("tax_id\tspecies\tfamily\trank\tname_class\n")
        for sp in sorted(resolved):
            if resolved[sp][1]:
                fh.write(f"{resolved[sp][1]}\t{sp}\t{name_to_family.get(sp, '')}\tspecies\tscientific name\n")

    synth = [(sp, v) for sp, v in sorted(resolved.items()) if v[0] != v[1] and v[1]]
    lines = [("tree", *d) for d in tree_decisions] + [("trait", *d) for d in decisions]
    lines += [("tax_id", sp, sp, "changed", v[2]) for sp, v in synth]
    with open(args.species_report, "w") as fh:
        n_removed = sum(1 for x in lines if x[0] != "tax_id" and x[3] == "removed")
        fh.write("# Species curation report\n")
        fh.write(f"# canonical species (curated tree tips): {len(resolved)}\n")
        fh.write(f"# trait rows: {len(decisions)}, kept: "
                 f"{sum(1 for d in decisions if d[2] != 'removed')}, "
                 f"removed: {sum(1 for d in decisions if d[2] == 'removed')}\n")
        fh.write(f"# tree tips removed: {sum(1 for d in tree_decisions if d[2] == 'removed')}\n")
        fh.write(f"# synthetic tax_ids: {len(synth)}\n")
        fh.write("source\toriginal_name\tcurated_name\tstatus\treason\n")
        order = {"removed": 0, "changed": 1, "maintained": 2}
        for x in sorted(lines, key=lambda x: (x[0], order[x[3]], x[1])):
            fh.write("\t".join(x) + "\n")
    return len(synth)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree",         required=True, help="Input newick tree")
    parser.add_argument("--ali-sp-names", required=True, help="Alignment species names (flat list)")
    parser.add_argument("--tax-id",       required=True, help="Taxid TSV/CSV mapping file")
    parser.add_argument("--output",       required=True, help="Output curated newick tree")
    parser.add_argument("--report",       required=True, help="Output TSV report")
    parser.add_argument("--traits",         default=None, help="Trait table (CSV/TSV) to curate")
    parser.add_argument("--sp-col",         default="species", help="Species column of --traits")
    parser.add_argument("--traits-out",     default=None, help="Output curated trait table")
    parser.add_argument("--species-table",  default=None, help="Output table of canonical species and tax_ids")
    parser.add_argument("--species-report", default=None, help="Output text report of the curation")
    parser.add_argument("--taxid-map",      default=None, help="Output curated tax_id map [tax_id, species]")
    args = parser.parse_args()
    if args.traits and not args.traits_out:
        parser.error("--traits needs --traits-out")
    if bool(args.species_table) != bool(args.species_report) or bool(args.species_table) != bool(args.taxid_map):
        parser.error("--species-table, --species-report and --taxid-map go together")

    name_to_taxid, taxid_to_names, name_to_family = load_tax_id_map(args.tax_id)

    tree = dendropy.Tree.get(path=args.tree, schema="newick",
                             preserve_underscores=True)
    if args.ali_sp_names == "NO_FILE":
        ali_sp = {normalize_label(t.label) for t in tree.taxon_namespace}
    else:
        ali_sp = load_ali_sp_names(args.ali_sp_names)

    # tax_id → alignment species name (the alphabetically first species wins when several share a tax_id).
    taxid_to_ali: dict[str, str] = {}
    for sp in sorted(ali_sp):
        taxid = name_to_taxid.get(sp)
        if taxid and taxid not in taxid_to_ali:
            taxid_to_ali[taxid] = sp

    sample_labels = [taxon.label for taxon in list(tree.taxon_namespace)[:5]]
    print(f"[tree_cleanup] Input tips={len(tree.taxon_namespace)} sample={sample_labels}", flush=True)

    taxa = list(tree.taxon_namespace)
    labels = {id(t): normalize_label(t.label) for t in taxa}
    fate: dict[int, tuple[str, str, str]] = {}  # taxon -> (curated name, fate, reason)
    used: set[str] = set()

    # Tips named like an alignment species stay. The others are translated through their
    # tax_id in name order, and a tip whose alignment species already has a tip is pruned
    # (a curated label stays unique).
    for t in taxa:
        if labels[id(t)] in ali_sp:
            fate[id(t)] = (labels[id(t)], "kept", "name matches an alignment species")
            used.add(labels[id(t)])
    for t in sorted((t for t in taxa if id(t) not in fate), key=lambda t: labels[id(t)]):
        tip = labels[id(t)]
        taxid = name_to_taxid.get(tip)
        ali_name = taxid_to_ali.get(taxid) if taxid else None
        if not ali_name:
            reason = (f"no alignment species with the tax_id {taxid}" if taxid
                      else "the name is not an alignment species and has no tax_id")
            fate[id(t)] = ("", "pruned", reason)
        elif ali_name in used:
            fate[id(t)] = ("", "pruned",
                           f"duplicate: the alignment species {ali_name} (tax_id {taxid}) already has a tip")
        else:
            fate[id(t)] = (ali_name, "renamed",
                           f"renamed to {ali_name} through the shared tax_id {taxid}")
            used.add(ali_name)

    taxa_to_prune = []
    report_rows: list[tuple[str, str, str]] = []
    tree_decisions: list[tuple[str, str, str, str]] = []
    status = {"kept": "maintained", "renamed": "changed", "pruned": "removed"}
    for t in taxa:
        curated, kind, reason = fate[id(t)]
        report_rows.append((labels[id(t)], curated, kind))
        tree_decisions.append((labels[id(t)], curated, status[kind], reason))
        if kind == "pruned":
            taxa_to_prune.append(t)
        else:
            t.label = curated
    kept = sum(1 for f in fate.values() if f[1] == "kept")
    renamed = sum(1 for f in fate.values() if f[1] == "renamed")

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

    if args.species_table and kept + renamed > 0:
        tips = {f[0] for f in fate.values() if f[1] != "pruned"}
        resolved = resolve_tip_taxids(sorted(tips), name_to_taxid, set(taxid_to_names))
        header, rows, sep, kept_rows, decisions = None, [], ",", [], []
        if args.traits:
            header, rows, sep = read_table(args.traits)
            if args.sp_col not in header:
                print(f"[tree_cleanup] ERROR: column '{args.sp_col}' not found in {args.traits}",
                      file=sys.stderr)
                sys.exit(1)
            kept_rows, decisions = curate_traits(header, rows, args.sp_col, tips, name_to_taxid, resolved)
        n_synth = write_species_outputs(args, header, kept_rows, decisions, tree_decisions, resolved, sep,
                                        name_to_family)
        n_rem = sum(1 for d in decisions if d[2] == "removed")
        n_chg = sum(1 for d in decisions if d[2] == "changed")
        print(f"[tree_cleanup] trait species: {len(rows)} -> {len(kept_rows)} "
              f"(changed={n_chg}, removed={n_rem}); synthetic tax_ids={n_synth}", flush=True)
        if n_rem:
            print(f"[tree_cleanup] WARNING: {n_rem} trait species removed; reasons in "
                  f"{args.species_report}", file=sys.stderr, flush=True)

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
