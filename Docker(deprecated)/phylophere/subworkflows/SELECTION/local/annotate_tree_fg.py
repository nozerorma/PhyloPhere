#!/usr/bin/env python3
"""
annotate_tree_fg.py
───────────────────
Read a caastools traitfile and a Newick tree, then write an annotated Newick
where every foreground leaf node carries the `{Foreground}` label consumed by
HyPhy FADE's `--branches Foreground` option.

Foreground is defined by contrast_group == 1 in the traitfile.

Usage
-----
    python annotate_tree_fg.py \\
        --traitfile  traitfile.tab \\
        --tree       tree.nwk \\
        --output     annotated.nwk \\
        [--flip]           # swap FG/BG (use for BOTTOM direction tests)
        [--label Foreground]   # HyPhy branch label (default: Foreground)

Input traitfile format (3-col TSV, no header)
─────────────────────────────────────────────
    species  contrast_group  pair
    Homo_sapiens  1  1
    Pan_troglodytes  0  1
    ...
"""

import argparse
import sys
import re
import dendropy


def parse_traitfile(path: str, flip: bool) -> set:
    """Return set of foreground species names from the traitfile."""
    fg_species = set()
    with open(path) as fh:
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
            # contrast_group == 1 → foreground (high-trait extreme)
            # --flip inverts: group == 0 becomes FG (low-trait extreme for BOTTOM tests)
            is_fg = (group == 0) if flip else (group == 1)
            if is_fg:
                fg_species.add(species)
    return fg_species


def annotate_newick(newick_str: str, fg_species: set, label: str) -> str:
    """
    Annotate foreground leaf branch names in a Newick string with {label}.

    Strategy: replace bare leaf labels that match FG species names with
    `name{label}`. We use a regex that matches word-boundary-delimited taxon
    names to avoid partial replacements.
    """
    annotated = newick_str

    # Sort by decreasing length to avoid partial-match issues (e.g., "Sp1" vs "Sp10")
    for sp in sorted(fg_species, key=len, reverse=True):
        # Escape any special regex characters in the species name
        escaped = re.escape(sp)
        # Match the species name when followed by one of:,  )  :  [  {  ; (typical Newick delimiters)
        # but NOT already followed by {  to avoid double-annotation
        pattern = r"(?<![a-zA-Z0-9_])" + escaped + r"(?!\{)(?=[,):;\[\s]|$)"
        replacement = f"{sp}{{{label}}}"
        annotated = re.sub(pattern, replacement, annotated)

    return annotated


def main():
    parser = argparse.ArgumentParser(description="Annotate foreground branches for HyPhy FADE")
    parser.add_argument("--traitfile", required=True, help="caastools traitfile (3-col TSV)")
    parser.add_argument("--tree",      required=True, help="Newick tree file")
    parser.add_argument("--output",    required=True, help="Output annotated Newick file")
    parser.add_argument("--flip",      action="store_true",
                        help="Swap FG/BG: use contrast_group==0 as foreground (for BOTTOM direction)")
    parser.add_argument("--label",     default="Foreground",
                        help="HyPhy branch label to apply (default: Foreground)")
    args = parser.parse_args()

    fg_species = parse_traitfile(args.traitfile, flip=args.flip)
    if not fg_species:
        sys.exit("ERROR annotate_tree_fg: no foreground species found in traitfile "
                 f"(flip={args.flip}). Check contrast_group column values.")

    with open(args.tree) as fh:
        newick_str = fh.read().strip()

    # Validate that the tree is parseable
    try:
        tree = dendropy.Tree.get(data=newick_str, schema="newick",
                                 preserve_underscores=True)
    except Exception as exc:
        sys.exit(f"ERROR annotate_tree_fg: cannot parse tree: {exc}")

    # Check coverage: warn about species in traitfile not in tree
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    missing_fg = fg_species - tree_taxa
    if missing_fg:
        sys.stderr.write(
            f"WARNING annotate_tree_fg: {len(missing_fg)} FG species not found in tree "
            f"(will be ignored): {', '.join(sorted(missing_fg))}\n"
        )
    # Only annotate species that are actually in the tree
    fg_in_tree = fg_species & tree_taxa
    if not fg_in_tree:
        sys.exit("ERROR annotate_tree_fg: none of the foreground species are present in the tree.")

    annotated = annotate_newick(newick_str, fg_in_tree, label=args.label)

    with open(args.output, "w") as fh:
        fh.write(annotated)
        if not annotated.endswith("\n"):
            fh.write("\n")

    sys.stderr.write(
        f"INFO annotate_tree_fg: annotated {len(fg_in_tree)} foreground branches "
        f"with label '{{{args.label}}}' (flip={args.flip})\n"
    )


if __name__ == "__main__":
    main()
