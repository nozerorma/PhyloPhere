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


def annotate_newick(newick_str: str, fg_species: set, label: str) -> tuple:
    """
    Annotate foreground leaf branch names in a Newick string with {label}.

    Handles two taxon-name formats:
      1. Single-quoted (Dendropy re-serialisation): 'Name'  → 'Name'{label}
         The label is placed AFTER the closing quote so HyPhy reads it correctly.
      2. Bare / unquoted:                           Name    → Name{label}

    Returns (annotated_string, substitution_count).
    """
    annotated = newick_str
    total_subs = 0

    # Sort by decreasing length to avoid partial-match issues (e.g., "Sp1" vs "Sp10")
    for sp in sorted(fg_species, key=len, reverse=True):
        # Escape any special regex characters in the species name
        escaped = re.escape(sp)

        # --- Pattern 1: single-quoted name (Dendropy output) ---
        # Matches: 'SpeciesName' followed by a Newick delimiter
        # Replacement puts {label} AFTER the closing quote: 'Name'{label}
        pattern_quoted = r"'" + escaped + r"'(?!\{)(?=[,):;\[\s]|$)"
        replacement_quoted = f"'{sp}'{{{label}}}"
        annotated, n = re.subn(pattern_quoted, replacement_quoted, annotated)

        if n == 0:
            # --- Pattern 2: unquoted / bare name ---
            # Negative lookbehind excludes [a-zA-Z0-9_'] so we don't match
            # a substring inside a quoted name.
            pattern_bare = r"(?<![a-zA-Z0-9_'])" + escaped + r"(?!\{)(?=[,):;\[\s]|$)"
            replacement_bare = f"{sp}{{{label}}}"
            annotated, n = re.subn(pattern_bare, replacement_bare, annotated)

        total_subs += n

    return annotated, total_subs


def parse_fasta_taxa(path: str) -> set:
    """Return the set of taxon names from a FASTA file (header lines, stripped of '>')."""
    taxa = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                taxa.add(line[1:].split()[0])
    return taxa


def write_filtered_fasta(in_path: str, out_path: str, keep_taxa: set) -> int:
    """
    Write a copy of *in_path* containing only sequences whose header name
    (first whitespace-delimited token after '>') is present in *keep_taxa*.
    Returns the number of sequences written.
    """
    written = 0
    include = False
    with open(in_path) as fh_in, open(out_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                taxon = line[1:].strip().split()[0]
                include = taxon in keep_taxa
                if include:
                    written += 1
            if include:
                fh_out.write(line)
    return written


def main():
    parser = argparse.ArgumentParser(description="Annotate foreground branches for HyPhy FADE")
    parser.add_argument("--traitfile", required=True, help="caastools traitfile (3-col TSV)")
    parser.add_argument("--tree",      required=True, help="Newick tree file")
    parser.add_argument("--output",    required=True, help="Output annotated Newick file")
    parser.add_argument("--fasta",     default=None,
                        help="FASTA alignment; if provided the tree is pruned to only "
                             "taxa present in the alignment before annotation.")
    parser.add_argument("--fasta_out", default=None,
                        help="If provided, write a FASTA filtered to only sequences "
                             "present in the (pruned) tree.  Required when FASTA taxa "
                             "exceed tree tips to prevent HyPhy alignment/tree mismatch.")
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

    # Optional: prune tree to taxa present in the FASTA alignment.
    # This prevents FADE from failing when the alignment has fewer sequences
    # than tree tips (e.g. genes absent in some species).
    if args.fasta:
        fasta_taxa = parse_fasta_taxa(args.fasta)
        tree_taxa_all = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
        taxa_to_prune = tree_taxa_all - fasta_taxa
        if taxa_to_prune:
            sys.stderr.write(
                f"INFO annotate_tree_fg: pruning {len(taxa_to_prune)} tip(s) not present "
                f"in the FASTA alignment (tree: {len(tree_taxa_all)}, "
                f"fasta: {len(fasta_taxa)}).\n"
            )
            taxa_to_retain = dendropy.TaxonNamespace(
                [t for t in tree.taxon_namespace if t.label in fasta_taxa]
            )
            tree.retain_taxa(taxa_to_retain)
            tree.purge_taxon_namespace()
            newick_str = tree.as_string(schema="newick",
                                        suppress_rooting=True).strip()

    # Check coverage: warn about species in traitfile not in tree
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter()}

    # Optional: write a FASTA filtered to sequences present in the (pruned) tree.
    # This handles the reverse case: alignment has MORE sequences than tree tips,
    # which would otherwise cause HyPhy FADE to abort with a tip-count mismatch.
    if args.fasta_out and args.fasta:
        fasta_taxa_check = parse_fasta_taxa(args.fasta)
        extra_in_fasta = fasta_taxa_check - tree_taxa
        if extra_in_fasta:
            sys.stderr.write(
                f"INFO annotate_tree_fg: filtering {len(extra_in_fasta)} sequence(s) "
                f"present in FASTA but absent from tree — writing filtered FASTA "
                f"({len(fasta_taxa_check) - len(extra_in_fasta)} seqs) to "
                f"{args.fasta_out}\n"
            )
        n_written = write_filtered_fasta(args.fasta, args.fasta_out, tree_taxa)
        sys.stderr.write(
            f"INFO annotate_tree_fg: wrote {n_written} sequence(s) to filtered FASTA "
            f"{args.fasta_out}\n"
        )

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

    annotated, n_subs = annotate_newick(newick_str, fg_in_tree, label=args.label)

    if n_subs == 0:
        sys.exit(
            f"ERROR annotate_tree_fg: regex found 0 substitutions despite {len(fg_in_tree)} "
            f"foreground species in tree. Check that taxon names in the tree match the traitfile "
            f"exactly (including quoting style).")

    with open(args.output, "w") as fh:
        fh.write(annotated)
        if not annotated.endswith("\n"):
            fh.write("\n")

    sys.stderr.write(
        f"INFO annotate_tree_fg: annotated {n_subs} foreground branches "
        f"with label '{{{args.label}}}' (flip={args.flip})\n"
    )


if __name__ == "__main__":
    main()
