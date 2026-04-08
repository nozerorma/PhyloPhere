#!/usr/bin/env python3
"""
Annotate a Newick tree with `{Foreground}` labels using a plain species list.
"""

import argparse
import re
import sys

import dendropy


def parse_species_file(path: str) -> set:
    fg_species = set()
    with open(path) as fh:
        for line in fh:
            species = line.strip()
            if species and not species.startswith("#"):
                fg_species.add(species)
    return fg_species


def annotate_newick(newick_str: str, fg_species: set, label: str) -> tuple:
    annotated = newick_str
    total_subs = 0

    for sp in sorted(fg_species, key=len, reverse=True):
        escaped = re.escape(sp)

        pattern_quoted = r"'" + escaped + r"'(?!\{)(?=[,):;\[\s]|$)"
        replacement_quoted = f"'{sp}'{{{label}}}"
        annotated, n = re.subn(pattern_quoted, replacement_quoted, annotated)

        if n == 0:
            pattern_bare = r"(?<![a-zA-Z0-9_'])" + escaped + r"(?!\{)(?=[,):;\[\s]|$)"
            replacement_bare = f"{sp}{{{label}}}"
            annotated, n = re.subn(pattern_bare, replacement_bare, annotated)

        total_subs += n

    return annotated, total_subs


def parse_fasta_taxa(path: str) -> set:
    taxa = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                taxa.add(line[1:].split()[0])
    return taxa


def write_filtered_fasta(in_path: str, out_path: str, keep_taxa: set) -> int:
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
    parser.add_argument(
        "--species-file",
        required=True,
        help="Plain-text species list, one species name per line",
    )
    parser.add_argument("--tree", required=True, help="Newick tree file")
    parser.add_argument("--output", required=True, help="Output annotated Newick file")
    parser.add_argument(
        "--fasta",
        default=None,
        help="FASTA alignment; if provided the tree is pruned to only taxa present in the alignment before annotation.",
    )
    parser.add_argument(
        "--fasta_out",
        default=None,
        help="If provided, write a FASTA filtered to only sequences present in the (pruned) tree.",
    )
    parser.add_argument(
        "--label",
        default="Foreground",
        help="HyPhy branch label to apply (default: Foreground)",
    )
    args = parser.parse_args()

    fg_species = parse_species_file(args.species_file)
    if not fg_species:
        sys.exit(
            f"ERROR annotate_tree_fg: no foreground species found in '{args.species_file}'"
        )

    with open(args.tree) as fh:
        newick_str = fh.read().strip()

    try:
        tree = dendropy.Tree.get(
            data=newick_str, schema="newick", preserve_underscores=True
        )
    except Exception as exc:
        sys.exit(f"ERROR annotate_tree_fg: cannot parse tree: {exc}")

    if args.fasta:
        fasta_taxa = parse_fasta_taxa(args.fasta)
        tree_taxa_all = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
        taxa_to_prune = tree_taxa_all - fasta_taxa
        if taxa_to_prune:
            sys.stderr.write(
                f"INFO annotate_tree_fg: pruning {len(taxa_to_prune)} tip(s) not present "
                f"in the FASTA alignment (tree: {len(tree_taxa_all)}, fasta: {len(fasta_taxa)}).\n"
            )
            taxa_to_retain = dendropy.TaxonNamespace(
                [t for t in tree.taxon_namespace if t.label in fasta_taxa]
            )
            tree.retain_taxa(taxa_to_retain)
            tree.purge_taxon_namespace()
            newick_str = tree.as_string(
                schema="newick", suppress_rooting=True
            ).strip()

    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter()}

    if args.fasta_out and args.fasta:
        fasta_taxa_check = parse_fasta_taxa(args.fasta)
        extra_in_fasta = fasta_taxa_check - tree_taxa
        if extra_in_fasta:
            sys.stderr.write(
                f"INFO annotate_tree_fg: filtering {len(extra_in_fasta)} sequence(s) "
                f"present in FASTA but absent from tree; writing filtered FASTA to "
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

    fg_in_tree = fg_species & tree_taxa
    if not fg_in_tree:
        sys.stderr.write(
            "WARNING annotate_tree_fg: none of the foreground species are present "
            "after tree/alignment reconciliation; skipping this gene.\n"
        )
        return

    annotated, n_subs = annotate_newick(newick_str, fg_in_tree, args.label)
    if n_subs == 0:
        sys.stderr.write(
            f"WARNING annotate_tree_fg: regex found 0 substitutions despite {len(fg_in_tree)} "
            "foreground species in tree; skipping this gene.\n"
        )
        return

    # Remove all single quotes around species names (dendropy may add them).
    # This includes patterns like 'species'{Foreground} and 'species':branch_length
    annotated = re.sub(r"'([^']+)'(\{|:)", r"\1\2", annotated)

    with open(args.output, "w") as fh:
        fh.write(annotated)
        if not annotated.endswith("\n"):
            fh.write("\n")

    sys.stderr.write(
        f"INFO annotate_tree_fg: annotated {n_subs} foreground branches "
        f"with label '{{{args.label}}}'\n"
    )


if __name__ == "__main__":
    main()
