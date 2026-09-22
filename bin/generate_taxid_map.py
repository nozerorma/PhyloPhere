#!/usr/bin/env python3
"""
generate_taxid_map.py  —  Auto-generate a tax_id/species mapping from a
species tree's tip labels, using NCBI taxonomy for resolution.

Replaces the need for a user-supplied --tax_id file in the common case
(species names already match, or differ only by underscore/space).
Resolution is exact-match only: a tip label is looked up as its scientific
name (underscores -> spaces) against NCBI taxonomy. Anything that doesn't
resolve exactly is reported in --unresolved rather than guessed via fuzzy
or synonym matching, so the user can supply corrections deliberately -- with
one narrow exception: a label with more than two whitespace-separated
tokens (e.g. "Cyperus eragrostis FM208065" -- a species name plus a
disambiguating suffix, the shape multi-accession-per-species fixtures use to
keep sibling tips unique) falls back to an exact match on just its first two
tokens ("Cyperus eragrostis"). This is still exact match, just against a
shorter, deterministically-derived candidate name -- not fuzzy/synonym
matching -- and multiple tips correctly collapsing to the same species-level
tax_id is expected in that case (they *are* the same species).

Resolution source: live NCBI eutils first (always current), falling back to
ete3's local cached NCBI taxonomy dump only if eutils itself can't be
reached at all (network outage / NCBI down) -- not merely because a
particular name didn't resolve. The local dump is a point-in-time snapshot
and can disagree with live NCBI for recently-updated names (observed for
e.g. "Machaerina articulata"), so it's the fallback, not the primary source.

Output (--output): TSV with columns tax_id, species — the exact schema
required by subworkflows/TRAIT_ANALYSIS/local/src/phylo.R and
subworkflows/RERCONVERGE/local/rer_master_tree.R.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import dendropy

_EUTILS_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def normalize_label(label: str) -> str:
    label = (label or "").strip().strip("'\"")
    return re.sub(r"\s+", "_", label)


def load_tip_labels(tree_path: str) -> list:
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    labels = [normalize_label(leaf.taxon.label) for leaf in tree.leaf_node_iter() if leaf.taxon]
    seen = set()
    unique = []
    for label in labels:
        if label not in seen:
            seen.add(label)
            unique.append(label)
    return unique


def _resolve_names_live(names: list[str]) -> dict[str, int]:
    """One eutils esearch per name, exact scientific-name match. Raises
    (urllib.error.URLError, TimeoutError, OSError) if NCBI can't be reached
    at all -- the caller decides whether to fall back to the local dump."""
    out: dict[str, int] = {}
    for name in names:
        term = urllib.parse.quote(f"{name}[Scientific Name]")
        with urllib.request.urlopen(f"{_EUTILS_ESEARCH}?db=taxonomy&retmode=json&term={term}",
                                     timeout=15) as r:
            ids = json.load(r)["esearchresult"]["idlist"]
        if len(ids) == 1:
            out[name] = int(ids[0])
        time.sleep(0.34)  # eutils default rate limit, no API key configured
    return out


def _resolve_names_local(names: list[str]) -> dict[str, int]:
    """Same exact-match resolution against ete3's local cached NCBI taxonomy
    dump -- used only when live eutils is unreachable. May disagree with
    live NCBI if the local dump predates a taxonomic update."""
    from ete3 import NCBITaxa

    name2taxid = NCBITaxa().get_name_translator(names)
    return {name: hits[0] for name, hits in name2taxid.items() if hits}


def resolve_taxids(labels: list[str]) -> tuple[dict[str, int], list[str]]:
    """label -> NCBI tax_id, live NCBI eutils first, falling back to ete3's
    local taxonomy dump only if NCBI itself can't be reached (network/API
    failure, not just an unresolved name). Within whichever source is used,
    tries the full label first, then the genus+species fallback described
    above for multi-token labels. Returns (resolved, unresolved).
    """
    query_names = {label: label.replace("_", " ") for label in labels}

    try:
        resolve_fn = _resolve_names_live
        name2taxid = resolve_fn(list(query_names.values()))
        source = "live NCBI"
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        print(f"NCBI eutils unreachable ({exc}) -- falling back to local ete3 taxonomy",
              file=sys.stderr)
        resolve_fn = _resolve_names_local
        name2taxid = resolve_fn(list(query_names.values()))
        source = "local ete3"

    resolved: dict[str, int] = {}
    fallback_query: dict[str, str] = {}
    for label, query_name in query_names.items():
        tax_id = name2taxid.get(query_name)
        if tax_id:
            resolved[label] = tax_id
            continue
        tokens = query_name.split()
        if len(tokens) > 2:
            fallback_query[label] = " ".join(tokens[:2])

    if fallback_query:
        fb2taxid = resolve_fn(list(set(fallback_query.values())))
        for label, fb_name in fallback_query.items():
            tax_id = fb2taxid.get(fb_name)
            if tax_id:
                resolved[label] = tax_id
                print(f"  resolved via genus+species fallback ({source}): "
                      f"{label!r} -> {fb_name!r}", file=sys.stderr)

    unresolved = [label for label in labels if label not in resolved]
    return resolved, unresolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Newick species tree")
    parser.add_argument("--output", required=True, help="Output TSV: tax_id, species")
    parser.add_argument("--unresolved", required=True,
                         help="Output TSV of tip labels that could not be resolved exactly")
    args = parser.parse_args()

    labels = load_tip_labels(args.tree)
    if not labels:
        print(f"Error: no tip labels found in {args.tree}", file=sys.stderr)
        sys.exit(1)

    resolved, unresolved = resolve_taxids(labels)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["tax_id", "species"])
        for label in labels:
            if label in resolved:
                writer.writerow([resolved[label], label])

    with open(args.unresolved, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["species"])
        for species in unresolved:
            writer.writerow([species])

    print(f"Resolved {len(resolved)}/{len(labels)} tip labels against NCBI taxonomy.",
          file=sys.stderr)
    if unresolved:
        print(f"WARNING: {len(unresolved)} tip labels unresolved, see {args.unresolved}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
