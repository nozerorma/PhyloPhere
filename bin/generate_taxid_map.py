#!/usr/bin/env python3
# generate_taxid_map.py — Resolve the tip labels of a species tree to NCBI tax_ids.
# PhyloPhere | bin/

"""
GenerateTaxidMap: builds the tax_id/species table of the pipeline from the tip labels of a
Newick tree, for runs where the user supplies no tax_id file.

Resolution is exact-match only: a label, with underscores read as spaces, is looked up as an NCBI
scientific name. A name without an exact match is written to --unresolved, never guessed by
fuzzy or synonym matching. One narrow exception: a label with more than two tokens (e.g.
"Cyperus eragrostis FM208065", a species name plus a suffix that keeps sibling tips unique)
is retried as its first two tokens. This is still an exact match on a deterministically
derived name; several tips then share one tax_id, as they are the same species.

Source: live NCBI eutils, queried one name at a time in three tiers, so that a transient
failure on one name does not affect the others:
  1. each name is retried up to LIVE_ATTEMPTS times with exponential backoff on network
     errors, timeouts, HTTP errors or malformed replies;
  2. the names still failing are retried once as a group after SECOND_PASS_PAUSE seconds
     (rate limiting, short outages);
  3. only the names that never got a live answer fall back to the local NCBI taxonomy dump of
     ete3, each reported on stderr. The dump is a snapshot and can disagree with live NCBI
     for recently updated names, so it is the last resort.
A live answer of "no exact match" is final and not retried.

Called by:  resolve_core_inputs.py (when tax_id is blank and a tree is given)
Inputs:     --tree        Newick species tree
Outputs:    --output      TSV with columns tax_id, species (the label), the schema read by
                          subworkflows/TRAIT_ANALYSIS/local/src/phylo.R and
                          subworkflows/RERCONVERGE/local/rer_master_tree.R; unresolved labels are absent
            --unresolved  TSV with one column, species, listing the unresolved labels
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
LIVE_ATTEMPTS = 3          # tier 1: attempts per name
BACKOFF_BASE_S = 1.0       # tier 1: waits of 1, 2, 4 ... seconds between attempts
SECOND_PASS_PAUSE = 30.0   # tier 2: pause before re-querying the names that failed
RATE_LIMIT_S = 0.34        # eutils default rate limit, no API key configured
_TRANSIENT = (urllib.error.URLError, TimeoutError, OSError, ValueError, KeyError)


def normalize_label(label: str) -> str:
    """Strip quotes and surrounding blanks from a tip label and turn inner whitespace into underscores."""
    label = (label or "").strip().strip("'\"")
    return re.sub(r"\s+", "_", label)


def load_tip_labels(tree_path: str) -> list:
    """Unique normalized tip labels of the tree, in tree order."""
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    labels = [normalize_label(leaf.taxon.label) for leaf in tree.leaf_node_iter() if leaf.taxon]
    seen = set()
    unique = []
    for label in labels:
        if label not in seen:
            seen.add(label)
            unique.append(label)
    return unique


def _esearch_once(name: str) -> list[str]:
    """One eutils esearch for an exact scientific name; returns the list of tax_ids found.

    Raises one of _TRANSIENT when NCBI gives no usable answer (network error, timeout,
    HTTP error, malformed JSON).
    """
    term = urllib.parse.quote(f"{name}[Scientific Name]")
    with urllib.request.urlopen(f"{_EUTILS_ESEARCH}?db=taxonomy&retmode=json&term={term}",
                                 timeout=15) as r:
        return json.load(r)["esearchresult"]["idlist"]


def _resolve_names_live(names: list[str], attempts: int = LIVE_ATTEMPTS
                        ) -> tuple[dict[str, int], list[str]]:
    """Tier 1: live eutils per name with retries. Returns (resolved, failed).

    `failed` holds the names for which NCBI never gave a usable answer; a name answered
    without a unique exact match is in neither.
    """
    out: dict[str, int] = {}
    failed: list[str] = []
    for name in names:
        for attempt in range(attempts):
            try:
                ids = _esearch_once(name)
            except _TRANSIENT as exc:
                if attempt + 1 < attempts:
                    time.sleep(BACKOFF_BASE_S * 2 ** attempt)
                    continue
                print(f"  NCBI eutils failed for {name!r} after {attempts} attempts ({exc})",
                      file=sys.stderr)
                failed.append(name)
                break
            if len(ids) == 1:
                out[name] = int(ids[0])
            break
        time.sleep(RATE_LIMIT_S)
    return out, failed


def _resolve_names_local(names: list[str]) -> dict[str, int]:
    """Tier 3: the same exact-match resolution against the local NCBI taxonomy dump of ete3.

    May disagree with live NCBI if the dump predates a taxonomic update.
    """
    from ete3 import NCBITaxa

    name2taxid = NCBITaxa().get_name_translator(names)
    return {name: hits[0] for name, hits in name2taxid.items() if hits}


def _resolve_names(names: list[str]) -> dict[str, int]:
    """Three-tier resolution of `names` (see the module docstring); returns name → tax_id."""
    resolved, failed = _resolve_names_live(names)
    if failed:
        print(f"  {len(failed)} name(s) without a live NCBI answer; retrying in "
              f"{SECOND_PASS_PAUSE:.0f} s", file=sys.stderr)
        time.sleep(SECOND_PASS_PAUSE)
        second, failed = _resolve_names_live(failed)
        resolved.update(second)
    if failed:
        try:
            local = _resolve_names_local(failed)
        except Exception as exc:  # ete3 missing or its dump unavailable
            print(f"  local ete3 taxonomy unavailable ({exc}); "
                  f"{len(failed)} name(s) left unresolved", file=sys.stderr)
            local = {}
        for name in failed:
            if name in local:
                print(f"  resolved from local ete3 dump (no live NCBI answer): {name!r} -> "
                      f"{local[name]}", file=sys.stderr)
        resolved.update(local)
    return resolved


def resolve_taxids(labels: list[str]) -> tuple[dict[str, int], list[str]]:
    """Label → NCBI tax_id via _resolve_names: the full label first, then the
    genus + species fallback for labels with more than two tokens. Returns (resolved, unresolved)."""
    query_names = {label: label.replace("_", " ") for label in labels}
    name2taxid = _resolve_names(list(query_names.values()))

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
        fb2taxid = _resolve_names(list(set(fallback_query.values())))
        for label, fb_name in fallback_query.items():
            tax_id = fb2taxid.get(fb_name)
            if tax_id:
                resolved[label] = tax_id
                print(f"  resolved via genus+species fallback: {label!r} -> {fb_name!r}",
                      file=sys.stderr)

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
