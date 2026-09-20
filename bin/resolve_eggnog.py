#!/usr/bin/env python3
"""
resolve_eggnog.py  —  Resolve --egg_members_file / --egg_annotations_file
before invoking Nextflow, when either is left blank: fetch a fresh copy of
the eggNOG5 Primates-level (taxid 9443) orthogroup members/annotations pair,
filtered down to the human-relevant subset build_position_gmt.py actually
reads, falling back to the vendored copies (assets/eggnog/) when the fetch
fails (offline node, no network, upstream host down).

Fetching is tried first deliberately — offline safety is what the vendored
copies are for, not staleness avoidance, so a live run should always get the
current upstream orthogroups when it can.

Why this runs outside main.nf, same reasoning as bin/resolve_core_inputs.py /
bin/resolve_gmts.py: params.egg_members_file / params.egg_annotations_file
are read directly by subworkflows/ENRICHMENT/posenrich.nf's process inputs,
and Nextflow (25.x) enforces single-assignment on params keys, so a
params.egg_members_file = ... set inside workflow{} after
conf/enrichment.config's own default has already run would be silently
ignored.

Filtering (why): build_position_gmt.py:586-605 only ever reads the
orthogroup id (members col 2 / annotations col 2), the description
(annotations col 4), and members whose taxon prefix is "9606.ENSP" (human) —
every other member and annotation row is discarded on load. Of the 23,677
Primates-level orthogroups in the full upstream files, only the ones with at
least one human member matter; keeping only those rows (and only the human
members within each row) shrinks the pair from ~1.8MB to ~470KB compressed
without changing what the pipeline actually uses.

Usage
-----
    resolve_eggnog.py --output-dir <dir> [--vendored-dir <assets/eggnog>] \
        [--egg-members-file <existing --egg_members_file value>] \
        [--egg-annotations-file <existing --egg_annotations_file value>] \
        [--timeout 30]

Prints two lines to stdout, shell-sourceable:
    EGG_MEMBERS_FILE=<path>
    EGG_ANNOTATIONS_FILE=<path>

If both --egg-members-file and --egg-annotations-file are already set, they
are echoed back unchanged (no fetch attempted — an explicit override always
wins). Resolution otherwise always regenerates both together: the members and
annotations files are paired by orthogroup id, so resolving only one of them
if the other is blank would silently mismatch releases.
"""

import argparse
import gzip
import os
import shutil
import sys
import urllib.request

_MEMBERS_URL = "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_members.tsv.gz"
_ANNOTATIONS_URL = "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/9443/9443_annotations.tsv.gz"

_MEMBERS_FNAME = "9443_members_human.tsv.gz"
_ANNOTATIONS_FNAME = "9443_annotations_human.tsv.gz"

_DEFAULT_VENDORED_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets", "eggnog"
)


def _fetch(url: str, timeout: int) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read()


def _filter_and_write(members_raw: bytes, annotations_raw: bytes,
                       members_dest: str, annotations_dest: str) -> None:
    human_ogs = set()
    members_lines = []
    for line in gzip.decompress(members_raw).decode("utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        human_members = [m for m in fields[4].split(",") if m.startswith("9606.ENSP")]
        if human_members:
            human_ogs.add(fields[1])
            members_lines.append("\t".join([fields[0], fields[1], fields[2], fields[3],
                                             ",".join(human_members)]))

    annotations_lines = []
    for line in gzip.decompress(annotations_raw).decode("utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        if fields[1] in human_ogs:
            annotations_lines.append("\t".join(fields[:4]))

    with gzip.open(members_dest, "wt") as fh:
        fh.write("\n".join(members_lines) + "\n")
    with gzip.open(annotations_dest, "wt") as fh:
        fh.write("\n".join(annotations_lines) + "\n")


def resolve_eggnog(output_dir: str, vendored_dir: str, timeout: int) -> tuple:
    os.makedirs(output_dir, exist_ok=True)
    members_dest = os.path.join(output_dir, _MEMBERS_FNAME)
    annotations_dest = os.path.join(output_dir, _ANNOTATIONS_FNAME)

    fetched = False
    try:
        members_raw = _fetch(_MEMBERS_URL, timeout)
        annotations_raw = _fetch(_ANNOTATIONS_URL, timeout)
        _filter_and_write(members_raw, annotations_raw, members_dest, annotations_dest)
        fetched = True
        print(f"Fetched and filtered eggNOG 9443 members/annotations from {_MEMBERS_URL}",
              file=sys.stderr)
    except Exception as exc:
        print(f"WARN: eggNOG fetch/filter failed ({exc}); "
              "falling back to the vendored copy.", file=sys.stderr)

    if not fetched:
        vendored_members = os.path.join(vendored_dir, _MEMBERS_FNAME)
        vendored_annotations = os.path.join(vendored_dir, _ANNOTATIONS_FNAME)
        if os.path.isfile(vendored_members) and os.path.isfile(vendored_annotations):
            shutil.copy(vendored_members, members_dest)
            shutil.copy(vendored_annotations, annotations_dest)
            print("Using vendored copy of eggNOG 9443 members/annotations", file=sys.stderr)
        else:
            print("WARN: no vendored eggNOG fallback found; leaving files unresolved.",
                  file=sys.stderr)
            return "", ""

    return members_dest, annotations_dest


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--vendored-dir", default=_DEFAULT_VENDORED_DIR)
    parser.add_argument("--egg-members-file", default="", help="Existing --egg_members_file value, if any")
    parser.add_argument("--egg-annotations-file", default="", help="Existing --egg_annotations_file value, if any")
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    if args.egg_members_file and args.egg_annotations_file:
        print(f"EGG_MEMBERS_FILE={args.egg_members_file}")
        print(f"EGG_ANNOTATIONS_FILE={args.egg_annotations_file}")
        return

    members_file, annotations_file = resolve_eggnog(args.output_dir, args.vendored_dir, args.timeout)
    print(f"EGG_MEMBERS_FILE={members_file}")
    print(f"EGG_ANNOTATIONS_FILE={annotations_file}")


if __name__ == "__main__":
    main()
