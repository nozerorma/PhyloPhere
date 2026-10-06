#!/usr/bin/env python3
# resolve_eggnog.py — Resolve the eggNOG orthogroup pair (members, annotations) read by POSENRICH.
# PhyloPhere | bin/

"""
ResolveEggnog: provides the eggNOG members and annotations files used when
--egg_members_file / --egg_annotations_file are left blank.

Two sources, chosen explicitly:

  default   the copy versioned in subworkflows/ENRICHMENT/dat/ (eggNOG 5.0, Primates level 9443, human members only).
            It exists for tax level 9443 with reference species 9606; any other pair is an error unless --fetch is given.
  --fetch   download the pair of the requested tax level from eggnog5.embl.de, keep the rows with a member of the
            reference species, and write them. A failed download is an error: the versioned copy is never substituted
            silently, because the two can hold different orthogroups.

Both modes write eggnog_source.json next to the pair: mode, tax level, reference species, the SHA-256 of each output
file and, for --fetch, the source URLs, the UTC retrieval time and the SHA-256 of each download.

Filtering: build_position_gmt.py reads the orthogroup id (members column 2, annotations column 2), the description
(annotations column 4) and the members whose taxon prefix is "<ref_taxid>." (human: "9606.ENSP"). Every other member
and annotation row is discarded on load, so only orthogroups with at least one reference member are kept, with only
those members in each row. For Primates this reduces ~1.8 MB to ~470 KB compressed and changes nothing the pipeline uses.

Called by:  RESOLVE_EGGNOG Nextflow process (subworkflows/ENRICHMENT/eggnog_resolution.nf → resolve_eggnog.py)
Inputs:     --output-dir, --tax-level (default 9443), --ref-taxid (default 9606), --fetch, --timeout,
            --versioned-dir (default subworkflows/ENRICHMENT/dat); with both --egg-members-file and
            --egg-annotations-file the two paths are echoed back unchanged
Outputs:    <output-dir>/<tax>_members_<ref>.tsv.gz, <tax>_annotations_<ref>.tsv.gz and eggnog_source.json;
            stdout, shell-sourceable:  EGG_MEMBERS_FILE=<path>  and  EGG_ANNOTATIONS_FILE=<path>

The two files are always produced together: the pair is joined by orthogroup id, and files of
different releases would mismatch.

Usage:
    resolve_eggnog.py --output-dir <dir> [--tax-level 9443] [--ref-taxid 9606] [--fetch] [--timeout 30] \
        [--versioned-dir <dir>] [--egg-members-file F --egg-annotations-file F]
"""

import argparse
import datetime
import gzip
import hashlib
import json
import os
import shutil
import sys
import urllib.request

_EGGNOG_URL = "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/{tax}/{tax}_{kind}.tsv.gz"

# The versioned pair covers one (tax level, reference species) combination.
_VERSIONED_TAX_LEVEL = "9443"
_VERSIONED_REF_TAXID = "9606"
_VERSIONED_MEMBERS = "9443_members_human.tsv.gz"
_VERSIONED_ANNOTATIONS = "9443_annotations_human.tsv.gz"

_DEFAULT_VERSIONED_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "subworkflows", "ENRICHMENT", "dat"
)


def _sha256_bytes(data: bytes) -> str:
    """Hex SHA-256 of a byte string."""
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: str) -> str:
    """Hex SHA-256 of a file, read in 1 MiB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch(url: str, timeout: int) -> bytes:
    """Body of an HTTP GET; errors propagate to the caller."""
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read()


def _filter_and_write(members_raw: bytes, annotations_raw: bytes,
                      members_dest: str, annotations_dest: str,
                      ref_taxid: str = "9606") -> None:
    """Keep the orthogroups with a member of `ref_taxid`, and only those members, in both files."""
    matched_ogs = set()
    members_lines = []
    prefix = f"{ref_taxid}."
    for line in gzip.decompress(members_raw).decode("utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        matched_members = [m for m in fields[4].split(",") if m.startswith(prefix)]
        if matched_members:
            matched_ogs.add(fields[1])
            members_lines.append("\t".join([fields[0], fields[1], fields[2], fields[3], ",".join(matched_members)]))

    annotations_lines = []
    for line in gzip.decompress(annotations_raw).decode("utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        if fields[1] in matched_ogs:
            annotations_lines.append("\t".join(fields[:4]))

    with gzip.open(members_dest, "wt") as fh:
        fh.write("\n".join(members_lines) + "\n")
    with gzip.open(annotations_dest, "wt") as fh:
        fh.write("\n".join(annotations_lines) + "\n")


def resolve_eggnog(output_dir: str, versioned_dir: str, timeout: int, tax_level: str = "9443",
                   ref_taxid: str = "9606", fetch: bool = False) -> tuple:
    """Write the pair and `eggnog_source.json` into `output_dir`; return the two paths. Raises RuntimeError."""
    os.makedirs(output_dir, exist_ok=True)
    members_dest = os.path.join(output_dir, f"{tax_level}_members_{ref_taxid}.tsv.gz")
    annotations_dest = os.path.join(output_dir, f"{tax_level}_annotations_{ref_taxid}.tsv.gz")
    source = {"tax_level": tax_level, "ref_taxid": ref_taxid}

    if fetch:
        urls = {kind: _EGGNOG_URL.format(tax=tax_level, kind=kind) for kind in ("members", "annotations")}
        try:
            raw = {kind: _fetch(url, timeout) for kind, url in urls.items()}
            _filter_and_write(raw["members"], raw["annotations"], members_dest, annotations_dest, ref_taxid=ref_taxid)
        except Exception as exc:
            raise RuntimeError(f"eggNOG download or filtering failed ({exc}); the versioned copy is not substituted. "
                               "Fix the network or set auto_fetch_eggnog=false.") from exc
        source.update(mode="fetched", urls=urls, retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                      download_sha256={kind: _sha256_bytes(data) for kind, data in raw.items()})
    else:
        if (tax_level, ref_taxid) != (_VERSIONED_TAX_LEVEL, _VERSIONED_REF_TAXID):
            raise RuntimeError(f"no versioned eggNOG copy for tax level {tax_level} and reference species {ref_taxid} "
                               f"(only {_VERSIONED_TAX_LEVEL} and {_VERSIONED_REF_TAXID}). Set auto_fetch_eggnog=true "
                               "(needs network) or give egg_members_file and egg_annotations_file.")
        members_src = os.path.join(versioned_dir, _VERSIONED_MEMBERS)
        annotations_src = os.path.join(versioned_dir, _VERSIONED_ANNOTATIONS)
        for path in (members_src, annotations_src):
            if not os.path.isfile(path):
                raise RuntimeError(f"versioned eggNOG file missing: {path}")
        shutil.copy(members_src, members_dest)
        shutil.copy(annotations_src, annotations_dest)
        source.update(mode="versioned", versioned_dir=versioned_dir)

    source["output_sha256"] = {"members": _sha256_file(members_dest), "annotations": _sha256_file(annotations_dest)}
    with open(os.path.join(output_dir, "eggnog_source.json"), "w") as fh:
        json.dump(source, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return members_dest, annotations_dest


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--versioned-dir", default=_DEFAULT_VERSIONED_DIR)
    parser.add_argument("--tax-level", default=_VERSIONED_TAX_LEVEL, help="eggNOG clade taxon ID")
    parser.add_argument("--ref-taxid", default=_VERSIONED_REF_TAXID, help="reference species NCBI taxon ID")
    parser.add_argument("--fetch", action="store_true", help="download instead of using the versioned copy")
    parser.add_argument("--egg-members-file", default="", help="existing --egg_members_file value, if any")
    parser.add_argument("--egg-annotations-file", default="", help="existing --egg_annotations_file value, if any")
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    if args.egg_members_file and args.egg_annotations_file:
        print(f"EGG_MEMBERS_FILE={args.egg_members_file}")
        print(f"EGG_ANNOTATIONS_FILE={args.egg_annotations_file}")
        return

    try:
        members_file, annotations_file = resolve_eggnog(args.output_dir, args.versioned_dir, args.timeout,
                                                        tax_level=args.tax_level, ref_taxid=args.ref_taxid, fetch=args.fetch)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    print(f"EGG_MEMBERS_FILE={members_file}")
    print(f"EGG_ANNOTATIONS_FILE={annotations_file}")


if __name__ == "__main__":
    main()
