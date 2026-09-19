#!/usr/bin/env python3
"""
resolve_gmts.py  —  Resolve --gmt_dir before invoking Nextflow, when left
blank: fetch fresh copies of the default GMT set, falling back to the
vendored copies (assets/gmt/) per-file when a fetch fails (offline node, no
network, upstream host down).

Fetching is tried first deliberately — offline safety is what the vendored
copies are for, not staleness avoidance, so a live run should always get the
current upstream gene sets when it can.

Why this runs outside main.nf, same reasoning as bin/resolve_core_inputs.py:
params.gmt_dir is read directly in multiple places inside
subworkflows/ENRICHMENT/fcs.nf, and Nextflow (25.x) enforces single-
assignment on params keys, so a params.gmt_dir = ... set inside workflow{}
after conf/enrichment.config's own default has already run would be
silently ignored.

Usage
-----
    resolve_gmts.py --output-dir <dir> [--vendored-dir <assets/gmt>] \
        [--gmt-dir <existing --gmt_dir value>] [--timeout 30]

Prints one line to stdout, shell-sourceable:
    GMT_DIR=<path>

If --gmt-dir is already set, it's echoed back unchanged (no fetch attempted
— an explicit override always wins).
"""

import argparse
import os
import re
import shutil
import sys
import urllib.request

_WIKIPATHWAYS_INDEX = "https://data.wikipathways.org/current/gmt/"

_SOURCES = {
    "go_biological_process.gmt":
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023",
    "go_molecular_function.gmt":
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2023",
    # wikipathways.gmt: resolved dynamically (see _resolve_wikipathways_url) —
    # the filename is date-stamped (wikipathways-YYYYMMDD-gmt-Homo_sapiens.gmt)
    # and changes each release, no stable "current" alias exists.
    #
    # Reactome ships as a zip archive (multi-file), not a plain-text GMT URL,
    # so it isn't included in the live-fetch set at all — the vendored copy is
    # always used for it. Refresh assets/gmt/reactome_pathways.gmt manually
    # (see assets/gmt/README.md) rather than teaching this script to unzip.
}


def _resolve_wikipathways_url(timeout: int):
    try:
        with urllib.request.urlopen(_WIKIPATHWAYS_INDEX, timeout=timeout) as resp:
            listing = resp.read().decode("utf-8", errors="replace")
    except Exception as exc:
        print(f"WARN: could not list {_WIKIPATHWAYS_INDEX} ({exc})", file=sys.stderr)
        return None
    match = re.search(r'wikipathways-\d{8}-gmt-Homo_sapiens\.gmt', listing)
    if not match:
        print("WARN: no Homo_sapiens GMT filename found in WikiPathways index", file=sys.stderr)
        return None
    return _WIKIPATHWAYS_INDEX + match.group(0)

_DEFAULT_VENDORED_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets", "gmt"
)


def resolve_gmts(output_dir: str, vendored_dir: str, timeout: int) -> None:
    os.makedirs(output_dir, exist_ok=True)

    vendored_files = {
        f for f in os.listdir(vendored_dir)
        if f.endswith(".gmt") and os.path.isfile(os.path.join(vendored_dir, f))
    } if os.path.isdir(vendored_dir) else set()

    all_fnames = vendored_files | set(_SOURCES) | {"wikipathways.gmt"}
    for fname in sorted(all_fnames):
        dest = os.path.join(output_dir, fname)
        if fname == "wikipathways.gmt":
            url = _resolve_wikipathways_url(timeout)
        else:
            url = _SOURCES.get(fname)
        fetched = False
        if url:
            try:
                with urllib.request.urlopen(url, timeout=timeout) as resp:
                    data = resp.read()
                if data and data.lstrip().startswith(b"<") is False and len(data) > 100:
                    with open(dest, "wb") as fh:
                        fh.write(data)
                    fetched = True
                    print(f"Fetched {fname} from {url}", file=sys.stderr)
            except Exception as exc:
                print(f"WARN: fetch failed for {fname} ({exc}); "
                      "falling back to the vendored copy.", file=sys.stderr)

        if not fetched:
            vendored_src = os.path.join(vendored_dir, fname)
            if os.path.isfile(vendored_src):
                shutil.copy(vendored_src, dest)
                print(f"Using vendored copy of {fname}", file=sys.stderr)
            else:
                print(f"WARN: no vendored fallback for {fname}; skipping.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--vendored-dir", default=_DEFAULT_VENDORED_DIR)
    parser.add_argument("--gmt-dir", default="", help="Existing --gmt_dir value, if any")
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    if args.gmt_dir:
        print(f"GMT_DIR={args.gmt_dir}")
        return

    resolve_gmts(args.output_dir, args.vendored_dir, args.timeout)
    print(f"GMT_DIR={args.output_dir}")


if __name__ == "__main__":
    main()
