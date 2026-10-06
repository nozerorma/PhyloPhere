#!/usr/bin/env python3
# resolve_gmts.py — Resolve gmt_dir, the directory of gene sets (*.gmt) read by the FCS enrichment.
# PhyloPhere | bin/

"""
ResolveGmts: provides the gene-set directory used when --gmt_dir is left blank.

Two modes, chosen explicitly:

  default   copy the *.gmt files versioned in subworkflows/ENRICHMENT/dat/ into the output directory.
  --fetch   the same copy, plus the current GO Biological Process and Molecular Function libraries (Enrichr) and the
            current WikiPathways human set, downloaded under the names go_biological_process.gmt,
            go_molecular_function.gmt and wikipathways.gmt. These are different collections from the versioned GO and
            WikiPathways files, so they are added to the set, not substituted for it. A failed download is an error:
            nothing is written and the versioned set is not offered as a replacement.

Both modes write gmt_source.json in the output directory: for each file its origin (versioned or fetched), its SHA-256
and, for a download, the URL and the UTC retrieval time.

Why it runs outside main.nf: params.gmt_dir is read directly in subworkflows/ENRICHMENT/fcs.nf, and Nextflow keeps the
first assignment of a params key, so a value set inside workflow {} after conf/enrichment.config would be ignored.

Called by:  the generated run script (gui/generation/templates/run_single.sh.j2), before nextflow
Inputs:     --output-dir, --fetch, --timeout, --versioned-dir (default subworkflows/ENRICHMENT/dat);
            --gmt-dir, an existing --gmt_dir value, wins and is echoed back unchanged
Outputs:    the *.gmt files and gmt_source.json in --output-dir;
            stdout, shell-sourceable:  GMT_DIR=<path>

Usage:
    resolve_gmts.py --output-dir <dir> [--fetch] [--versioned-dir <dir>] [--gmt-dir <existing --gmt_dir value>] [--timeout 30]
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import sys
import urllib.request

_WIKIPATHWAYS_INDEX = "https://data.wikipathways.org/current/gmt/"

# Downloaded files and where they come from. The WikiPathways file name carries the release date, so its URL is resolved
# from the index page.
_SOURCES = {
    "go_biological_process.gmt":
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023",
    "go_molecular_function.gmt":
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2023",
}
_WIKIPATHWAYS_FILE = "wikipathways.gmt"

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULT_VERSIONED_DIR = os.path.join(_REPO_ROOT, "subworkflows", "ENRICHMENT", "dat")


def _download(url: str, timeout: int) -> bytes:
    """Body of an HTTP GET; errors propagate to the caller."""
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read()


def _resolve_wikipathways_url(timeout: int) -> str:
    """URL of the current Homo_sapiens GMT, read from the WikiPathways index page (the file name carries the release date)."""
    try:
        listing = _download(_WIKIPATHWAYS_INDEX, timeout).decode("utf-8", errors="replace")
    except Exception as exc:
        raise RuntimeError(f"could not list {_WIKIPATHWAYS_INDEX} ({exc})") from exc
    match = re.search(r"wikipathways-\d{8}-gmt-Homo_sapiens\.gmt", listing)
    if not match:
        raise RuntimeError(f"no Homo_sapiens GMT file in the WikiPathways index {_WIKIPATHWAYS_INDEX}")
    return _WIKIPATHWAYS_INDEX + match.group(0)


def _sha256(data: bytes) -> str:
    """Hex SHA-256 of a byte string."""
    return hashlib.sha256(data).hexdigest()


def _download_all(timeout: int) -> dict:
    """{file name: (url, bytes)} of every downloaded set; RuntimeError if any one fails or is not a gene-set table."""
    urls = dict(_SOURCES)
    urls[_WIKIPATHWAYS_FILE] = _resolve_wikipathways_url(timeout)
    out = {}
    for name, url in urls.items():
        try:
            data = _download(url, timeout)
        except Exception as exc:
            raise RuntimeError(f"download of {name} from {url} failed ({exc})") from exc
        if len(data) <= 100 or data.lstrip().startswith(b"<"):
            raise RuntimeError(f"{url} did not return a gene-set table (got {len(data)} bytes, starts with {data[:20]!r})")
        out[name] = (url, data)
    return out


def resolve_gmts(output_dir: str, versioned_dir: str, timeout: int, fetch: bool = False) -> None:
    """Write the versioned *.gmt (and, with `fetch`, the downloaded ones) and `gmt_source.json` into `output_dir`."""
    if not os.path.isdir(versioned_dir):
        raise RuntimeError(f"versioned gene-set directory missing: {versioned_dir}")
    versioned = sorted(f for f in os.listdir(versioned_dir)
                       if f.endswith(".gmt") and os.path.isfile(os.path.join(versioned_dir, f)))
    downloaded = _download_all(timeout) if fetch else {}
    retrieved = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    os.makedirs(output_dir, exist_ok=True)
    files = {}
    for name in versioned:
        shutil.copy(os.path.join(versioned_dir, name), os.path.join(output_dir, name))
        with open(os.path.join(output_dir, name), "rb") as fh:
            files[name] = {"origin": "versioned", "sha256": _sha256(fh.read())}
    for name, (url, data) in downloaded.items():
        with open(os.path.join(output_dir, name), "wb") as fh:
            fh.write(data)
        files[name] = {"origin": "fetched", "sha256": _sha256(data), "url": url, "retrieved_utc": retrieved}

    with open(os.path.join(output_dir, "gmt_source.json"), "w") as fh:
        json.dump({"mode": "fetched" if fetch else "versioned", "versioned_dir": versioned_dir, "files": files},
                  fh, indent=2, sort_keys=True)
        fh.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--versioned-dir", default=_DEFAULT_VERSIONED_DIR)
    parser.add_argument("--gmt-dir", default="", help="existing --gmt_dir value, if any")
    parser.add_argument("--fetch", action="store_true", help="also download the current GO and WikiPathways sets")
    parser.add_argument("--timeout", type=int, default=30)
    args = parser.parse_args()

    if args.gmt_dir:
        print(f"GMT_DIR={args.gmt_dir}")
        return

    try:
        resolve_gmts(args.output_dir, args.versioned_dir, args.timeout, fetch=args.fetch)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    print(f"GMT_DIR={args.output_dir}")


if __name__ == "__main__":
    main()
