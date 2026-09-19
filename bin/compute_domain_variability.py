#!/usr/bin/env python3
"""
compute_domain_variability.py  —  Auto-generate ENRICHMENT's
--domain_variability_file (Pfam domain-per-gene annotation used to build
pfam_domains.gmt/pfam_clans.gmt in build_position_gmt.py) directly from the
alignment, via a cached Pfam-A + hmmscan.

Note this is a DIFFERENT schema from ortholog_characterizator's own
domain_variability.tsv (map_domain_variability.py's output: per-domain-hit
*variability statistics* — mean/max/n_analyzed). phylophere's file instead
needs Pfam *metadata* per hit (target_name, description, clan_acc,
clan_name — see build_position_gmt.py:334-335) to build gene-set GMTs, no
variability numbers at all. So this reuses map_domain_variability.py's
domtblout parser (parse_domtblout) and reference-sequence extraction
(get_ref_seq) verbatim — that IS re-derived from ortholog_characterizator —
but the Pfam-A.clans.tsv metadata join below is new: phylophere doesn't
already have anything that does it, in either repo.

hmmscan is run ONCE across all genes' reference sequences pooled into one
FASTA (the pattern map_domain_variability.py's own docstring documents for
efficiency), not once per gene.

Cache ("cache large" pattern, not committed): --cache-dir (default
~/.cache/phylophere/pfam/) holds Pfam-A.hmm (+ hmmpress binary index) and
Pfam-A.clans.tsv, downloaded once from EBI's InterPro FTP and reused across
runs/genes.

Usage
-----
    compute_domain_variability.py --alignment-dir <dir> --output-dir <dir> \
        [--cache-dir ~/.cache/phylophere/pfam] [--ref-species Homo_sapiens] \
        [--evalue-threshold 0.01] [--alignment-format fasta]

Output schema (--output-dir/domain_variability.tsv), per
build_position_gmt.py:334-335:
    gene, pfam_id, target_name, description, clan_acc, clan_name, ali_start, ali_end
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from map_domain_variability import get_ref_seq, parse_domtblout, read_fasta  # noqa: E402

_PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
_PFAM_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"


def _download_gz(url: str, dest: str) -> None:
    print(f"Downloading {url} -> {dest}", file=sys.stderr)
    tmp_gz = dest + ".gz.tmp"
    urllib.request.urlretrieve(url, tmp_gz)
    with gzip.open(tmp_gz, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(tmp_gz)


def ensure_pfam_cache(cache_dir: str) -> tuple:
    os.makedirs(cache_dir, exist_ok=True)
    hmm_path = os.path.join(cache_dir, "Pfam-A.hmm")
    clans_path = os.path.join(cache_dir, "Pfam-A.clans.tsv")

    if not os.path.exists(hmm_path):
        _download_gz(_PFAM_HMM_URL, hmm_path)
    if not os.path.exists(hmm_path + ".h3f"):
        print(f"Running hmmpress on {hmm_path}", file=sys.stderr)
        subprocess.run(["hmmpress", hmm_path], check=True)

    if not os.path.exists(clans_path):
        _download_gz(_PFAM_CLANS_URL, clans_path)

    return hmm_path, clans_path


def load_clan_metadata(clans_tsv: str) -> dict:
    """pfam_acc (no version) -> {target_name, description, clan_acc, clan_name}.

    Pfam-A.clans.tsv columns (no header): pfam_acc, clan_acc, clan_id,
    pfam_id (short name), pfam_description.
    """
    meta = {}
    with open(clans_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            pfam_acc, clan_acc, clan_id, pfam_id, description = parts[:5]
            meta[pfam_acc] = {
                "target_name": pfam_id,
                "description": description,
                "clan_acc": clan_acc or "NA",
                "clan_name": clan_id or "NA",
            }
    return meta


def build_reference_fasta(alignment_dir: str, alignment_format: str, ref_species: str,
                          out_fasta: str) -> list:
    """Extract each gene's reference-species sequence into one combined FASTA.

    Only 'fasta' alignments are supported (map_domain_variability.py's
    read_fasta() is a hand-rolled FASTA-only parser, same constraint as the
    rest of this Valdar-variability toolchain).
    """
    if alignment_format != "fasta":
        sys.exit("Error: compute_domain_variability.py only supports --alignment-format fasta "
                  "(map_domain_variability.py's read_fasta() is FASTA-only).")

    genes = []
    with open(out_fasta, "w") as out:
        for fname in sorted(os.listdir(alignment_dir)):
            path = os.path.join(alignment_dir, fname)
            if not os.path.isfile(path):
                continue
            gene = os.path.splitext(fname)[0]
            ref = get_ref_seq(path, ref_species)
            if ref is None:
                print(f"WARN: reference species '{ref_species}' not found in {fname}",
                      file=sys.stderr)
                continue
            _, seq = ref
            ungapped = seq.replace("-", "").replace("X", "")
            out.write(f">{gene}\n{ungapped}\n")
            genes.append(gene)
    return genes


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cache-dir", default=os.path.expanduser("~/.cache/phylophere/pfam"))
    parser.add_argument("--ref-species", default="Homo_sapiens")
    parser.add_argument("--evalue-threshold", type=float, default=0.01)
    parser.add_argument("--alignment-format", default="fasta")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    hmm_path, clans_path = ensure_pfam_cache(args.cache_dir)
    clan_meta = load_clan_metadata(clans_path)

    ref_fasta = os.path.join(args.output_dir, "reference_seqs.fa")
    genes = build_reference_fasta(args.alignment_dir, args.alignment_format,
                                  args.ref_species, ref_fasta)
    if not genes:
        sys.exit(f"Error: no reference ('{args.ref_species}') sequences found in "
                  f"{args.alignment_dir}")

    domtblout = os.path.join(args.output_dir, "hmmscan.domtblout")
    print(f"Running hmmscan for {len(genes)} genes against {hmm_path}...", file=sys.stderr)
    subprocess.run(
        ["hmmscan", "--domtblout", domtblout, "--cpu", str(os.cpu_count() or 1),
         hmm_path, ref_fasta],
        check=True, stdout=subprocess.DEVNULL,
    )

    hits = parse_domtblout(domtblout, args.evalue_threshold)

    out_tsv = os.path.join(args.output_dir, "domain_variability.tsv")
    n_written = 0
    with open(out_tsv, "w") as fh:
        fh.write("gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\n")
        for hit in hits:
            meta = clan_meta.get(hit["pfam_id"])
            if not meta:
                print(f"WARN: no Pfam-A.clans.tsv entry for {hit['pfam_id']} "
                      f"(gene {hit['gene']}) — skipping", file=sys.stderr)
                continue
            fh.write(f"{hit['gene']}\t{hit['pfam_id']}\t{meta['target_name']}\t"
                     f"{meta['description']}\t{meta['clan_acc']}\t{meta['clan_name']}\t"
                     f"{hit['ali_start']}\t{hit['ali_end']}\n")
            n_written += 1

    print(f"Wrote {n_written} domain hits ({len(genes)} genes scanned) to {out_tsv}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
