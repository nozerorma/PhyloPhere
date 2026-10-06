#!/usr/bin/env python3
# compute_domain_variability.py — Pfam domain hits per gene, with Pfam metadata, from the reference sequences of an alignment directory.
# PhyloPhere | bin/

"""
ComputeDomainVariability: scans the reference-species sequence of every gene against Pfam-A
with hmmscan and writes the domain table that build_position_gmt.py turns into the Pfam
domain and clan gene sets.

The table holds Pfam metadata per hit (name, description, clan), not variability statistics, so
its schema differs from the domain_variability.tsv of map_domain_variability.py. The hmmscan
domtblout parser and the reference-sequence extraction are imported from that script
(parse_domtblout, get_ref_seq); the join with Pfam-A.clans.tsv is done here. hmmscan runs once
on the reference sequences of all genes pooled in one FASTA.

Cache: --cache-dir (default ~/.cache/phylophere/pfam) holds Pfam-A.hmm (with its hmmpress index)
and Pfam-A.clans.tsv, downloaded once from the EBI Pfam FTP and reused. Requires hmmscan and
hmmpress on the PATH.

Called by:  COMPUTE_DOMAIN_VARIABILITY Nextflow process (subworkflows/ENRICHMENT/domain_variability_generation.nf → compute_domain_variability.py)
Inputs:     --alignment-dir     directory of FASTA alignments, one per gene (gene = file name without
                                extension); only fasta is supported
            --ref-species       species whose sequence is scanned (default Homo_sapiens); a gene
                                without it is skipped with a warning
            --evalue-threshold  maximum i-evalue of a domain hit (default 0.01)
Outputs:    <output-dir>/domain_variability.tsv  gene, pfam_id, target_name, description, clan_acc,
                                clan_name, ali_start, ali_end; ali_start and ali_end are positions
                                in the ungapped reference sequence, as reported by hmmscan.
                                Hits without a Pfam-A.clans.tsv entry are dropped with a warning.
            <output-dir>/reference_seqs.fa, hmmscan.domtblout  intermediate files

Usage:
    compute_domain_variability.py --alignment-dir <dir> --output-dir <dir> \
        [--cache-dir ~/.cache/phylophere/pfam] [--ref-species Homo_sapiens] \
        [--evalue-threshold 0.01] [--alignment-format fasta]
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
    """Download a gzip file and write it decompressed to dest."""
    print(f"Downloading {url} -> {dest}", file=sys.stderr)
    tmp_gz = dest + ".gz.tmp"
    urllib.request.urlretrieve(url, tmp_gz)
    with gzip.open(tmp_gz, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(tmp_gz)


def ensure_pfam_cache(cache_dir: str) -> tuple:
    """Make sure Pfam-A.hmm (hmmpress-indexed) and Pfam-A.clans.tsv are in cache_dir; returns their paths."""
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
    """Write the reference-species sequence of each gene, gaps removed, into one FASTA; returns the genes written.

    Only 'fasta' alignments are supported, because read_fasta() of map_domain_variability.py
    parses FASTA only.
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
