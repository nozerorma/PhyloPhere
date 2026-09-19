#!/usr/bin/env python3
"""
Map Pfam domain hits (from hmmscan) onto per-position Valdar variability data.

Prerequisite: run hmmscan on reference sequences extracted from PROT_RAW/:
  # Extract Homo_sapiens sequences from all PROT_RAW/*.fa files
  for f in PROT_RAW/*.fa; do
      gene=$(basename $f .fa)
      awk -v gene="$gene" '/^>Homo_sapiens/{found=1; print ">"gene; next}
           found && /^>/{exit} found{print}' "$f"
  done > ref_seqs.fa

  hmmscan --domtblout hmmscan.domtblout --cpu 8 /path/to/Pfam-A.hmm ref_seqs.fa

Output: domain_variability.tsv
  gene, pfam_id, domain_instance, ali_start, ali_end, n_positions, n_analyzed,
  mean_variability, max_variability, mean_gap_fraction

Coordinate system:
  - hmmscan ali_from/ali_to = 1-based ungapped sequence positions in the reference protein
  - seq_pos_to_ali_col maps those → 0-based MSA alignment columns
  - entropy.tsv position column = 1-based MSA column (0-based + 1)
"""

import argparse
import glob
import math
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_fasta(path):
    seqs = []
    header, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(buf)))
                header = line[1:]
                buf = []
            else:
                buf.append(line.upper())
    if header is not None:
        seqs.append((header, ''.join(buf)))
    return seqs


# ---------------------------------------------------------------------------
# hmmscan domtblout parser
# ---------------------------------------------------------------------------

def parse_domtblout(domtblout_path, evalue_threshold=0.01):
    """
    Parse hmmscan --domtblout output.
    Relevant columns (0-indexed, space-separated):
      0:  target_name  (HMM NAME field, e.g. 'Ig_2' — human-readable but not stable)
      1:  accession    (HMM ACC field, e.g. 'PF13895.13' → strip to 'PF13895')
      3:  query_name   (gene name)
      17: ali_from     (1-based start in query sequence)
      18: ali_to       (1-based end in query sequence)
      11: i-evalue
    Returns list of dicts; domain_instance is per-(gene, pfam_id) counter.
    """
    hits = []
    domain_counter = defaultdict(int)  # (gene, pfam_id) → instance count

    with open(domtblout_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 19:
                continue
            try:
                i_evalue = float(fields[11])
            except ValueError:
                continue
            if i_evalue > evalue_threshold:
                continue

            acc_raw = fields[1]
            pfam_id = acc_raw.split('.')[0]  # strip version: PF13895.13 → PF13895
            gene = fields[3]
            try:
                ali_start = int(fields[17])
                ali_end = int(fields[18])
            except ValueError:
                continue

            key = (gene, pfam_id)
            domain_counter[key] += 1
            hits.append({
                'pfam_id': pfam_id,
                'gene': gene,
                'ali_start': ali_start,
                'ali_end': ali_end,
                'i_evalue': i_evalue,
                'domain_instance': domain_counter[key],
            })

    return hits


# ---------------------------------------------------------------------------
# Reference sequence helpers
# ---------------------------------------------------------------------------

def get_ref_seq(fasta_path, ref_species='Homo_sapiens'):
    """
    Find the sequence whose header contains ref_species (case-insensitive substring).
    Returns (header, seq) or None.
    """
    seqs = read_fasta(fasta_path)
    pattern = ref_species.lower()
    for header, seq in seqs:
        if pattern in header.lower():
            return header, seq
    return None


def seq_pos_to_ali_col(prot_seq):
    """
    Build mapping: 1-based sequence position (non-gap char count) → 0-based alignment column.
    Gap characters: '-', 'X'.
    E.g. for 'A-BC': {1: 0, 2: 2, 3: 3}
    """
    mapping = {}
    seq_pos = 0
    for col, aa in enumerate(prot_seq):
        if aa not in ('-', 'X'):
            seq_pos += 1
            mapping[seq_pos] = col
    return mapping


# ---------------------------------------------------------------------------
# Entropy TSV loader
# ---------------------------------------------------------------------------

def load_entropy_tsv(tsv_path):
    """
    Read <gene>.entropy.tsv.
    Returns dict: {position_1based: {'t': float, 'r': float, 'g': float,
                                      'C_trident': float, 'variability': float}}
    position is the 1-based MSA column as written by compute_variability.py.
    """
    data = {}
    with open(tsv_path) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            try:
                pos = int(parts[1])
                data[pos] = {
                    't': float(parts[2]),
                    'r': float(parts[3]),
                    'g': float(parts[4]),
                    'C_trident': float(parts[5]),
                    'variability': float(parts[6]),
                }
            except (ValueError, IndexError):
                continue
    return data


# ---------------------------------------------------------------------------
# Domain statistics
# ---------------------------------------------------------------------------

def compute_domain_stats(ali_cols_0based, entropy_data):
    """
    Collect entropy rows for the given 0-based alignment columns.
    entropy_data uses 1-based keys → convert: key = ali_col_0based + 1.
    Returns dict with n_positions, n_analyzed, mean_variability,
    max_variability, mean_gap_fraction.
    """
    variabilities = []
    gap_fracs = []
    for col0 in ali_cols_0based:
        pos1 = col0 + 1
        if pos1 in entropy_data:
            row = entropy_data[pos1]
            variabilities.append(row['variability'])
            gap_fracs.append(row['g'])

    n_positions = len(ali_cols_0based)
    n_analyzed = len(variabilities)
    if n_analyzed == 0:
        return {
            'n_positions': n_positions, 'n_analyzed': 0,
            'mean_variability': float('nan'), 'max_variability': float('nan'),
            'mean_gap_fraction': float('nan'),
        }
    return {
        'n_positions': n_positions,
        'n_analyzed': n_analyzed,
        'mean_variability': sum(variabilities) / n_analyzed,
        'max_variability': max(variabilities),
        'mean_gap_fraction': sum(gap_fracs) / n_analyzed,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Map Pfam domain hits to per-position variability scores.')
    parser.add_argument('--var_dir', required=True,
                        help='PROT_VAR_RAW/ directory containing *.entropy.tsv files')
    parser.add_argument('--prot_raw_dir', required=True,
                        help='PROT_RAW/ directory containing protein MSA FASTA files')
    parser.add_argument('--hmmscan_domtblout', required=True,
                        help='hmmscan --domtblout output file')
    parser.add_argument('--out_tsv', required=True,
                        help='Output TSV path (domain_variability.tsv)')
    parser.add_argument('--ref_species', default='Homo_sapiens',
                        help='Reference species name for coordinate mapping (default: Homo_sapiens)')
    parser.add_argument('--evalue_threshold', type=float, default=0.01,
                        help='i-evalue cutoff for domain hits (default 0.01)')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out_tsv) or '.', exist_ok=True)

    hits = parse_domtblout(args.hmmscan_domtblout, args.evalue_threshold)
    if not hits:
        print(f'WARNING: no domain hits found in {args.hmmscan_domtblout}', file=sys.stderr)

    # Group hits by gene for efficient processing
    hits_by_gene = defaultdict(list)
    for h in hits:
        hits_by_gene[h['gene']].append(h)

    rows = []
    for gene, gene_hits in sorted(hits_by_gene.items()):
        prot_fa = os.path.join(args.prot_raw_dir, f'{gene}.fa')
        ent_tsv = os.path.join(args.var_dir, f'{gene}.entropy.tsv')

        if not os.path.isfile(prot_fa):
            print(f'WARN {gene}: PROT_RAW FASTA not found at {prot_fa}', file=sys.stderr)
            continue
        if not os.path.isfile(ent_tsv):
            print(f'WARN {gene}: entropy TSV not found at {ent_tsv}', file=sys.stderr)
            continue

        ref = get_ref_seq(prot_fa, args.ref_species)
        if ref is None:
            print(f'WARN {gene}: reference species "{args.ref_species}" not found in {prot_fa}',
                  file=sys.stderr)
            continue
        _, ref_seq = ref

        pos_map = seq_pos_to_ali_col(ref_seq)  # 1-based seq pos → 0-based ali col
        entropy_data = load_entropy_tsv(ent_tsv)

        for hit in gene_hits:
            ali_start = hit['ali_start']
            ali_end = hit['ali_end']
            # Collect 0-based alignment columns for this domain range
            ali_cols = [pos_map[p] for p in range(ali_start, ali_end + 1) if p in pos_map]
            if not ali_cols:
                print(f'WARN {gene} {hit["pfam_id"]}: no alignment columns for '
                      f'domain range {ali_start}-{ali_end}', file=sys.stderr)
                continue

            stats = compute_domain_stats(ali_cols, entropy_data)
            mean_v = stats['mean_variability']
            max_v  = stats['max_variability']
            mean_g = stats['mean_gap_fraction']

            rows.append((
                gene,
                hit['pfam_id'],
                hit['domain_instance'],
                ali_start, ali_end,
                stats['n_positions'],
                stats['n_analyzed'],
                f'{mean_v:.6f}' if not math.isnan(mean_v) else 'NA',
                f'{max_v:.6f}'  if not math.isnan(max_v)  else 'NA',
                f'{mean_g:.6f}' if not math.isnan(mean_g) else 'NA',
            ))

    with open(args.out_tsv, 'w') as fh:
        fh.write('gene\tpfam_id\tdomain_instance\tali_start\tali_end\t'
                 'n_positions\tn_analyzed\tmean_variability\tmax_variability\t'
                 'mean_gap_fraction\n')
        for row in rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')

    print(f'Written: {args.out_tsv} ({len(rows)} domain hits)', flush=True)


if __name__ == '__main__':
    main()
