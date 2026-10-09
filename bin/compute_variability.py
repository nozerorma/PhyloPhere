#!/usr/bin/env python3
# compute_variability.py — Valdar (2002) C_trident positional conservation of one protein alignment.
# PhyloPhere | bin/

"""
ComputeVariability: scores every column of a protein alignment with the three-pronged Valdar
conservation score and summarizes variability per clade.

C_trident(x) = (1 - t(x))^alpha * (1 - r(x))^beta * (1 - g(x))^gamma   (Valdar 2002, Eq. 49)
  t(x): symbol diversity, weighted Shannon entropy over 20 amino acids + gap (Eq. 50-52)
  r(x): stereochemical diversity, mean distance from the consensus in BLOSUM62 20D space (Eq. 54-56)
  g(x): gap fraction
Variability = 1 - C_trident. Sequence weights are those of Henikoff and Henikoff (1994).

Called by:  compute_alignment_entropy.py (one call per gene); imported by aggregate_ucr.py
            (clade_variability, load_taxonomy)
Inputs:     --prot_ali         FASTA protein alignment; the gene name is basename.replace('.fa', '')
            --taxid_tsv        taxonomy table: species in column 2, family in column 3, name class
                               in column 5 (only "scientific name" rows are used)
            --family_order_tsv optional family → order table (otherwise the order is "Unknown")
Outputs:    <out_dir>/<var_subdir>/<gene>.entropy.tsv        per-position t, r, g, C_trident, variability
            <out_dir>/<var_subdir>/<gene>.fa                 columns with variability > 0 only
            <out_dir>/<var_subdir>/<gene>.clade_entropy.tsv  mean and max variability per order and
                                                             family (clades with at least 3 sequences
                                                             whose FASTA header matches a species)
            stdout: gene, n_seqs, n_cols, n_variable, mean_var, max_var (tab-separated)
"""

import argparse
import math
import os
import sys
from collections import defaultdict

# ── BLOSUM62 coordinates ──────────────────────────────────────────────────────

# 20x20 symmetric matrix, rows and columns in BLOSUM_AAS order (source: NCBI BLAST blosum62).
BLOSUM_AAS = list('ARNDCQEGHILKMFPSTWYV')
_BLOSUM62_RAW = [
# A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
[ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
[-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
[-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
[-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
[ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
[-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
[-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
[ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
[-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
[-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
[-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
[-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],  # K
[-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
[-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
[-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
[ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
[ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
[-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
[-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
[ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
]

AA_IDX = {aa: i for i, aa in enumerate(BLOSUM_AAS)}

import statistics as _stats

def _build_blosum_coords():
    """20D coordinate of each amino acid: its BLOSUM62 row with every column centered and divided by its std."""
    raw = _BLOSUM62_RAW
    n = len(BLOSUM_AAS)
    coords = {}
    col_means = [sum(raw[i][j] for i in range(n)) / n for j in range(n)]
    col_vars  = [sum((raw[i][j] - col_means[j])**2 for i in range(n)) / n for j in range(n)]
    col_stds  = [v**0.5 if v > 0 else 1.0 for v in col_vars]
    for i, aa in enumerate(BLOSUM_AAS):
        coords[aa] = tuple((raw[i][j] - col_means[j]) / col_stds[j] for j in range(n))
    return coords

BLOSUM_COORDS = _build_blosum_coords()

# Normalizer of r: mean Euclidean distance over all pairs of amino-acid coordinates.
def _max_mean_dist():
    aas = BLOSUM_AAS
    n = len(aas)
    total = 0.0
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((BLOSUM_COORDS[aas[i]][k] - BLOSUM_COORDS[aas[j]][k])**2
                              for k in range(n)))
            total += d
            cnt += 1
    return total / cnt if cnt else 1.0

_R_MAX = _max_mean_dist()

# ── I/O helpers ───────────────────────────────────────────────────────────────


def read_fasta(path):
    """Read a FASTA file into [(header, upper-cased sequence)]."""
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


def write_fasta(path, seqs):
    with open(path, 'w') as fh:
        for header, seq in seqs:
            fh.write(f'>{header}\n{seq}\n')


def load_taxonomy(taxid_tsv, family_order_tsv):
    """Map lower-cased species name → {'order', 'family'} from the taxid table.

    Rows with fewer than 5 columns or a name class other than "scientific name" are skipped,
    so a 2-column tax_id file gives an empty map. The order comes from family_order_tsv
    (family, order) and is 'Unknown' without it.
    """
    family_to_order = {}
    if family_order_tsv and os.path.isfile(family_order_tsv):
        with open(family_order_tsv) as fh:
            next(fh)
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    family_to_order[parts[0].strip()] = parts[1].strip()

    tax = {}
    with open(taxid_tsv) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            if parts[4].strip() != 'scientific name':
                continue
            species = parts[1].strip()
            family = parts[2].strip()
            order = family_to_order.get(family, 'Unknown')
            tax[species.lower()] = {'order': order, 'family': family}
    if not tax:
        print(f"WARNING: no usable row in {taxid_tsv} (rows need at least 5 columns, with the name class "
              "'scientific name' in column 5): the clade variability is empty.", file=sys.stderr)
    return tax

# ── Henikoff & Henikoff (1994) sequence weights ───────────────────────────────


def henikoff_weights(seqs_matrix):
    """Henikoff and Henikoff sequence weights, normalized to sum to 1 (uniform for an invariant alignment)."""
    n_seqs = len(seqs_matrix)
    if n_seqs == 0:
        return []
    n_cols = len(seqs_matrix[0])
    weights = [0.0] * n_seqs
    for col in range(n_cols):
        col_res = [seqs_matrix[i][col] for i in range(n_seqs)]
        distinct = set(col_res)
        r = len(distinct)
        if r <= 1:
            continue
        counts = defaultdict(int)
        for aa in col_res:
            counts[aa] += 1
        for i, aa in enumerate(col_res):
            weights[i] += 1.0 / (r * counts[aa])
    total = sum(weights)
    if total == 0:
        return [1.0 / n_seqs] * n_seqs
    return [w / total for w in weights]

# ── Three Valdar prongs ───────────────────────────────────────────────────────

_VALID_AAS = set(BLOSUM_AAS)
_GAP_CHARS = set('-X')


def symbol_diversity_t(col_residues, seq_weights):
    """Prong 1: weighted Shannon entropy over 20 amino acids + gap as 21st symbol (Eq. 50-52).

    Normalized by lambda_t = 1 / log2(min(N, 21)), N being the number of sequences. Returns t in [0, 1].
    Characters outside the 20 amino acids count as gap.
    """
    freq = defaultdict(float)
    for aa, w in zip(col_residues, seq_weights):
        symbol = aa if (aa in _VALID_AAS or aa in _GAP_CHARS) else '-'
        freq[symbol] += w

    # renormalize: the weights of a sub-alignment need not sum to 1
    total = sum(freq.values())
    if total == 0:
        return 0.0

    H = -sum((p / total) * math.log2(p / total) for p in freq.values() if p > 0)
    n_symbols = min(len(col_residues), 21)
    lambda_t = 1.0 / math.log2(n_symbols) if n_symbols > 1 else 1.0
    return min(H * lambda_t, 1.0)


def stereochem_diversity_r(col_residues, seq_weights):
    """Prong 2: mean weighted distance from the consensus point in BLOSUM62 20D space (Eq. 54-56).

    Only the 20 amino acids contribute; gaps are excluded from this prong. Returns r in [0, 1].
    """
    # total weight of each amino-acid type
    type_weight = defaultdict(float)
    for aa, w in zip(col_residues, seq_weights):
        if aa in _VALID_AAS:
            type_weight[aa] += w

    if not type_weight:
        return 0.0

    n_dim = len(BLOSUM_AAS)
    total_w = sum(type_weight.values())

    # consensus: weighted mean point
    consensus = [0.0] * n_dim
    for aa, w in type_weight.items():
        v = BLOSUM_COORDS[aa]
        frac = w / total_w
        for k in range(n_dim):
            consensus[k] += frac * v[k]

    # weighted mean Euclidean distance of the types from the consensus
    mean_dist = 0.0
    for aa, w in type_weight.items():
        v = BLOSUM_COORDS[aa]
        d = math.sqrt(sum((v[k] - consensus[k])**2 for k in range(n_dim)))
        mean_dist += (w / total_w) * d

    return min(mean_dist / _R_MAX, 1.0)


def gap_fraction_g(col_residues):
    """Prong 3: fraction of gap or unknown characters (anything but the 20 amino acids)."""
    n_gaps = sum(1 for aa in col_residues if aa in _GAP_CHARS or aa not in _VALID_AAS)
    return n_gaps / len(col_residues) if col_residues else 0.0


def valdar_column(col_residues, seq_weights, alpha=1.0, beta=1.0, gamma=1.0):
    """(t, r, g, C_trident, variability) of one column; variability = 1 - C_trident."""
    t = symbol_diversity_t(col_residues, seq_weights)
    r = stereochem_diversity_r(col_residues, seq_weights)
    g = gap_fraction_g(col_residues)
    C = ((1 - t) ** alpha) * ((1 - r) ** beta) * ((1 - g) ** gamma)
    return t, r, g, C, 1.0 - C

# ── Per-column and clade computation ──────────────────────────────────────────


def compute_per_column(seqs, seq_weights, alpha=1.0, beta=1.0, gamma=1.0):
    """Per-column (t, r, g, C, var) tuples and the 0-based indices of the columns with variability > 0."""
    n_cols = len(seqs[0][1]) if seqs else 0
    per_col = []
    variable_cols = []
    for col in range(n_cols):
        col_res = [s[1][col] for s in seqs]
        row = valdar_column(col_res, seq_weights, alpha, beta, gamma)
        per_col.append(row)
        if row[4] > 0.0:
            variable_cols.append(col)
    return per_col, variable_cols


def clade_variability(seqs, alpha=1.0, beta=1.0, gamma=1.0):
    """Mean and max variability over the variable columns of a sub-alignment.

    Sequence weights and all three prongs are recomputed on the sub-alignment alone.
    """
    seq_matrix = [s[1] for s in seqs]
    weights = henikoff_weights(seq_matrix)
    per_col, var_cols = compute_per_column(seqs, weights, alpha, beta, gamma)
    vals = [per_col[c][4] for c in var_cols]
    mean_v = sum(vals) / len(vals) if vals else 0.0
    max_v = max(vals) if vals else 0.0
    return mean_v, max_v

# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prot_ali', required=True)
    parser.add_argument('--taxid_tsv', required=True)
    parser.add_argument('--family_order_tsv', default=None)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--alpha', type=float, default=1.0,
                        help='Exponent for symbol-diversity prong t (default 1)')
    parser.add_argument('--beta',  type=float, default=1.0,
                        help='Exponent for stereochem-diversity prong r (default 1)')
    parser.add_argument('--gamma', type=float, default=1.0,
                        help='Exponent for gap-fraction prong g (default 1)')
    parser.add_argument('--var_subdir', default='PROT_VAR',
                        help='Subdirectory under --out_dir for variability outputs '
                             '(default: PROT_VAR; use PROT_VAR_RAW for raw branch)')
    args = parser.parse_args()

    gene = os.path.basename(args.prot_ali).replace('.fa', '')
    var_dir = os.path.join(args.out_dir, args.var_subdir)
    os.makedirs(var_dir, exist_ok=True)

    seqs = read_fasta(args.prot_ali)
    if not seqs:
        print(f'{gene}\t0\t0\t0\t0.000000\t0.000000', flush=True)
        return

    n_seqs = len(seqs)
    n_cols = len(seqs[0][1])

    seq_matrix = [s[1] for s in seqs]
    seq_weights = henikoff_weights(seq_matrix)
    per_col, variable_cols = compute_per_column(seqs, seq_weights,
                                                args.alpha, args.beta, args.gamma)

    # per-position table (positions are 1-based alignment columns)
    ent_path = os.path.join(var_dir, f'{gene}.entropy.tsv')
    with open(ent_path, 'w') as fh:
        fh.write('gene\tposition\tt\tr\tg\tC_trident\tvariability\tn_seqs\n')
        for col, (t, r, g, C, var) in enumerate(per_col):
            fh.write(f'{gene}\t{col + 1}\t{t:.6f}\t{r:.6f}\t{g:.6f}'
                     f'\t{C:.6f}\t{var:.6f}\t{n_seqs}\n')

    # alignment restricted to the variable columns
    var_seqs = [(h, ''.join(seq[c] for c in variable_cols)) for h, seq in seqs]
    write_fasta(os.path.join(var_dir, f'{gene}.fa'), var_seqs)

    # gene-level summary
    n_variable = len(variable_cols)
    vals = [per_col[c][4] for c in variable_cols]
    mean_var = sum(vals) / len(vals) if vals else 0.0
    max_var = max(vals) if vals else 0.0

    # sequence headers are matched to species names with spaces written as underscores
    taxonomy = load_taxonomy(args.taxid_tsv, args.family_order_tsv)

    def norm_key(name):
        return name.replace(' ', '_').lower()

    header_to_tax = {}
    for header, _ in seqs:
        info = taxonomy.get(norm_key(header))
        if info:
            header_to_tax[header] = info

    # per-clade variability at order and family level
    clade_rows = []
    for level in ('order', 'family'):
        clade_to_indices = defaultdict(list)
        for i, (header, _) in enumerate(seqs):
            info = header_to_tax.get(header)
            if info:
                clade_to_indices[info[level]].append(i)
        for clade_name, indices in sorted(clade_to_indices.items()):
            if len(indices) < 3:
                continue
            sub_seqs = [seqs[i] for i in indices]
            c_mean, c_max = clade_variability(sub_seqs, args.alpha, args.beta, args.gamma)
            clade_rows.append((gene, level, clade_name, len(indices),
                               f'{c_mean:.6f}', f'{c_max:.6f}'))

    clade_path = os.path.join(var_dir, f'{gene}.clade_entropy.tsv')
    with open(clade_path, 'w') as fh:
        fh.write('gene\tclade_level\tclade_name\tn_seqs\tmean_variability\tmax_variability\n')
        for row in clade_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')

    print(f'{gene}\t{n_seqs}\t{n_cols}\t{n_variable}\t{mean_var:.6f}\t{max_var:.6f}', flush=True)


if __name__ == '__main__':
    main()
