#!/usr/bin/env python3
"""
Detect Ultra-Conserved Regions (UCRs) in a protein alignment from per-position
Valdar C_trident scores (computed by compute_variability.py).

Three complementary detection methods are applied:

  absolute  -- consecutive positions with C_trident >= threshold and g <= max_gap.
               Answers: "what is absolutely conserved across the full phylogeny?"

  relative  -- consecutive positions whose within-gene C_trident z-score >= threshold.
               Answers: "what is unusually conserved for this specific gene?"
               Uniformly conserved genes yield no hits (no within-gene contrast).

  sliding   -- sliding window of k positions where mean C_trident >= threshold and
               every position has g <= max_pos_gap; overlapping windows are merged.
               Answers: "are there structurally interesting stretches with locally
               elevated conservation, tolerant of isolated gappy positions?"

Each detected region gains ±flank_size flanking positions (clamped to alignment
boundaries), stored as flank_start / flank_end.

Output: <gene>.ucr.tsv  (file is not written when no regions are detected, so that
Nextflow's `optional: true` can suppress empty per-gene jobs downstream)

Stdout: one summary line per gene  gene<TAB>n_absolute<TAB>n_relative<TAB>n_sliding
"""

import argparse
import math
import os
import sys


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_entropy_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip().split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                continue
            d = dict(zip(header, parts))
            rows.append({
                'position':    int(d['position']),
                'C_trident':   float(d['C_trident']),
                'variability': float(d['variability']),
                'g':           float(d['g']),
            })
    return rows


# ---------------------------------------------------------------------------
# Record helpers
# ---------------------------------------------------------------------------

def _stats(positions):
    n        = len(positions)
    mean_C   = sum(p['C_trident']   for p in positions) / n
    min_C    = min(p['C_trident']   for p in positions)
    mean_var = sum(p['variability'] for p in positions) / n
    mean_g   = sum(p['g']           for p in positions) / n
    return mean_C, min_C, mean_var, mean_g


def _make(gene, method, counter, positions, n_total, flank_size):
    start_pos  = positions[0]['position']
    end_pos    = positions[-1]['position']
    mean_C, min_C, mean_var, mean_g = _stats(positions)
    return {
        'gene':             gene,
        'ucr_id':           f'{gene}_{method}_{counter}',
        'method':           method,
        'start_pos':        start_pos,
        'end_pos':          end_pos,
        'flank_start':      max(1, start_pos - flank_size),
        'flank_end':        min(n_total, end_pos + flank_size),
        'n_positions':      len(positions),
        'mean_C_trident':   round(mean_C,   6),
        'min_C_trident':    round(min_C,    6),
        'mean_variability': round(mean_var, 6),
        'mean_gap':         round(mean_g,   6),
    }


# ---------------------------------------------------------------------------
# Detection methods
# ---------------------------------------------------------------------------

def detect_absolute(gene, rows, threshold, min_len, max_gap, flank_size):
    """Consecutive positions with C_trident >= threshold and g <= max_gap."""
    ucrs, counter, run = [], 1, []
    for p in rows:
        if p['C_trident'] >= threshold and p['g'] <= max_gap:
            run.append(p)
        else:
            if len(run) >= min_len:
                ucrs.append(_make(gene, 'absolute', counter, run, len(rows), flank_size))
                counter += 1
            run = []
    if len(run) >= min_len:
        ucrs.append(_make(gene, 'absolute', counter, run, len(rows), flank_size))
    return ucrs


def detect_relative(gene, rows, zscore_thresh, min_len, flank_size):
    """
    Consecutive positions with C_trident z-score >= zscore_thresh (within-gene).
    Returns empty list for uniformly conserved genes (std < 1e-10).
    """
    if len(rows) < 3:
        return []
    vals  = [p['C_trident'] for p in rows]
    mean_c = sum(vals) / len(vals)
    std_c  = math.sqrt(sum((v - mean_c) ** 2 for v in vals) / len(vals))
    if std_c < 1e-10:
        return []
    ucrs, counter, run = [], 1, []
    for p in rows:
        if (p['C_trident'] - mean_c) / std_c >= zscore_thresh:
            run.append(p)
        else:
            if len(run) >= min_len:
                ucrs.append(_make(gene, 'relative', counter, run, len(rows), flank_size))
                counter += 1
            run = []
    if len(run) >= min_len:
        ucrs.append(_make(gene, 'relative', counter, run, len(rows), flank_size))
    return ucrs


def detect_sliding(gene, rows, window_size, threshold, max_pos_gap, flank_size):
    """
    Sliding window: qualify windows with mean C_trident >= threshold and all
    positions having g <= max_pos_gap.  Overlapping windows are merged into
    contiguous spans.
    """
    n = len(rows)
    if n < window_size:
        return []
    covered = [False] * n
    for i in range(n - window_size + 1):
        win = rows[i:i + window_size]
        if all(p['g'] <= max_pos_gap for p in win):
            if sum(p['C_trident'] for p in win) / window_size >= threshold:
                for j in range(i, i + window_size):
                    covered[j] = True
    ucrs, counter, run = [], 1, []
    for i, ok in enumerate(covered):
        if ok:
            run.append(rows[i])
        else:
            if run:
                ucrs.append(_make(gene, 'sliding', counter, run, n, flank_size))
                counter += 1
                run = []
    if run:
        ucrs.append(_make(gene, 'sliding', counter, run, n, flank_size))
    return ucrs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

COLUMNS = [
    'gene', 'ucr_id', 'method',
    'start_pos', 'end_pos', 'flank_start', 'flank_end',
    'n_positions', 'mean_C_trident', 'min_C_trident', 'mean_variability', 'mean_gap',
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--entropy_tsv',          required=True,
                    help='<gene>.entropy.tsv from compute_variability.py')
    ap.add_argument('--gene',                 required=True)
    ap.add_argument('--out_tsv',              required=True)
    # absolute
    ap.add_argument('--abs_threshold',        type=float, default=0.9,
                    help='C_trident floor for absolute UCR detection (default: 0.9)')
    ap.add_argument('--abs_min_len',          type=int,   default=5,
                    help='Minimum run length in absolute method (default: 5)')
    ap.add_argument('--abs_max_gap',          type=float, default=0.3,
                    help='Max gap fraction for absolute method (default: 0.3)')
    # relative
    ap.add_argument('--rel_zscore',           type=float, default=1.5,
                    help='C_trident within-gene z-score threshold (default: 1.5)')
    ap.add_argument('--rel_min_len',          type=int,   default=5,
                    help='Minimum run length in relative method (default: 5)')
    # sliding window
    ap.add_argument('--window_size',          type=int,   default=10,
                    help='Sliding window size in aa (default: 10)')
    ap.add_argument('--window_threshold',     type=float, default=0.85,
                    help='Mean C_trident threshold for sliding window (default: 0.85)')
    ap.add_argument('--window_max_pos_gap',   type=float, default=0.3,
                    help='Max per-position gap fraction in sliding window (default: 0.3)')
    # shared
    ap.add_argument('--flank_size',           type=int,   default=5,
                    help='Flanking aa on each side of core window (default: 5)')
    args = ap.parse_args()

    rows = load_entropy_tsv(args.entropy_tsv)
    if not rows:
        print(f'{args.gene}\t0\t0\t0', flush=True)
        return

    ucrs = (
        detect_absolute(args.gene, rows, args.abs_threshold, args.abs_min_len,
                        args.abs_max_gap, args.flank_size)
        + detect_relative(args.gene, rows, args.rel_zscore, args.rel_min_len,
                          args.flank_size)
        + detect_sliding(args.gene, rows, args.window_size, args.window_threshold,
                         args.window_max_pos_gap, args.flank_size)
    )

    n_abs = sum(1 for u in ucrs if u['method'] == 'absolute')
    n_rel = sum(1 for u in ucrs if u['method'] == 'relative')
    n_sli = sum(1 for u in ucrs if u['method'] == 'sliding')
    print(f'{args.gene}\t{n_abs}\t{n_rel}\t{n_sli}', flush=True)

    if not ucrs:
        return  # no file → Nextflow optional: true handles absence

    with open(args.out_tsv, 'w') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for ucr in ucrs:
            fh.write('\t'.join(str(ucr[k]) for k in COLUMNS) + '\n')


if __name__ == '__main__':
    main()
