#!/usr/bin/env python3
"""
Aggregate UCR windows (from detect_ucr.py) across all genes and produce:

  ucr_windows.tsv            One row per detected UCR (all genes, all methods).
  ucr_positions.tsv          Per-position table within each UCR's full window
                             (core + flanks), enriched with optional Pfam domain
                             overlap and positive-selection calls.
  ucr_clade_variability.tsv  Per-UCR × clade mean variability for core and flanks
                             (requires PROT_RAW FASTA files + taxid_tsv).
  ucr_selection_summary.tsv  Per-UCR selection count summary (written only when
                             --fubar_sites or --fel_results is provided).

Position cross-referencing for positive selection:
  HyPhy site numbers are 1-based positions in the BMGE-trimmed protein alignment,
  while entropy.tsv positions are 1-based raw protein alignment columns.  When
  --map_dir is provided (*.map.tsv files from parse_bmge_track.py), the mapping
  ori_codon_col → trim_codon_col is used to translate raw positions to HyPhy site
  numbers.  Without MAP files, raw and trim positions are assumed equal (accurate
  for genes where BMGE removed very few columns).

Usage:
  python aggregate_ucr.py \\
      --ucr_dir     UCR_RAW         \\
      --entropy_dir PROT_VAR_RAW    \\
      --prot_dir    PROT_RAW        \\
      --taxid_tsv   taxid_...tsv    \\
      --out_dir     .               \\
    [ --flank_size        5 ]       \\
    [ --family_order_tsv  ... ]     \\
    [ --map_dir           MAP ]     \\
    [ --domain_tsv        domain_variability.tsv ]  \\
    [ --fubar_sites       fubar_sites.tsv ]         \\
    [ --fel_results       fel_results.tsv ]         \\
    [ --alpha 1.0 --beta 1.0 --gamma 1.0 ]
"""

import argparse
import glob
import math
import os
import sys
from collections import defaultdict

# Import Valdar helpers from sibling script (co-located in subworkflows/variability/local/).
# __file__ resolves to the absolute projectDir path when called via Nextflow.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_variability import clade_variability, load_taxonomy


# ---------------------------------------------------------------------------
# Sentinels
# ---------------------------------------------------------------------------

def _is_sentinel(path):
    """Return True for Nextflow sentinel files like NO_FUBAR, NO_DOMAIN_TSV, …"""
    return path is None or os.path.basename(path).startswith('NO_')


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_fasta(path):
    seqs, hdr, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if hdr is not None:
                    seqs.append((hdr, ''.join(buf)))
                hdr, buf = line[1:], []
            else:
                buf.append(line.upper())
    if hdr is not None:
        seqs.append((hdr, ''.join(buf)))
    return seqs


def load_ucr_tsvs(ucr_dir):
    """Returns dict  gene -> list[ucr_record_dict]."""
    result = {}
    for path in sorted(glob.glob(os.path.join(ucr_dir, '*.ucr.tsv'))):
        gene = os.path.basename(path).replace('.ucr.tsv', '')
        rows = []
        with open(path) as fh:
            header = fh.readline().rstrip().split('\t')
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < len(header):
                    continue
                d = dict(zip(header, parts))
                rows.append({
                    'gene':             d['gene'],
                    'ucr_id':           d['ucr_id'],
                    'method':           d['method'],
                    'start_pos':        int(d['start_pos']),
                    'end_pos':          int(d['end_pos']),
                    'flank_start':      int(d['flank_start']),
                    'flank_end':        int(d['flank_end']),
                    'n_positions':      int(d['n_positions']),
                    'mean_C_trident':   float(d['mean_C_trident']),
                    'min_C_trident':    float(d['min_C_trident']),
                    'mean_variability': float(d['mean_variability']),
                    'mean_gap':         float(d['mean_gap']),
                })
        if rows:
            result[gene] = rows
    return result


def load_entropy_index(entropy_dir):
    """Returns dict  gene -> {position(int): {C_trident, variability, g}}."""
    result = {}
    for path in sorted(glob.glob(os.path.join(entropy_dir, '*.entropy.tsv'))):
        gene = os.path.basename(path).replace('.entropy.tsv', '')
        pos_map = {}
        with open(path) as fh:
            header = fh.readline().rstrip().split('\t')
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < len(header):
                    continue
                d = dict(zip(header, parts))
                pos = int(d['position'])
                pos_map[pos] = {
                    'C_trident':   float(d['C_trident']),
                    'variability': float(d['variability']),
                    'g':           float(d['g']),
                }
        result[gene] = pos_map
    return result


def load_map_files(map_dir):
    """
    Returns dict  gene -> {ori_codon_col(int): trim_codon_col(int)}.
    Skips sentinel files.
    """
    result = {}
    for path in sorted(glob.glob(os.path.join(map_dir, '*.map.tsv'))):
        if _is_sentinel(path):
            continue
        gene = os.path.basename(path).replace('.map.tsv', '')
        ori_to_trim = {}
        with open(path) as fh:
            fh.readline()  # header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                ori_col, status, trim_col = parts[0], parts[1], parts[2]
                if status == 'selected' and trim_col not in ('NA', ''):
                    try:
                        ori_to_trim[int(ori_col)] = int(trim_col)
                    except ValueError:
                        pass
        result[gene] = ori_to_trim
    return result


def load_fubar_sites(path):
    """Returns dict  (gene, site_int) -> {fubar_prob_pos, fubar_prob_neg, fubar_is_pos, fubar_is_neg}."""
    result = {}
    with open(path) as fh:
        header = fh.readline().rstrip().split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                continue
            d = dict(zip(header, parts))
            try:
                key = (d['gene'], int(d['site']))
                result[key] = {
                    'fubar_prob_pos': _safe_float(d.get('prob_pos')),
                    'fubar_prob_neg': _safe_float(d.get('prob_neg')),
                    'fubar_is_pos':   d.get('is_pos_hit', '0') == '1',
                    'fubar_is_neg':   d.get('is_neg_hit', '0') == '1',
                }
            except (ValueError, KeyError):
                pass
    return result


def load_fel_results(path):
    """Returns dict  (gene, site_int) -> {dN, dS, p_fel}."""
    result = {}
    with open(path) as fh:
        header = fh.readline().rstrip().split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                continue
            d = dict(zip(header, parts))
            try:
                key = (d['gene'], int(d['site']))
                result[key] = {
                    'dN':    _safe_float(d.get('dN')),
                    'dS':    _safe_float(d.get('dS')),
                    'p_fel': _safe_float(d.get('p_fel')),
                }
            except (ValueError, KeyError):
                pass
    return result


def load_domain_tsv(path):
    """Returns dict  gene -> list[{pfam_id, domain_instance, ali_start, ali_end}]."""
    result = {}
    with open(path) as fh:
        header = fh.readline().rstrip().split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                continue
            d = dict(zip(header, parts))
            gene = d.get('gene', '')
            if not gene:
                continue
            try:
                hit = {
                    'pfam_id':         d['pfam_id'],
                    'domain_instance': int(d['domain_instance']),
                    'ali_start':       int(d['ali_start']),
                    'ali_end':         int(d['ali_end']),
                }
                result.setdefault(gene, []).append(hit)
            except (ValueError, KeyError):
                pass
    return result


def _safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def region_type(pos, ucr):
    if pos < ucr['start_pos']:
        return 'flank_up'
    if pos > ucr['end_pos']:
        return 'flank_down'
    return 'core'


def find_domain(domain_list, pos):
    """Return pfam_id of the first overlapping domain hit, or None."""
    for hit in domain_list:
        if hit['ali_start'] <= pos <= hit['ali_end']:
            return hit['pfam_id']
    return None


def _fmt(v, digits=6):
    if isinstance(v, float) and math.isnan(v):
        return 'NA'
    if isinstance(v, float):
        return f'{v:.{digits}f}'
    return str(v)


# ---------------------------------------------------------------------------
# Clade variability within UCR windows
# ---------------------------------------------------------------------------

def _norm_key(name):
    return name.replace(' ', '_').lower()


def compute_ucr_clade_variability(gene, ucrs, prot_dir, taxonomy, alpha, beta, gamma):
    """
    For each UCR window in the gene, compute per-clade mean variability
    separately for the core positions and the flank positions.
    Returns a list of row dicts.
    """
    fa_path = os.path.join(prot_dir, f'{gene}.fa')
    if not os.path.isfile(fa_path):
        return []
    seqs = read_fasta(fa_path)
    if not seqs:
        return []

    header_to_tax = {}
    for header, _ in seqs:
        info = taxonomy.get(_norm_key(header))
        if info:
            header_to_tax[header] = info

    rows = []
    for ucr in ucrs:
        core_cols  = list(range(ucr['start_pos'] - 1, ucr['end_pos']))       # 0-based
        flank_cols = (
            list(range(ucr['flank_start'] - 1, ucr['start_pos'] - 1))
            + list(range(ucr['end_pos'],          ucr['flank_end']))
        )

        for clade_level in ('order', 'family'):
            clade_to_idx = defaultdict(list)
            for i, (header, _) in enumerate(seqs):
                info = header_to_tax.get(header)
                if info:
                    clade_to_idx[info[clade_level]].append(i)

            for clade_name, indices in sorted(clade_to_idx.items()):
                if len(indices) < 3:
                    continue
                clade_seqs = [seqs[i] for i in indices]

                if core_cols:
                    c_sub = [(h, ''.join(s[c] for c in core_cols if c < len(s)))
                             for h, s in clade_seqs]
                    mean_var_core, _ = clade_variability(c_sub, alpha, beta, gamma)
                else:
                    mean_var_core = float('nan')

                if flank_cols:
                    f_sub = [(h, ''.join(s[c] for c in flank_cols if c < len(s)))
                             for h, s in clade_seqs]
                    mean_var_flank, _ = clade_variability(f_sub, alpha, beta, gamma)
                else:
                    mean_var_flank = float('nan')

                rows.append({
                    'gene':           gene,
                    'ucr_id':         ucr['ucr_id'],
                    'method':         ucr['method'],
                    'clade_level':    clade_level,
                    'clade_name':     clade_name,
                    'n_seqs':         len(indices),
                    'mean_var_core':  mean_var_core,
                    'mean_var_flank': mean_var_flank,
                })
    return rows


# ---------------------------------------------------------------------------
# Selection class helper
# ---------------------------------------------------------------------------

def sel_class(fubar_row, fel_row):
    """Return 'positive', 'negative', or 'neutral'."""
    is_pos = fubar_row.get('fubar_is_pos', False)
    is_neg = fubar_row.get('fubar_is_neg', False)
    if not is_pos and fel_row:
        p = fel_row.get('p_fel', 1.0)
        if p < 0.05:
            is_pos = fel_row.get('dN', 0) > fel_row.get('dS', 0)
            is_neg = not is_pos
    if not is_neg and fel_row:
        p = fel_row.get('p_fel', 1.0)
        if p < 0.05 and not is_pos:
            is_neg = True
    if is_pos:
        return 'positive'
    if is_neg:
        return 'negative'
    return 'neutral'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ucr_dir',          required=True,
                    help='Directory containing *.ucr.tsv files from detect_ucr.py')
    ap.add_argument('--entropy_dir',      required=True,
                    help='Directory containing *.entropy.tsv files (PROT_VAR_RAW)')
    ap.add_argument('--prot_dir',         required=True,
                    help='Directory containing *.fa raw protein MSA files (PROT_RAW)')
    ap.add_argument('--taxid_tsv',        required=True)
    ap.add_argument('--out_dir',          required=True)
    ap.add_argument('--flank_size',       type=int,   default=5)
    ap.add_argument('--family_order_tsv', default=None)
    ap.add_argument('--map_dir',          default=None,
                    help='Directory of *.map.tsv files (MAP/) for raw↔trim position mapping')
    ap.add_argument('--domain_tsv',       default=None,
                    help='domain_variability.tsv from map_domain_variability.py')
    ap.add_argument('--fubar_sites',      default=None,
                    help='fubar_sites.tsv from aggregate_positive_selection.py')
    ap.add_argument('--fel_results',      default=None,
                    help='fel_results.tsv from aggregate_positive_selection.py')
    ap.add_argument('--alpha',  type=float, default=1.0)
    ap.add_argument('--beta',   type=float, default=1.0)
    ap.add_argument('--gamma',  type=float, default=1.0)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load core data ─────────────────────────────────────────────────────
    all_ucrs = load_ucr_tsvs(args.ucr_dir)
    if not all_ucrs:
        print('No UCR TSV files found — writing empty ucr_windows.tsv.', file=sys.stderr)
        open(os.path.join(args.out_dir, 'ucr_windows.tsv'), 'w').write(
            'gene\tucr_id\tmethod\tstart_pos\tend_pos\tflank_start\tflank_end\t'
            'n_positions\tmean_C_trident\tmin_C_trident\tmean_variability\tmean_gap\n'
        )
        open(os.path.join(args.out_dir, 'ucr_positions.tsv'), 'w').write(
            'gene\tucr_id\tmethod\tposition\tregion_type\tC_trident\tvariability\tg\n'
        )
        return

    entropy_idx = load_entropy_index(args.entropy_dir)

    has_maps   = args.map_dir and not _is_sentinel(args.map_dir)
    gene_maps  = load_map_files(args.map_dir) if has_maps else {}

    has_domain  = args.domain_tsv  and not _is_sentinel(args.domain_tsv)
    has_fubar   = args.fubar_sites and not _is_sentinel(args.fubar_sites)
    has_fel     = args.fel_results and not _is_sentinel(args.fel_results)

    domain_hits = load_domain_tsv(args.domain_tsv)  if has_domain else {}
    fubar_data  = load_fubar_sites(args.fubar_sites) if has_fubar  else {}
    fel_data    = load_fel_results(args.fel_results) if has_fel    else {}

    taxonomy = load_taxonomy(args.taxid_tsv, args.family_order_tsv)

    # ── ucr_windows.tsv ────────────────────────────────────────────────────
    win_cols = [
        'gene', 'ucr_id', 'method',
        'start_pos', 'end_pos', 'flank_start', 'flank_end',
        'n_positions', 'mean_C_trident', 'min_C_trident', 'mean_variability', 'mean_gap',
    ]
    win_path = os.path.join(args.out_dir, 'ucr_windows.tsv')
    n_win = 0
    with open(win_path, 'w') as fh:
        fh.write('\t'.join(win_cols) + '\n')
        for gene in sorted(all_ucrs):
            for ucr in all_ucrs[gene]:
                fh.write('\t'.join(str(ucr[c]) for c in win_cols) + '\n')
                n_win += 1
    print(f'Written: {win_path} ({n_win} UCR windows)', flush=True)

    # ── ucr_positions.tsv ─────────────────────────────────────────────────
    pos_base = ['gene', 'ucr_id', 'method', 'position', 'region_type',
                'C_trident', 'variability', 'g']
    pos_opt_dom  = ['pfam_domain']  if has_domain else []
    pos_opt_fub  = ['fubar_prob_pos', 'fubar_prob_neg', 'fubar_is_pos', 'fubar_is_neg'] \
                   if has_fubar else []
    pos_opt_fel  = ['dN', 'dS', 'p_fel'] if has_fel else []
    pos_opt_sel  = ['sel_class'] if (has_fubar or has_fel) else []
    pos_cols = pos_base + pos_opt_dom + pos_opt_fub + pos_opt_fel + pos_opt_sel

    pos_path = os.path.join(args.out_dir, 'ucr_positions.tsv')
    n_pos = 0
    with open(pos_path, 'w') as fh:
        fh.write('\t'.join(pos_cols) + '\n')
        for gene in sorted(all_ucrs):
            ent  = entropy_idx.get(gene, {})
            gmap = gene_maps.get(gene, {})
            doms = domain_hits.get(gene, [])
            for ucr in all_ucrs[gene]:
                for pos in range(ucr['flank_start'], ucr['flank_end'] + 1):
                    ent_row = ent.get(pos)
                    if ent_row is None:
                        continue
                    hphy = gmap.get(pos, pos)  # fallback: raw pos == trim pos
                    fubar_row = fubar_data.get((gene, hphy), {})
                    fel_row   = fel_data.get((gene, hphy),   {})

                    row = {
                        'gene':        gene,
                        'ucr_id':      ucr['ucr_id'],
                        'method':      ucr['method'],
                        'position':    pos,
                        'region_type': region_type(pos, ucr),
                        'C_trident':   _fmt(ent_row['C_trident']),
                        'variability': _fmt(ent_row['variability']),
                        'g':           _fmt(ent_row['g']),
                    }
                    if has_domain:
                        row['pfam_domain'] = find_domain(doms, pos) or 'NA'
                    if has_fubar:
                        row['fubar_prob_pos'] = _fmt(fubar_row.get('fubar_prob_pos', float('nan')), 4)
                        row['fubar_prob_neg'] = _fmt(fubar_row.get('fubar_prob_neg', float('nan')), 4)
                        row['fubar_is_pos']   = int(fubar_row.get('fubar_is_pos', False))
                        row['fubar_is_neg']   = int(fubar_row.get('fubar_is_neg', False))
                    if has_fel:
                        row['dN']    = _fmt(fel_row.get('dN',    float('nan')), 4)
                        row['dS']    = _fmt(fel_row.get('dS',    float('nan')), 4)
                        row['p_fel'] = _fmt(fel_row.get('p_fel', float('nan')), 6)
                    if has_fubar or has_fel:
                        row['sel_class'] = sel_class(fubar_row, fel_row)

                    fh.write('\t'.join(str(row[c]) for c in pos_cols) + '\n')
                    n_pos += 1
    print(f'Written: {pos_path} ({n_pos} position rows)', flush=True)

    # ── ucr_clade_variability.tsv ─────────────────────────────────────────
    clad_cols = ['gene', 'ucr_id', 'method', 'clade_level', 'clade_name',
                 'n_seqs', 'mean_var_core', 'mean_var_flank']
    clad_path = os.path.join(args.out_dir, 'ucr_clade_variability.tsv')
    n_clad = 0
    with open(clad_path, 'w') as fh:
        fh.write('\t'.join(clad_cols) + '\n')
        for gene in sorted(all_ucrs):
            for row in compute_ucr_clade_variability(
                    gene, all_ucrs[gene], args.prot_dir, taxonomy,
                    args.alpha, args.beta, args.gamma):
                fh.write('\t'.join([
                    row['gene'], row['ucr_id'], row['method'],
                    row['clade_level'], row['clade_name'],
                    str(row['n_seqs']),
                    _fmt(row['mean_var_core']),
                    _fmt(row['mean_var_flank']),
                ]) + '\n')
                n_clad += 1
    print(f'Written: {clad_path} ({n_clad} clade×UCR rows)', flush=True)

    # ── ucr_selection_summary.tsv (only when psel data is present) ─────────
    if not (has_fubar or has_fel):
        return

    sel_cols = [
        'gene', 'ucr_id', 'method',
        'n_positions_core', 'n_positions_flank',
        'n_fubar_pos_core',  'n_fubar_neg_core',
        'n_fubar_pos_flank', 'n_fubar_neg_flank',
        'pct_core_pos_sel',  'pct_core_neg_sel',
        'n_fel_pos_core',    'n_fel_neg_core',
    ]
    sel_path = os.path.join(args.out_dir, 'ucr_selection_summary.tsv')
    n_sel = 0
    with open(sel_path, 'w') as fh:
        fh.write('\t'.join(sel_cols) + '\n')
        for gene in sorted(all_ucrs):
            ent  = entropy_idx.get(gene, {})
            gmap = gene_maps.get(gene, {})
            for ucr in all_ucrs[gene]:
                n_core = n_flank = 0
                n_fub_pc = n_fub_nc = n_fub_pf = n_fub_nf = 0
                n_fel_pc = n_fel_nc = 0
                for pos in range(ucr['flank_start'], ucr['flank_end'] + 1):
                    if pos not in ent:
                        continue
                    hphy    = gmap.get(pos, pos)
                    is_core = ucr['start_pos'] <= pos <= ucr['end_pos']
                    if is_core:
                        n_core  += 1
                    else:
                        n_flank += 1

                    if has_fubar:
                        fr = fubar_data.get((gene, hphy), {})
                        if is_core:
                            n_fub_pc += int(fr.get('fubar_is_pos', False))
                            n_fub_nc += int(fr.get('fubar_is_neg', False))
                        else:
                            n_fub_pf += int(fr.get('fubar_is_pos', False))
                            n_fub_nf += int(fr.get('fubar_is_neg', False))

                    if has_fel and is_core:
                        fe = fel_data.get((gene, hphy), {})
                        p  = fe.get('p_fel', 1.0) if not math.isnan(fe.get('p_fel', float('nan'))) else 1.0
                        if p < 0.05:
                            if fe.get('dN', 0) > fe.get('dS', 0):
                                n_fel_pc += 1
                            else:
                                n_fel_nc += 1

                pct_pos = round(100.0 * n_fub_pc / n_core, 2) if n_core else 0.0
                pct_neg = round(100.0 * n_fub_nc / n_core, 2) if n_core else 0.0
                row = {
                    'gene': gene, 'ucr_id': ucr['ucr_id'], 'method': ucr['method'],
                    'n_positions_core':  n_core,   'n_positions_flank': n_flank,
                    'n_fubar_pos_core':  n_fub_pc, 'n_fubar_neg_core':  n_fub_nc,
                    'n_fubar_pos_flank': n_fub_pf, 'n_fubar_neg_flank': n_fub_nf,
                    'pct_core_pos_sel':  pct_pos,  'pct_core_neg_sel':  pct_neg,
                    'n_fel_pos_core':    n_fel_pc, 'n_fel_neg_core':    n_fel_nc,
                }
                fh.write('\t'.join(str(row[c]) for c in sel_cols) + '\n')
                n_sel += 1
    print(f'Written: {sel_path} ({n_sel} UCR selection summary rows)', flush=True)


if __name__ == '__main__':
    main()
