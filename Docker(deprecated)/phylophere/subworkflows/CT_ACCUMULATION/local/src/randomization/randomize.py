#!/usr/bin/env python3
"""
Randomization / permutation significance test — adapted for CT_ACCUMULATION in PhyloPhere.

Changes vs. original to_integrate version:
  - _remap_caas_df(): normalises global_meta_caas.tsv column names to the internal schema
    (Gene→gene, Position→msa_pos, sig_both→is_significant, IsConserved→is_stable,
     Pattern→pattern_type, Tag→tag, CAAP_Group→iscaap derivation).
    Falls back transparently when the original column schema is present.
  - Z-score p-values removed: PValueZ_* columns and scipy.stats.norm dependency removed.
    Only empirical p-values (and their BH-FDR corrections) are retained.
  - Gene list export: after each category CSV is written, two gene-list .txt files are
    exported to gene_lists/ — one using FDR-corrected p, one using raw empirical p.
  - --fdr-threshold arg added (default 0.05).
"""

import argparse
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pyarrow as pa
import pyarrow.parquet as pq
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import random
import multiprocessing.shared_memory as shm
from statsmodels.stats.multitest import multipletests


# ---------------------------
# Column remapping
# ---------------------------

def _remap_caas_df(df):
    """Normalise a CAAS DataFrame to the internal column schema.

    Handles two source formats:
      1. filtered_discovery.tsv / global_meta_caas.tsv produced by CT_POSTPROC / CT_SIGNIFICATION
         Columns include: Gene, Position, tag, is_significant, sig_both, Pattern,
         CAAP_Group, change_side, asr_is_conserved, is_conserved_meta, ...
      2. Original caas_csv used by the standalone tool:
         gene, msa_pos, is_stable, is_significant, change_side, isAntiCAAS,
         ASR_ConsAntiCAAS, tag, pattern_type, iscaap, caap_group
    """
    cols = set(df.columns)

    # Detect signification format: has Gene/Position but not yet msa_pos
    is_signification_fmt = (
        ('Gene' in cols or 'gene' in cols) and
        ('Position' in cols or 'position' in cols) and
        'msa_pos' not in cols
    )

    if not is_signification_fmt:
        logging.info("CAAS CSV: using original column schema")
        return df

    logging.info("CAAS CSV: detected filtered_discovery.tsv / signification format — remapping columns")

    # Rename Gene→gene, Position→msa_pos
    rename_map = {}
    if 'Gene' in cols:     rename_map['Gene'] = 'gene'
    if 'Position' in cols: rename_map['Position'] = 'msa_pos'
    df = df.rename(columns=rename_map)
    cols = set(df.columns)  # refresh after rename

    # is_significant: prefer sig_both (passes both hyp + perm tests) over raw is_significant.
    # Always coerce to Python bool regardless of source dtype.
    _bool = lambda x: str(x).strip().lower() in {'true', 't', '1', 'yes', 'y'}
    for cand in ['sig_both', 'is_significant', 'isSignificant']:
        if cand in cols:
            df['is_significant'] = df[cand].map(_bool)
            break
    else:
        df['is_significant'] = False

    # is_stable: prefer asr_is_conserved (filtered_discovery.tsv), then is_conserved_meta,
    # then IsConserved (older global_meta_caas.tsv); default True if absent.
    if 'is_stable' not in cols:
        for cand in ['asr_is_conserved', 'is_conserved_meta', 'IsConserved']:
            if cand in cols:
                df['is_stable'] = df[cand].map(_bool)
                break
        else:
            df['is_stable'] = True

    # tag — already lowercase in filtered_discovery.tsv; fall back to Tag
    if 'tag' not in cols:
        df['tag'] = df['Tag'] if 'Tag' in cols else ''

    # pattern_type: from Pattern column
    if 'pattern_type' not in cols:
        df['pattern_type'] = df['Pattern'] if 'Pattern' in cols else ''

    # caap_group: preserve the raw CAAP_Group label (GS0–GS4, US) — needed for
    # the GS0-exclusion logic in the convergent filter downstream.
    if 'caap_group' not in cols:
        if 'CAAP_Group' in cols:
            df['caap_group'] = df['CAAP_Group'].astype(str).str.strip().str.upper()
        else:
            df['caap_group'] = ''

    # iscaap: any labelled group other than 'US' (includes GS0–GS4)
    if 'iscaap' not in cols:
        df['iscaap'] = df['caap_group'].apply(lambda x: str(x).strip().upper() != 'US')

    # change_side: present in filtered_discovery.tsv; warn only if truly absent
    if 'change_side' not in cols:
        logging.warning(
            "change_side not found in CAAS CSV. Defaulting to 'top' for all entries. "
            "Top/bottom/both category lists will reflect this simplification."
        )
        df['change_side'] = 'top'

    # AntiCAAS flags — not propagated via signification pipeline; safe defaults
    if 'isAntiCAAS' not in cols:
        df['isAntiCAAS'] = False
    if 'ASR_ConsAntiCAAS' not in cols:
        df['ASR_ConsAntiCAAS'] = True

    return df


# ---------------------------
# Worker class
# ---------------------------

class RandomizationWorker:
    def __init__(
        self,
        randomization_type,
        decile_bins,
        export_individual,
        output_dir,
        global_seed,
        precompute_masks,
        n_rows,
        n_genes,
        shm_names,
        shapes,
        dtypes,
        extra_key_sizes,
        caas_data,
        position_to_tag,
        actual_counts_all,
        actual_sig,
        actual_top,
        actual_bottom,
        actual_both,
        actual_convergent,
        actual_convergent_caap,
        actual_convergent_caap_top,
        actual_convergent_caap_bottom,
        actual_convergent_caap_both,
    ):
        self.randomization_type = randomization_type
        self.decile_bins = np.array(decile_bins) if decile_bins is not None else None
        self.export_individual = export_individual
        self.output_dir = output_dir
        self.global_seed = global_seed
        self.precompute_masks = precompute_masks
        self.n_rows = n_rows
        self.n_genes = n_genes
        self.caas_data = caas_data
        self.position_to_tag = position_to_tag or {}
        self.actual_all    = actual_counts_all.astype(np.int64, copy=False)
        self.actual_sig    = actual_sig.astype(np.int64, copy=False)
        self.actual_top    = actual_top.astype(np.int64, copy=False)
        self.actual_bottom = actual_bottom.astype(np.int64, copy=False)
        self.actual_both   = actual_both.astype(np.int64, copy=False)
        self.actual_convergent         = actual_convergent.astype(np.int64, copy=False)
        self.actual_convergent_caap    = actual_convergent_caap.astype(np.int64, copy=False)
        self.actual_convergent_caap_top    = actual_convergent_caap_top.astype(np.int64, copy=False)
        self.actual_convergent_caap_bottom = actual_convergent_caap_bottom.astype(np.int64, copy=False)
        self.actual_convergent_caap_both   = actual_convergent_caap_both.astype(np.int64, copy=False)

        if self.global_seed is not None:
            np.random.seed(self.global_seed)
            random.seed(self.global_seed)

        self.caas_by_key_counts = defaultdict(lambda: {
            'n_all': 0, 'n_sig': 0,
            'n_top': 0, 'n_bottom': 0, 'n_both': 0,
            'n_convergent': 0, 'n_convergent_caap': 0,
            'n_convergent_caap_top': 0, 'n_convergent_caap_bottom': 0,
            'n_convergent_caap_both': 0,
        })
        for k, payload in extra_key_sizes.items():
            if k not in self.caas_by_key_counts:
                self.caas_by_key_counts[k]
            for name, val in payload.items():
                self.caas_by_key_counts[k][name] = int(val)

        if self.export_individual:
            byk = defaultdict(list)
            for rec in self.caas_data:
                byk[rec['key']].append(rec)
            self.caas_by_key_lists = dict(byk)

        self._open_shared(shm_names, shapes, dtypes)
        self._prepare_keyed_eligibles()

    def _open_shared(self, shm_names, shapes, dtypes):
        self._shm_positions = shm.SharedMemory(name=shm_names['positions'])
        self.positions = np.ndarray(shapes['positions'], dtype=dtypes['positions'], buffer=self._shm_positions.buf)
        self._shm_masked = shm.SharedMemory(name=shm_names['masked'])
        self.masked = np.ndarray(shapes['masked'], dtype=dtypes['masked'], buffer=self._shm_masked.buf)
        self._shm_cons = shm.SharedMemory(name=shm_names['cons_idx'])
        self.cons_idx = np.ndarray(shapes['cons_idx'], dtype=dtypes['cons_idx'], buffer=self._shm_cons.buf)
        self._shm_genes = shm.SharedMemory(name=shm_names['genes'])
        self.genes = np.ndarray(shapes['genes'], dtype=dtypes['genes'], buffer=self._shm_genes.buf)
        self._shm_caas = shm.SharedMemory(name=shm_names['iscaas'])
        self.caas = np.ndarray(shapes['iscaas'], dtype=dtypes['iscaas'], buffer=self._shm_caas.buf)
        self.row_indices = np.arange(self.n_rows, dtype=np.int32)

    def _prepare_keyed_eligibles(self):
        self.eligible_by_key = {}
        all_eligible = self.row_indices[~self.masked]
        if self.randomization_type == 'naive':
            self.eligible_by_key['global'] = all_eligible
            print(f"Total eligible positions for 'naive' randomization: {len(all_eligible)}")
            return
        if self.randomization_type == 'cons_decile':
            assert self.decile_bins is not None
            dec = np.digitize(self.cons_idx, bins=self.decile_bins[:-1], right=False)
            for d in range(len(self.decile_bins) - 1):
                mask = dec == d
                self.eligible_by_key[d] = self.row_indices[mask]
            return
        raise ValueError(f"Unknown randomization_type: {self.randomization_type}")

    def process_chunk(self, chunk_size, chunk_idx):
        logging.info(f"Starting chunk {chunk_idx} with {chunk_size} randomizations in PID {os.getpid()}")
        if self.global_seed is not None:
            np.random.seed(((self.global_seed or 0) + 9973 * chunk_idx) & 0x7FFFFFFF)

        # Vectorized accumulators per category.
        # Must be a list (not a tuple) so that S[0] += c in-place assignment works.
        def _zeros():
            return [
                np.zeros(self.n_genes, dtype=np.int64),
                np.zeros(self.n_genes, dtype=np.float64),
                np.zeros(self.n_genes, dtype=np.int64),
            ]
        S_all    = _zeros(); S_sig    = _zeros()
        S_top    = _zeros(); S_bottom = _zeros(); S_both = _zeros()
        S_conv   = _zeros(); S_conv_c = _zeros()
        S_cc_top = _zeros(); S_cc_bot = _zeros(); S_cc_both = _zeros()

        schema = None
        writer = None
        batch_rows = []
        if self.export_individual:
            schema = pa.schema([
                ('rand_id', pa.string()),
                ('caas_position', pa.int64()),
                ('tag', pa.string()),
                ('randomized_position', pa.int64()),
                ('randomized_gene_id', pa.int32()),
            ])
            os.makedirs(self.output_dir, exist_ok=True)
            chunk_file = os.path.join(self.output_dir, f"chunk_{chunk_idx}_{os.getpid()}.parquet")
            writer = pq.ParquetWriter(chunk_file, schema, compression='snappy')

        bincount = np.bincount

        for r in range(chunk_size):
            c_all = np.zeros(self.n_genes, dtype=np.int64)
            c_sig = np.zeros(self.n_genes, dtype=np.int64)
            c_top = np.zeros(self.n_genes, dtype=np.int64)
            c_bot = np.zeros(self.n_genes, dtype=np.int64)
            c_both= np.zeros(self.n_genes, dtype=np.int64)
            c_conv = np.zeros(self.n_genes, dtype=np.int64)
            c_cc   = np.zeros(self.n_genes, dtype=np.int64)
            c_cct  = np.zeros(self.n_genes, dtype=np.int64)
            c_ccb  = np.zeros(self.n_genes, dtype=np.int64)
            c_ccboth = np.zeros(self.n_genes, dtype=np.int64)
            uid = f"{os.getpid()}_{chunk_idx}_{r}" if self.export_individual else None

            for key, cinfo in self.caas_by_key_counts.items():
                elig = self.eligible_by_key.get(key)
                if elig is None or elig.size == 0:
                    continue
                cat_map = [
                    ('n_all', c_all), ('n_sig', c_sig),
                    ('n_top', c_top), ('n_bottom', c_bot), ('n_both', c_both),
                    ('n_convergent', c_conv), ('n_convergent_caap', c_cc),
                    ('n_convergent_caap_top', c_cct), ('n_convergent_caap_bottom', c_ccb),
                    ('n_convergent_caap_both', c_ccboth),
                ]
                for cname, carr in cat_map:
                    n = cinfo.get(cname, 0)
                    if n > 0:
                        idx = np.random.choice(elig, size=n, replace=True)
                        carr += bincount(self.genes[idx], minlength=self.n_genes)

                if self.export_individual:
                    caas_list = self.caas_by_key_lists.get(key, [])
                    if caas_list:
                        idx_ind = np.random.choice(elig, size=len(caas_list), replace=True)
                        pos_rand  = self.positions[idx_ind]
                        gene_rand = self.genes[idx_ind]
                        for rec, pr, gid in zip(caas_list, pos_rand, gene_rand):
                            tag = rec.get('tag', self.position_to_tag.get(rec['position'], ''))
                            batch_rows.append({
                                'rand_id': uid,
                                'caas_position': int(rec['position']),
                                'tag': tag,
                                'randomized_position': int(pr),
                                'randomized_gene_id': int(gid),
                            })
                    if len(batch_rows) >= 1000 and writer is not None:
                        dfb = pd.DataFrame(batch_rows)
                        writer.write_table(pa.Table.from_pandas(dfb, schema=schema))
                        batch_rows = []

            def _upd(S, c, act):
                S[0] += c
                S[1] += c.astype(np.float64) ** 2
                S[2] += (c >= act).astype(np.int64)
            _upd(S_all, c_all, self.actual_all)
            _upd(S_sig, c_sig, self.actual_sig)
            _upd(S_top, c_top, self.actual_top)
            _upd(S_bottom, c_bot, self.actual_bottom)
            _upd(S_both, c_both, self.actual_both)
            _upd(S_conv, c_conv, self.actual_convergent)
            _upd(S_conv_c, c_cc, self.actual_convergent_caap)
            _upd(S_cc_top, c_cct, self.actual_convergent_caap_top)
            _upd(S_cc_bot, c_ccb, self.actual_convergent_caap_bottom)
            _upd(S_cc_both, c_ccboth, self.actual_convergent_caap_both)

        if self.export_individual and batch_rows and writer is not None:
            writer.write_table(pa.Table.from_pandas(pd.DataFrame(batch_rows), schema=schema))
        if writer is not None:
            writer.close()

        def _pack(S):
            return {'sum': S[0], 'sum_sq': S[1], 'count_above': S[2]}
        return {
            'chunk_idx': chunk_idx, 'n_rands': chunk_size,
            'all':    _pack(S_all),    'sig':    _pack(S_sig),
            'top':    _pack(S_top),    'bottom': _pack(S_bottom), 'both': _pack(S_both),
            'convergent':          _pack(S_conv),
            'convergent_caap':     _pack(S_conv_c),
            'convergent_caap_top': _pack(S_cc_top),
            'convergent_caap_bottom': _pack(S_cc_bot),
            'convergent_caap_both':   _pack(S_cc_both),
        }


# ---------------------------
# Worker helpers (module-level for pickling)
# ---------------------------

def init_worker(
    randomization_type, decile_bins, export_individual, output_dir, global_seed,
    precompute_masks, n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
    caas_data, position_to_tag,
    actual_all, actual_sig, actual_top, actual_bottom, actual_both,
    actual_convergent, actual_convergent_caap,
    actual_convergent_caap_top, actual_convergent_caap_bottom, actual_convergent_caap_both,
):
    global worker
    worker = RandomizationWorker(
        randomization_type, decile_bins, export_individual, output_dir, global_seed,
        precompute_masks, n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
        caas_data, position_to_tag,
        actual_all, actual_sig, actual_top, actual_bottom, actual_both,
        actual_convergent, actual_convergent_caap,
        actual_convergent_caap_top, actual_convergent_caap_bottom, actual_convergent_caap_both,
    )

def process_wrapper(ch):
    size, idx = ch
    return worker.process_chunk(size, idx)


# ---------------------------
# Utility functions
# ---------------------------

def _compute_bins_from_series(series):
    return np.percentile(series.dropna(), np.arange(0, 101, 10))

def _build_caas_payload(merged_df, randomization_type, decile_bins):
    caas_data = []
    for _, row in merged_df.iterrows():
        if pd.isna(row['iscaas']) or not row['iscaas']:
            continue
        if randomization_type == 'naive':
            key = 'global'
        elif randomization_type == 'cons_decile':
            key = int(np.digitize([row['cons_idx']], bins=decile_bins[:-1], right=False)[0])
        else:
            continue
        caas_data.append({'position': row['position'], 'key': key})
    unique_keys = list({item['key'] for item in caas_data})
    return unique_keys, caas_data

def _write_gene_lists(df_cat, cat, gene_lists_dir, fdr_threshold):
    """Write two gene-list .txt files for a category: one FDR-corrected, one raw empirical."""
    os.makedirs(gene_lists_dir, exist_ok=True)
    for suffix, col in [
        ('fdr', f'PValueEmpirical_{cat}_FDR'),
        ('raw', f'PValueEmpirical_{cat}'),
    ]:
        sig = df_cat[
            (df_cat[f'ActualCount_{cat}'] > 0) &
            (df_cat[col] < fdr_threshold)
        ]['Gene']
        out_path = os.path.join(gene_lists_dir, f"{cat}_significant_{suffix}.txt")
        sig.to_csv(out_path, index=False, header=False)
        logging.info(f"Gene list ({cat}, {suffix}): {len(sig)} genes → {out_path}")


def _write_empty_outputs_and_exit(args):
    """Write empty-but-valid outputs when no rows are available for randomization."""
    categories = [
        'all', 'sig', 'top', 'bottom', 'both',
        'convergent', 'convergent_caap',
        'convergent_caap_top', 'convergent_caap_bottom', 'convergent_caap_both',
    ]

    base_dir = os.path.dirname(args.output_prefix) or '.'
    os.makedirs(base_dir, exist_ok=True)

    gene_lists_dir = os.path.join(base_dir, 'gene_lists')
    os.makedirs(gene_lists_dir, exist_ok=True)

    for cat in categories:
        df_cat = pd.DataFrame(columns=[
            'GeneID',
            'Gene',
            f'ActualCount_{cat}',
            f'RandsAbove_{cat}',
            f'RandsBelow_{cat}',
            f'PValueEmpirical_{cat}',
            f'MeanSim_{cat}',
            f'SD_{cat}',
            f'NumRands_{cat}',
            f'PValueEmpirical_{cat}_FDR',
        ])

        out_cat = f"{args.output_prefix}_{cat}_aggregated_results.csv"
        df_cat.to_csv(out_cat, index=False, compression='gzip' if args.compress else None)
        logging.info(f"Results ({cat}): {out_cat} [empty]")

        # Keep output contract: create empty gene list files for every category.
        for suffix in ['fdr', 'raw']:
            out_path = os.path.join(gene_lists_dir, f"{cat}_significant_{suffix}.txt")
            with open(out_path, 'w'):
                pass

    if args.randomization_type == 'cons_decile':
        pd.DataFrame(columns=['decile', 'NumPositions']).to_parquet(
            f"{args.output_prefix}_cons_decile_global.parquet", index=False
        )
        pd.DataFrame(columns=['gene']).to_parquet(
            f"{args.output_prefix}_cons_decile_per_gene.parquet", index=False
        )

    logging.info("Randomization analysis complete (empty input).")


# ---------------------------
# Main
# ---------------------------

def main(args):
    logging.info("Loading data")
    # global_csv now contains both positional data (cons_idx) and group data (masked, iscaas)
    global_df = pd.read_csv(args.global_csv)
    caas_raw  = pd.read_csv(args.caas_csv, sep=None, engine='python')

    # Normalise CAAS columns (handles global_meta_caas.tsv or original schema)
    caas_df = _remap_caas_df(caas_raw)

    # Keep only needed columns to save memory
    required_cols = ['gene', 'msa_pos', 'is_stable', 'is_significant',
                     'change_side', 'isAntiCAAS', 'ASR_ConsAntiCAAS',
                     'tag', 'pattern_type', 'iscaap', 'caap_group']
    available = [c for c in required_cols if c in caas_df.columns]
    caas_df = caas_df[available]
    logging.info(f"CAAS df: {len(caas_df)} rows, columns: {available}")

    # Dedup
    if global_df.duplicated('position').any():
        logging.warning('global_df has duplicate positions; keeping first occurrence.')
        global_df = global_df.sort_values('position').drop_duplicates('position', keep='first')

    merged_df = pd.merge(global_df, caas_df, how='left', on=['gene', 'msa_pos'])
    logging.info(f"Merged after CAAS join: {merged_df.shape}")

    if merged_df.empty:
        logging.warning(
            "No rows after global/CAAS join — global_csv may be empty (check aggregation step). "
            "Writing empty CT_ACCUMULATION randomization outputs and exiting cleanly."
        )
        _write_empty_outputs_and_exit(args)
        return

    unique_genes = merged_df['gene'].unique()
    gene_to_id   = {g: i for i, g in enumerate(unique_genes)}
    id_to_gene   = {i: g for g, i in gene_to_id.items()}
    merged_df['gene_id'] = merged_df['gene'].map(gene_to_id)
    genes_int = merged_df['gene_id'].values.astype(np.int32)
    n_genes = len(unique_genes)
    n_rows  = len(merged_df)

    positions     = merged_df['position'].values.astype(np.int64)
    masked        = merged_df['masked'].values.astype(bool)
    cons_idx_arr  = merged_df['cons_idx'].values.astype(np.float32)
    genes_int_arr = genes_int.astype(np.int32)
    caas_bool     = np.asarray(
        merged_df['iscaas'].replace({'TRUE': True, 'FALSE': False})
            .astype('boolean').fillna(False),
        dtype=bool
    )

    # Shared memory
    logging.info("Creating shared memory")
    SHM_positions = shm.SharedMemory(create=True, size=positions.nbytes)
    np.ndarray(positions.shape,     dtype=np.int64,   buffer=SHM_positions.buf)[:] = positions
    SHM_masked    = shm.SharedMemory(create=True, size=masked.nbytes)
    np.ndarray(masked.shape,        dtype=bool,        buffer=SHM_masked.buf)[:] = masked
    SHM_cons      = shm.SharedMemory(create=True, size=cons_idx_arr.nbytes)
    np.ndarray(cons_idx_arr.shape,  dtype=np.float32, buffer=SHM_cons.buf)[:] = cons_idx_arr
    SHM_genes     = shm.SharedMemory(create=True, size=genes_int_arr.nbytes)
    np.ndarray(genes_int_arr.shape, dtype=np.int32,   buffer=SHM_genes.buf)[:] = genes_int_arr
    SHM_caas      = shm.SharedMemory(create=True, size=caas_bool.nbytes)
    np.ndarray(caas_bool.shape,     dtype=bool,        buffer=SHM_caas.buf)[:] = caas_bool

    shm_names = {
        'positions': SHM_positions.name, 'masked': SHM_masked.name,
        'cons_idx':  SHM_cons.name,      'genes':  SHM_genes.name,
        'iscaas':    SHM_caas.name,
    }
    shapes = {
        'positions': positions.shape,     'masked': masked.shape,
        'cons_idx':  cons_idx_arr.shape,  'genes':  genes_int_arr.shape,
        'iscaas':    caas_bool.shape,
    }
    dtypes = {
        'positions': np.int64, 'masked': bool,
        'cons_idx':  np.float32, 'genes': np.int32, 'iscaas': bool,
    }

    # Actual counts
    caas_filter = (
        (merged_df['iscaas'] == True) &
        (merged_df['is_stable'] == True) &
        (merged_df['change_side'] != 'none') &
        (
            (merged_df['isAntiCAAS'] == False) |
            ((merged_df['isAntiCAAS'] == True) & (merged_df['ASR_ConsAntiCAAS'] == True))
        )
    )

    decile_bins = None
    if args.randomization_type == 'cons_decile':
        decile_bins = (
            _compute_bins_from_series(merged_df['cons_idx'])
            if args.decile_bins is None
            else np.array([float(x) for x in args.decile_bins.split(',')])
        )

    sig_mask    = caas_filter & (merged_df['is_significant'] == True)
    top_mask    = sig_mask & (merged_df['change_side'].isin(['top', 'both']))
    bottom_mask = sig_mask & (merged_df['change_side'].isin(['bottom', 'both']))
    both_mask   = sig_mask & (merged_df['change_side'] == 'both')

    convergent_patterns = ['convergent', 'parallel_convergence', 'parallel_mixed']
    # Exclude entries whose convergence support comes only from GS0.
    # Include GS1–GS4 (focal clade labelled groups) and rows where the same
    # position is convergent under US, GS1…GS4 but divergent under GS0 — those
    # appear as separate rows with a different caap_group tag.
    # NaN caap_group means a non-CAAS background row (left-join result); .fillna('')
    # ensures those rows are never excluded.
    not_gs0 = merged_df['caap_group'].fillna('').str.upper().ne('GS0')
    convergent_mask      = sig_mask & not_gs0 & (merged_df['pattern_type'].isin(convergent_patterns))
    convergent_caap_mask = sig_mask & not_gs0 & (
        (merged_df['pattern_type'].isin(convergent_patterns)) | (merged_df['iscaap'] == True)
    )
    convergent_caap_top_mask    = convergent_caap_mask & (merged_df['change_side'].isin(['top', 'both']))
    convergent_caap_bottom_mask = convergent_caap_mask & (merged_df['change_side'].isin(['bottom', 'both']))
    convergent_caap_both_mask   = convergent_caap_mask & (merged_df['change_side'] == 'both')

    def _bc(mask):
        return np.bincount(np.asarray(genes_int[np.asarray(mask, dtype=bool)], dtype=np.int32), minlength=n_genes)

    actual_all    = _bc(caas_filter)
    actual_sig    = _bc(sig_mask)
    actual_top    = _bc(top_mask)
    actual_bottom = _bc(bottom_mask)
    actual_both   = _bc(both_mask)
    actual_convergent         = _bc(convergent_mask)
    actual_convergent_caap    = _bc(convergent_caap_mask)
    actual_convergent_caap_top    = _bc(convergent_caap_top_mask)
    actual_convergent_caap_bottom = _bc(convergent_caap_bottom_mask)
    actual_convergent_caap_both   = _bc(convergent_caap_both_mask)

    _, caas_data = _build_caas_payload(merged_df, args.randomization_type, decile_bins)

    def _key_of(row):
        if args.randomization_type == 'naive':
            return 'global'
        return int(np.digitize([row['cons_idx']], bins=decile_bins[:-1], right=False)[0])

    extra_key_sizes = {}
    for name, m in [
        ('n_all', caas_filter), ('n_sig', sig_mask), ('n_top', top_mask),
        ('n_bottom', bottom_mask), ('n_both', both_mask),
        ('n_convergent', convergent_mask), ('n_convergent_caap', convergent_caap_mask),
        ('n_convergent_caap_top', convergent_caap_top_mask),
        ('n_convergent_caap_bottom', convergent_caap_bottom_mask),
        ('n_convergent_caap_both', convergent_caap_both_mask),
    ]:
        sub = merged_df.loc[m]
        for _, row in sub.iterrows():
            try:
                k = _key_of(row)
            except Exception:
                continue
            d = extra_key_sizes.setdefault(k, {})
            d[name] = d.get(name, 0) + 1

    position_to_tag = dict(zip(merged_df['position'].to_numpy(),
                               merged_df['tag'].astype(str).fillna('').to_numpy()))

    # Parallel plan
    total_rands = int(args.n_randomizations)
    workers     = args.workers if args.workers else os.cpu_count() or 1
    chunk_size  = max(1, total_rands // workers)
    chunks      = [(chunk_size, i) for i in range(total_rands // chunk_size)]
    if total_rands % chunk_size > 0:
        chunks.append((total_rands % chunk_size, len(chunks)))
    logging.info(f"Total randomizations: {total_rands} in {len(chunks)} chunks of ~{chunk_size}")

    output_dir = f"{args.output_prefix}_randomizations"
    results = []
    SHM_list = [SHM_positions, SHM_masked, SHM_cons, SHM_genes, SHM_caas]

    try:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=init_worker,
            initargs=(
                args.randomization_type, decile_bins, args.export_individual_rand,
                output_dir, args.global_seed, args.precompute_masks,
                n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
                caas_data, position_to_tag if args.export_individual_rand else None,
                actual_all, actual_sig, actual_top, actual_bottom, actual_both,
                actual_convergent, actual_convergent_caap,
                actual_convergent_caap_top, actual_convergent_caap_bottom,
                actual_convergent_caap_both,
            ),
        ) as executor:
            futures = [executor.submit(process_wrapper, ch) for ch in chunks]
            for fut in as_completed(futures):
                res = fut.result()
                results.append(res)
                logging.info(f"Completed chunk {res['chunk_idx']}")
    finally:
        for seg in SHM_list:
            try:
                seg.close(); seg.unlink()
            except FileNotFoundError:
                pass

    total_n_rands = sum(r['n_rands'] for r in results)
    fdr_thr = getattr(args, 'fdr_threshold', 0.05)
    gene_lists_dir = os.path.join(os.path.dirname(args.output_prefix) or '.', 'gene_lists')

    for cat, act in [
        ('all',    actual_all),
        ('sig',    actual_sig),
        ('top',    actual_top),
        ('bottom', actual_bottom),
        ('both',   actual_both),
        ('convergent',             actual_convergent),
        ('convergent_caap',        actual_convergent_caap),
        ('convergent_caap_top',    actual_convergent_caap_top),
        ('convergent_caap_bottom', actual_convergent_caap_bottom),
        ('convergent_caap_both',   actual_convergent_caap_both),
    ]:
        sum_cat    = np.sum([r[cat]['sum']         for r in results], axis=0)
        sumsq_cat  = np.sum([r[cat]['sum_sq']      for r in results], axis=0)
        cnt_ge_cat = np.sum([r[cat]['count_above'] for r in results], axis=0)

        mean_cat = sum_cat / total_n_rands if total_n_rands > 0 else np.zeros_like(sum_cat, dtype=float)
        p_emp    = (cnt_ge_cat + 1) / (total_n_rands + 1)

        var_cat  = np.clip((sumsq_cat / total_n_rands) - mean_cat ** 2, 0, None) if total_n_rands > 0 else np.zeros_like(sum_cat)
        sd_cat   = np.sqrt(var_cat)

        df_cat = pd.DataFrame({
            'GeneID': np.arange(n_genes, dtype=int),
            'Gene':   list(id_to_gene.values()),
            f'ActualCount_{cat}':     act,
            f'RandsAbove_{cat}':      cnt_ge_cat,
            f'RandsBelow_{cat}':      total_n_rands - cnt_ge_cat,
            f'PValueEmpirical_{cat}': p_emp,
            f'MeanSim_{cat}':         mean_cat,
            f'SD_{cat}':              sd_cat,
            f'NumRands_{cat}':        total_n_rands,
        })
        # Keep only genes that have at least one actual CAAS event; rows with
        # ActualCount == 0 carry no information and inflate FDR denominator.
        df_cat = df_cat[df_cat[f'ActualCount_{cat}'] > 0].reset_index(drop=True)

        # Apply BH-FDR on the remaining (all non-zero) empirical p-values.
        fdr_vals = np.ones(len(df_cat))
        if len(df_cat) > 0:
            fdr_vals = multipletests(
                df_cat[f'PValueEmpirical_{cat}'].values, method='fdr_bh'
            )[1]
        df_cat[f'PValueEmpirical_{cat}_FDR'] = fdr_vals

        base_dir = os.path.dirname(args.output_prefix) or '.'
        os.makedirs(base_dir, exist_ok=True)
        out_cat = f"{args.output_prefix}_{cat}_aggregated_results.csv"
        df_cat.to_csv(out_cat, index=False,
                      compression='gzip' if args.compress else None)
        logging.info(f"Results ({cat}): {out_cat}")

        _write_gene_lists(df_cat, cat, gene_lists_dir, fdr_thr)

    # Decile window info
    if args.randomization_type == 'cons_decile':
        df_elig = merged_df[merged_df['masked'] == False].copy()
        bins = decile_bins if decile_bins is not None else _compute_bins_from_series(
            merged_df.loc[caas_filter, 'cons_idx']
        )
        df_elig['decile'] = np.digitize(df_elig['cons_idx'], bins=bins[:-1], right=False)
        global_decile     = df_elig.groupby('decile').size().reset_index(name='NumPositions')
        per_gene_decile   = df_elig.groupby(['gene', 'decile']).size().reset_index(name='NumPositions')
        per_gene_pivot    = per_gene_decile.pivot(
            index='gene', columns='decile', values='NumPositions'
        ).fillna(0).reset_index()
        global_decile.to_parquet(f"{args.output_prefix}_cons_decile_global.parquet",   index=False)
        per_gene_pivot.to_parquet(f"{args.output_prefix}_cons_decile_per_gene.parquet", index=False)

    logging.info("Randomization analysis complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CT_ACCUMULATION: permutation significance test')
    parser.add_argument('--global_csv',   required=True)
    parser.add_argument('--caas_csv',     required=True)
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--randomization-type', choices=['naive', 'cons_decile'], required=True)
    parser.add_argument('--n-randomizations', type=int, default=100000)
    parser.add_argument('--workers',      type=int, default=None)
    parser.add_argument('--compress',     action='store_true')
    parser.add_argument('--export-individual-rand', action='store_true')
    parser.add_argument('--decile-bins',  type=str, default=None)
    parser.add_argument('--global-seed',  type=int, default=None)
    parser.add_argument('--fdr-threshold', type=float, default=0.05,
                        help='FDR threshold for gene list export (default: 0.05)')
    parser.add_argument('--precompute-masks', dest='precompute_masks', action='store_true')
    parser.add_argument('--no-precompute-masks', dest='precompute_masks', action='store_false')
    parser.set_defaults(precompute_masks=True)
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')
    main(args)
