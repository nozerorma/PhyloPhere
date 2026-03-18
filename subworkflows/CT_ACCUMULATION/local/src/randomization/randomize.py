#!/usr/bin/env python3
"""
Randomization / permutation significance test — adapted for CT_ACCUMULATION in PhyloPhere.

Changes vs. original to_integrate version:
  - _remap_caas_df(): remaps filtered_discovery.tsv column names to the internal schema
    (Gene→gene, Position→msa_pos, sig_both→is_significant,
     Pattern→pattern_type, CAAP_Group→caap_group/iscaap).
  - Gene list categories mirror CT_postproc.Rmd focal_sets chunk (both GS0-included and
    GS0-excluded variants):
      global, global_significant, top_significant, bottom_significant
      convergent_pool (+top/bottom)    — mirrors site_flags (GS0 included in gs_nondiv)
      divergent_pool  (+top/bottom)
      convergent_pool_no_gs0 (+top/bottom) — mirrors site_flags_no_gs0
      divergent_pool_no_gs0  (+top/bottom)
      us_significant (+top/bottom), us_nondiv (+top/bottom), us_div (+top/bottom)
      gs0_significant (+top/bottom), gs0_nondiv (+top/bottom), gs0_div (+top/bottom)
      gs{1-4}_significant (+top/bottom), gs{1-4}_nondiv (+top/bottom), gs{1-4}_div (+top/bottom)
  - GS0 is hard-excluded from valid_row (dubious-conserved + GS0 exclusions) but a
    separate valid_row_any mask (dubious-conserved only) gates GS0 and pool categories.
  - actual_counts dict replaces individual actual_* arrays; WorkerClass and init_worker
    take a single actual_counts dict; process_chunk uses dict-driven accumulators.
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
    """Normalise a filtered_discovery.tsv DataFrame to the internal column schema.

    Source-of-truth columns (tab-separated):
      Gene, Position, tag, caas, is_significant, Pvalue, pvalue_boot, Pattern,
      convergence_description, convergence_mode, CAAP_Group, amino_encoded,
      is_conserved_meta, conserved_pair, sig_hyp, sig_perm, sig_both,
      top_change_type, bottom_change_type, change_side, low_confidence_nodes,
      asr_is_conserved, comments, ..., Trait

    Produces internal schema columns:
      gene, msa_pos, is_significant (from sig_both), tag, pattern_type (from Pattern),
      caap_group (from CAAP_Group, uppercased), iscaap,
      change_side, is_conserved_meta, asr_is_conserved
    """
    # Rename Gene→gene, Position→msa_pos
    df = df.rename(columns={'Gene': 'gene', 'Position': 'msa_pos'})

    _bool = lambda x: str(x).strip().lower() in {'true', 't', '1', 'yes', 'y'}

    # is_significant: use sig_both (passes both hyp + perm tests)
    df['is_significant'] = df['sig_both'].map(_bool)

    # Coerce conserved-state columns to bool
    df['is_conserved_meta'] = df['is_conserved_meta'].map(_bool)
    df['asr_is_conserved']  = df['asr_is_conserved'].map(_bool)
<<<<<<< HEAD
=======
    if 'asr_root_conserved' in df.columns:
        df['asr_root_conserved'] = df['asr_root_conserved'].map(_bool)
    else:
        df['asr_root_conserved'] = False
>>>>>>> asr

    # pattern_type: from Pattern column
    df['pattern_type'] = df['Pattern']

    # caap_group: preserve the raw CAAP_Group label (GS0–GS4, US)
    df['caap_group'] = df['CAAP_Group'].astype(str).str.strip().str.upper()

    # iscaap: any labelled group other than 'US' (includes GS0–GS4)
    df['iscaap'] = df['caap_group'].apply(lambda x: x != 'US')

    logging.info("CAAS CSV: remapped filtered_discovery.tsv columns to internal schema")
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
        actual_counts,          # dict: category_name -> np.ndarray(n_genes, int64)
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
        self.actual_counts = {k: v.astype(np.int64, copy=False) for k, v in actual_counts.items()}

        if self.global_seed is not None:
            np.random.seed(self.global_seed)
            random.seed(self.global_seed)

        _zero_counts = {f'n_{cat}': 0 for cat in actual_counts}
        self.caas_by_key_counts = defaultdict(lambda: dict(_zero_counts))
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

        # Vectorized accumulators per category (dict-driven).
        # Each entry is a list [sum, sum_sq, count_above] — mutable so in-place += works.
        def _zeros():
            return [
                np.zeros(self.n_genes, dtype=np.int64),
                np.zeros(self.n_genes, dtype=np.float64),
                np.zeros(self.n_genes, dtype=np.int64),
            ]
        S = {cat: _zeros() for cat in self.actual_counts}

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
            z = lambda: np.zeros(self.n_genes, dtype=np.int64)
            c = {cat: z() for cat in self.actual_counts}
            uid = f"{os.getpid()}_{chunk_idx}_{r}" if self.export_individual else None

            for key, cinfo in self.caas_by_key_counts.items():
                elig = self.eligible_by_key.get(key)
                if elig is None or elig.size == 0:
                    continue
                for cat, carr in c.items():
                    n = cinfo.get(f'n_{cat}', 0)
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

            def _upd(Sv, cv, act):
                Sv[0] += cv
                Sv[1] += cv.astype(np.float64) ** 2
                Sv[2] += (cv >= act).astype(np.int64)
            for cat in S:
                _upd(S[cat], c[cat], self.actual_counts[cat])

        if self.export_individual and batch_rows and writer is not None:
            writer.write_table(pa.Table.from_pandas(pd.DataFrame(batch_rows), schema=schema))
        if writer is not None:
            writer.close()

        def _pack(Sv):
            return {'sum': Sv[0], 'sum_sq': Sv[1], 'count_above': Sv[2]}
        return {
            'chunk_idx': chunk_idx, 'n_rands': chunk_size,
            **{cat: _pack(S[cat]) for cat in S},
        }


# ---------------------------
# Worker helpers (module-level for pickling)
# ---------------------------

def init_worker(
    randomization_type, decile_bins, export_individual, output_dir, global_seed,
    precompute_masks, n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
    caas_data, position_to_tag,
    actual_counts,              # dict: category_name -> np.ndarray
):
    global worker
    worker = RandomizationWorker(
        randomization_type, decile_bins, export_individual, output_dir, global_seed,
        precompute_masks, n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
        caas_data, position_to_tag,
        actual_counts,
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

def _normalize_cat_for_filename(cat):
    """Normalize category token for export filenames.

    Convention:
      - *_significant* -> *_sig*
      - keep direction suffixes as-is (top/bottom) and use explicit nodir categories
    """
    return str(cat).replace('_significant', '_sig')


def _write_gene_lists(df_cat, cat, gene_lists_dir, fdr_threshold):
    """Write a single FDR-corrected gene-list .txt file for one category."""
    os.makedirs(gene_lists_dir, exist_ok=True)
    col = f'PValueEmpirical_{cat}_FDR'
    sig = df_cat[
        (df_cat[f'ActualCount_{cat}'] > 0) &
        (df_cat[col] < fdr_threshold)
    ]['Gene']
    file_cat = _normalize_cat_for_filename(cat)
    out_path = os.path.join(gene_lists_dir, f"{file_cat}_fdr.txt")
    sig.to_csv(out_path, index=False, header=False)
    logging.info(f"Gene list ({cat}, fdr): {len(sig)} genes → {out_path}")


def _write_empty_outputs_and_exit(args):
    """Write empty-but-valid outputs when no rows are available for randomization."""
    # Full category list mirrors CT_postproc focal_sets (both GS0-included and GS0-excluded).
    _grp_cats = lambda g: [
        f'{g}_significant', f'{g}_significant_top', f'{g}_significant_bottom',
        f'{g}_nondiv',      f'{g}_nondiv_top',      f'{g}_nondiv_bottom',
        f'{g}_div',         f'{g}_div_top',          f'{g}_div_bottom',
    ]
    categories = (
        ['global', 'global_significant', 'top_significant', 'bottom_significant']
        + ['convergent_pool', 'convergent_pool_top', 'convergent_pool_bottom']
        + ['divergent_pool',  'divergent_pool_top',  'divergent_pool_bottom']
        + ['convergent_pool_no_gs0', 'convergent_pool_no_gs0_top', 'convergent_pool_no_gs0_bottom']
        + ['divergent_pool_no_gs0',  'divergent_pool_no_gs0_top',  'divergent_pool_no_gs0_bottom']
        + _grp_cats('us')
        + _grp_cats('gs0')
        + _grp_cats('gs1')
        + _grp_cats('gs2')
        + _grp_cats('gs3')
        + _grp_cats('gs4')
    )

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

        # Keep output contract: create empty FDR gene list files for every category.
        file_cat = _normalize_cat_for_filename(cat)
        out_path = os.path.join(gene_lists_dir, f"{file_cat}_fdr.txt")
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
    required_cols = ['gene', 'msa_pos', 'is_significant',
                     'change_side', 'tag', 'pattern_type', 'iscaap', 'caap_group',
<<<<<<< HEAD
                     'is_conserved_meta', 'asr_is_conserved']
=======
                     'is_conserved_meta', 'asr_is_conserved', 'asr_root_conserved']
>>>>>>> asr
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
        (merged_df['change_side'] != 'none')
    )

    decile_bins = None
    if args.randomization_type == 'cons_decile':
        decile_bins = (
            _compute_bins_from_series(merged_df['cons_idx'])
            if args.decile_bins is None
            else np.array([float(x) for x in args.decile_bins.split(',')])
        )

    sig_mask    = caas_filter & (merged_df['is_significant'] == True)

    # --- Row-level exclusions: GS0 and dubious-conserved (is_conserved_meta=T & asr_is_conserved=F) ---
    _boolify = lambda s: str(s).strip().lower() in {'true', 't', '1', 'yes', 'y'}
    if 'is_conserved_meta' in merged_df.columns:
        _icm = merged_df['is_conserved_meta'].map(_boolify).fillna(False)
    else:
        _icm = pd.Series(False, index=merged_df.index)
    if 'asr_is_conserved' in merged_df.columns:
        _arc = merged_df['asr_is_conserved'].map(_boolify).fillna(True)
    else:
        _arc = pd.Series(True, index=merged_df.index)
<<<<<<< HEAD
    dubious_conserved = _icm & ~_arc
    is_gs0    = merged_df['caap_group'].fillna('').str.strip().str.upper().eq('GS0')
    valid_row     = ~dubious_conserved & ~is_gs0   # gate for US + GS1-4
    valid_row_any = ~dubious_conserved             # gate for GS0 and pool categories
    logging.info(
        f"Row exclusions: {dubious_conserved.sum()} dubious-conserved "
        f"(is_conserved_meta=T & asr_is_conserved=F), {is_gs0.sum()} GS0 — "
=======
    if getattr(args, 'use_all_mrca_filter', False) and 'asr_root_conserved' in merged_df.columns:
        _arrc = merged_df['asr_root_conserved'].map(_boolify).fillna(True)
    elif getattr(args, 'use_all_mrca_filter', False):
        _arrc = pd.Series(True, index=merged_df.index)
    else:
        _arrc = pd.Series(True, index=merged_df.index)
    dubious_conserved = _icm & ~_arc
    root_filtered = _icm & ~_arrc
    is_gs0    = merged_df['caap_group'].fillna('').str.strip().str.upper().eq('GS0')
    valid_row     = ~dubious_conserved & ~root_filtered & ~is_gs0   # gate for US + GS1-4
    valid_row_any = ~dubious_conserved & ~root_filtered             # gate for GS0 and pool categories
    logging.info(
        f"Row exclusions: {dubious_conserved.sum()} dubious-conserved "
        f"(is_conserved_meta=T & asr_is_conserved=F), "
        f"{root_filtered.sum()} root-filtered "
        f"(is_conserved_meta=T & asr_root_conserved=F), {is_gs0.sum()} GS0 — "
>>>>>>> asr
        f"{valid_row.sum()} rows remain (valid_row), {valid_row_any.sum()} (valid_row_any)"
    )

    # Propagate non-GS0 exclusions to base caas_filter and sig_mask.
    caas_filter   = caas_filter & valid_row
    sig_mask      = caas_filter & (merged_df['is_significant'] == True)

    # GS0-inclusive caas/sig masks (dubious-conserved still excluded).
    caas_filter_any = (
        (merged_df['iscaas'] == True) &
        (merged_df['change_side'] != 'none')
    ) & valid_row_any
    sig_mask_any = caas_filter_any & (merged_df['is_significant'] == True)

    # --- Per-row group and pattern helpers ---
    _gu      = merged_df['caap_group'].fillna('').str.strip().str.upper()
    is_nondiv = merged_df['pattern_type'].fillna('').ne('divergent')
    is_div    = merged_df['pattern_type'].fillna('').eq('divergent')
    _cs       = merged_df['change_side']

    def _top(m):    return m & _cs.isin(['top', 'both'])
    def _bottom(m): return m & _cs.isin(['bottom', 'both'])

    # --- Global: sig_both in any valid non-GS0 group ---
    global_mask     = caas_filter
    global_sig_mask = sig_mask
    top_sig_mask    = _top(global_sig_mask)
    bottom_sig_mask = _bottom(global_sig_mask)

    # --- US masks ---
    us_sig_mask        = sig_mask & (_gu == 'US')
    us_nondiv_mask     = us_sig_mask & is_nondiv
    us_div_mask        = us_sig_mask & is_div

    # --- GS0 masks (valid_row_any: dubious-conserved excluded, GS0 allowed) ---
    gs0_sig_mask       = sig_mask_any & (_gu == 'GS0')
    gs0_nondiv_mask    = gs0_sig_mask & is_nondiv
    gs0_div_mask       = gs0_sig_mask & is_div

    # --- GS1-4 masks ---
    gs1_sig_mask       = sig_mask & (_gu == 'GS1')
    gs1_nondiv_mask    = gs1_sig_mask & is_nondiv
    gs1_div_mask       = gs1_sig_mask & is_div

    gs2_sig_mask       = sig_mask & (_gu == 'GS2')
    gs2_nondiv_mask    = gs2_sig_mask & is_nondiv
    gs2_div_mask       = gs2_sig_mask & is_div

    gs3_sig_mask       = sig_mask & (_gu == 'GS3')
    gs3_nondiv_mask    = gs3_sig_mask & is_nondiv
    gs3_div_mask       = gs3_sig_mask & is_div

    gs4_sig_mask       = sig_mask & (_gu == 'GS4')
    gs4_nondiv_mask    = gs4_sig_mask & is_nondiv
    gs4_div_mask       = gs4_sig_mask & is_div

    # --- Pool masks: mirrors compute_site_flags() in CT_postproc ---
    # convergent_pool (GS0 included): us_nondiv | gs_nondiv  (gs = any non-US incl. GS0)
    convergent_pool_mask        = us_nondiv_mask | gs0_nondiv_mask | gs1_nondiv_mask | gs2_nondiv_mask | gs3_nondiv_mask | gs4_nondiv_mask
    divergent_pool_mask         = us_div_mask    | gs0_div_mask    | gs1_div_mask    | gs2_div_mask    | gs3_div_mask    | gs4_div_mask
    # convergent_pool_no_gs0: mirrors site_flags_no_gs0
    convergent_pool_no_gs0_mask = us_nondiv_mask | gs1_nondiv_mask | gs2_nondiv_mask | gs3_nondiv_mask | gs4_nondiv_mask
    divergent_pool_no_gs0_mask  = us_div_mask    | gs1_div_mask    | gs2_div_mask    | gs3_div_mask    | gs4_div_mask

    logging.info(
        f"Category sizes: global_sig={global_sig_mask.sum()}, "
        f"convergent_pool={convergent_pool_mask.sum()}, convergent_pool_no_gs0={convergent_pool_no_gs0_mask.sum()}, "
        f"us_sig={us_sig_mask.sum()}, us_nondiv={us_nondiv_mask.sum()}, us_div={us_div_mask.sum()}, "
        f"gs0_sig={gs0_sig_mask.sum()}, gs1_sig={gs1_sig_mask.sum()}, gs2_sig={gs2_sig_mask.sum()}, "
        f"gs3_sig={gs3_sig_mask.sum()}, gs4_sig={gs4_sig_mask.sum()}"
    )

    def _bc(mask):
        return np.bincount(np.asarray(genes_int[np.asarray(mask, dtype=bool)], dtype=np.int32), minlength=n_genes)

    # Ordered dict: category_name → actual count array.  Mirrors CT_postproc focal_sets tables.
    # Ordered dict (insertion order guaranteed in Python 3.7+): category_name → mask.
    # Mirrors CT_postproc focal_sets (both GS0-included and GS0-excluded variants).
    actual_counts_masks = dict([
        # ── Global pools ────────────────────────────────────────────────────────────
        ('global',                          global_mask),
        ('global_significant',              global_sig_mask),
        ('top_significant',                 top_sig_mask),
        ('bottom_significant',              bottom_sig_mask),
        # ── Convergent / divergent pools (GS0 included) ─────────────────────────────
        ('convergent_pool',                 convergent_pool_mask),
        ('convergent_pool_top',             _top(convergent_pool_mask)),
        ('convergent_pool_bottom',          _bottom(convergent_pool_mask)),
        ('divergent_pool',                  divergent_pool_mask),
        ('divergent_pool_top',              _top(divergent_pool_mask)),
        ('divergent_pool_bottom',           _bottom(divergent_pool_mask)),
        # ── Convergent / divergent pools (GS0 excluded) ─────────────────────────────
        ('convergent_pool_no_gs0',          convergent_pool_no_gs0_mask),
        ('convergent_pool_no_gs0_top',      _top(convergent_pool_no_gs0_mask)),
        ('convergent_pool_no_gs0_bottom',   _bottom(convergent_pool_no_gs0_mask)),
        ('divergent_pool_no_gs0',           divergent_pool_no_gs0_mask),
        ('divergent_pool_no_gs0_top',       _top(divergent_pool_no_gs0_mask)),
        ('divergent_pool_no_gs0_bottom',    _bottom(divergent_pool_no_gs0_mask)),
        # ── US ──────────────────────────────────────────────────────────────────────
        ('us_significant',                  us_sig_mask),
        ('us_significant_top',              _top(us_sig_mask)),
        ('us_significant_bottom',           _bottom(us_sig_mask)),
        ('us_nondiv',                       us_nondiv_mask),
        ('us_nondiv_top',                   _top(us_nondiv_mask)),
        ('us_nondiv_bottom',                _bottom(us_nondiv_mask)),
        ('us_div',                          us_div_mask),
        ('us_div_top',                      _top(us_div_mask)),
        ('us_div_bottom',                   _bottom(us_div_mask)),
        # ── GS0 ─────────────────────────────────────────────────────────────────────
        ('gs0_significant',                 gs0_sig_mask),
        ('gs0_significant_top',             _top(gs0_sig_mask)),
        ('gs0_significant_bottom',          _bottom(gs0_sig_mask)),
        ('gs0_nondiv',                      gs0_nondiv_mask),
        ('gs0_nondiv_top',                  _top(gs0_nondiv_mask)),
        ('gs0_nondiv_bottom',               _bottom(gs0_nondiv_mask)),
        ('gs0_div',                         gs0_div_mask),
        ('gs0_div_top',                     _top(gs0_div_mask)),
        ('gs0_div_bottom',                  _bottom(gs0_div_mask)),
        # ── GS1 ─────────────────────────────────────────────────────────────────────
        ('gs1_significant',                 gs1_sig_mask),
        ('gs1_significant_top',             _top(gs1_sig_mask)),
        ('gs1_significant_bottom',          _bottom(gs1_sig_mask)),
        ('gs1_nondiv',                      gs1_nondiv_mask),
        ('gs1_nondiv_top',                  _top(gs1_nondiv_mask)),
        ('gs1_nondiv_bottom',               _bottom(gs1_nondiv_mask)),
        ('gs1_div',                         gs1_div_mask),
        ('gs1_div_top',                     _top(gs1_div_mask)),
        ('gs1_div_bottom',                  _bottom(gs1_div_mask)),
        # ── GS2 ─────────────────────────────────────────────────────────────────────
        ('gs2_significant',                 gs2_sig_mask),
        ('gs2_significant_top',             _top(gs2_sig_mask)),
        ('gs2_significant_bottom',          _bottom(gs2_sig_mask)),
        ('gs2_nondiv',                      gs2_nondiv_mask),
        ('gs2_nondiv_top',                  _top(gs2_nondiv_mask)),
        ('gs2_nondiv_bottom',               _bottom(gs2_nondiv_mask)),
        ('gs2_div',                         gs2_div_mask),
        ('gs2_div_top',                     _top(gs2_div_mask)),
        ('gs2_div_bottom',                  _bottom(gs2_div_mask)),
        # ── GS3 ─────────────────────────────────────────────────────────────────────
        ('gs3_significant',                 gs3_sig_mask),
        ('gs3_significant_top',             _top(gs3_sig_mask)),
        ('gs3_significant_bottom',          _bottom(gs3_sig_mask)),
        ('gs3_nondiv',                      gs3_nondiv_mask),
        ('gs3_nondiv_top',                  _top(gs3_nondiv_mask)),
        ('gs3_nondiv_bottom',               _bottom(gs3_nondiv_mask)),
        ('gs3_div',                         gs3_div_mask),
        ('gs3_div_top',                     _top(gs3_div_mask)),
        ('gs3_div_bottom',                  _bottom(gs3_div_mask)),
        # ── GS4 ─────────────────────────────────────────────────────────────────────
        ('gs4_significant',                 gs4_sig_mask),
        ('gs4_significant_top',             _top(gs4_sig_mask)),
        ('gs4_significant_bottom',          _bottom(gs4_sig_mask)),
        ('gs4_nondiv',                      gs4_nondiv_mask),
        ('gs4_nondiv_top',                  _top(gs4_nondiv_mask)),
        ('gs4_nondiv_bottom',               _bottom(gs4_nondiv_mask)),
        ('gs4_div',                         gs4_div_mask),
        ('gs4_div_top',                     _top(gs4_div_mask)),
        ('gs4_div_bottom',                  _bottom(gs4_div_mask)),
    ])
    actual_counts = {cat: _bc(mask) for cat, mask in actual_counts_masks.items()}

    _, caas_data = _build_caas_payload(merged_df, args.randomization_type, decile_bins)

    def _key_of(row):
        if args.randomization_type == 'naive':
            return 'global'
        return int(np.digitize([row['cons_idx']], bins=decile_bins[:-1], right=False)[0])

    extra_key_sizes = {}
    for cat, m in actual_counts_masks.items():
        name = f'n_{cat}'
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
                actual_counts,
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

    for cat, act in actual_counts.items():
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
        # I'm removing this criteria for now but it is worth considering.
        # df_cat = df_cat[df_cat[f'ActualCount_{cat}'] > 0].reset_index(drop=True)

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
    parser.add_argument('--use-all-mrca-filter', action='store_true',
                        help='Exclude conserved-meta rows failing asr_root_conserved')
    parser.add_argument('--precompute-masks', dest='precompute_masks', action='store_true')
    parser.add_argument('--no-precompute-masks', dest='precompute_masks', action='store_false')
    parser.set_defaults(precompute_masks=True)
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')
    main(args)
