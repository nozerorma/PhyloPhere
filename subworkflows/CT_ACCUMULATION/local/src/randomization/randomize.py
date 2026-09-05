#!/usr/bin/env python3
"""
Randomization / permutation analysis for CT_ACCUMULATION in PhyloPhere.

Exports aggregated randomization outputs for five per-group tests:
  - us, gs1, gs2, gs3, gs4

Test design
-----------
Two randomization families, selected by --randomization-type:

'naive' / 'cons_decile' — an occupancy null. For each group, the OBSERVED
number of CAAS positions is redrawn uniformly (without replacement) from the
eligible position pool, and the per-gene counts are recounted. Because the
pool is position-level, gene length (strictly, each gene's number of eligible
positions) is controlled for by construction.

'permulation' — an excess null, opt-in (see run_permulation_null below). The
randomised CAAS set for a gene at cycle i is that gene's ACTUAL detections
replayed from a prior CAAS_PERMULATION run, not a synthetic draw from a pool;
the per-gene count is free per cycle rather than pinned to the observed count.
It holds the phenotype-tree confound but is NOT conservation-decile matched.

In every case the empirical p-value per gene is
    p = (#{null count >= observed count} + 1) / (N + 1)
so a gene is significant when it accumulates more CAAS than its null predicts.

Eligible pool
-------------
The pool is the set of positions an observed CAAS could actually have come from:
the positions CAAStools tested (background.output), restricted to the genes
surviving post-processing (cleaned_background_main.txt). This mirrors posenrich's
build_background(), and the restriction matters in both directions — a position
that failed the gap/missingness/pattern filters can never be a CAAS, and neither
can any position in a gene post-processing removed.

Getting this set right is what makes the per-gene p-values comparable to one
another. A gene's null expectation is proportional to its share of the pool, so
the pool must be proportional to each gene's genuine opportunity to carry a CAAS.
Widening it to every ungapped alignment column would break that proportionality
gene by gene: on cancer_complete_bm the tested/ungapped ratio ranges from 0.06
(p10) through 0.24 (median) to 0.51 (p90), so poorly covered genes would be
handed a ~4x too-high expectation (conservative) and well covered genes a ~2x
too-low one (anti-conservative), concentrating false positives in the
best-covered genes.

Per-group draws are independent
-------------------------------
Each group draws separately, which is what keeps each group's MARGINAL p-value
correctly calibrated. It does not make the five p-values independent of one
another: a single physical position can be a CAAS under several nested schemes,
so the OBSERVED counts — and hence the p-values — are positively correlated
whatever the null does. Evidence is therefore combined across groups with the
Cauchy Combination Test in scoring_compute.R, which is valid under arbitrary
dependence between the combined p-values.

No significance/convergence/divergence category families are exported.
No FDR or gene-list outputs are generated here.
"""

import argparse
import gzip
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import random
import multiprocessing.shared_memory as shm


# ---------------------------
# Eligible position pool
# ---------------------------

def load_universe_genes(path):
    """Post-processing-surviving gene list (cleaned_background_main.txt)."""
    genes = set()
    with open(path) as f:
        for line in f:
            g = line.strip()
            if g and g != "Gene":
                genes.add(g)
    return genes


def load_tested_positions(background_file, universe_genes=None):
    """Positions CAAStools actually tested, as a set of (gene, msa_pos) pairs.

    Format mirrors posenrich's build_background(): `gene<TAB>comma-separated
    positions`, one line per gene. The positions are in the SAME coordinate space
    as filtered_discovery's `Position` — i.e. `msa_pos` after _remap_caas_df, which
    is also the key the global/CAAS merge joins on. (Verified against real output:
    every observed CAAS position is present in this file.) Do NOT match on the
    global table's `position` column — that is a globally unique running row index,
    not an alignment coordinate.

    Restricting to `universe_genes` drops positions belonging to genes that
    post-processing removed, which can never contribute an observed CAAS.
    """
    tested = set()
    with open(background_file) as f:
        for line in f:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if len(parts) < 2:
                continue
            gene, positions = parts[0], parts[1]
            if gene == "Gene":
                continue
            if universe_genes is not None and gene not in universe_genes:
                continue
            if positions in ("", "NULL", "Position"):
                continue
            for p in positions.split(","):
                p = p.strip()
                if not p:
                    continue
                try:
                    tested.add((gene, int(p)))
                except ValueError:
                    continue
    return tested


# ---------------------------
# Column remapping
# ---------------------------

def _remap_caas_df(df):
    """Normalise a filtered_discovery.tsv DataFrame to the internal column schema.

    Source-of-truth columns (tab-separated):
      Gene, Position, tag, caas, recovery_boot, convergence_type, caap_group,
      amino_encoded, is_conserved_meta, conserved_pair, change_top, change_bottom,
      change_side, asr_is_conserved, ..., Trait

    Produces internal schema columns:
      gene, msa_pos, tag, convergence_type, caap_group (uppercased),
      iscaap, change_side, is_conserved_meta, asr_is_conserved
    """
    # Rename Gene→gene, Position→msa_pos (structural keys; concept columns are
    # already in disambiguation's canonical lowercase form).
    df = df.rename(columns={'Gene': 'gene', 'Position': 'msa_pos'})

    _bool = lambda x: str(x).strip().lower() in {'true', 't', '1', 'yes', 'y'}

    # Coerce conserved-state columns to bool
    df['is_conserved_meta'] = df['is_conserved_meta'].map(_bool)
    # convergence_type is already canonical in filtered_discovery.tsv — no remap needed.

    # caap_group: normalise casing of the raw label (GS1–GS4, US)
    df['caap_group'] = df['caap_group'].astype(str).str.strip().str.upper()

    # iscaap: any labelled group other than 'US' (includes GS1–GS4)
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
            # np.digitize returns 1-based indices [1, 10] for percentiles.
            # Clip and shift to 0-based indices [0, 9] mapping to deciles.
            dec = np.digitize(self.cons_idx, bins=self.decile_bins[:-1], right=False)
            dec = np.clip(dec - 1, 0, len(self.decile_bins) - 2)
            for d in range(len(self.decile_bins) - 1):
                mask = (dec == d) & ~self.masked
                self.eligible_by_key[d] = self.row_indices[mask]
            return
        raise ValueError(f"Unknown randomization_type: {self.randomization_type}")

    def process_chunk(self, chunk_size, chunk_idx):
        logging.info(f"Starting chunk {chunk_idx} with {chunk_size} randomizations in PID {os.getpid()}")
        seed = None
        if self.global_seed is not None:
            seed = ((self.global_seed or 0) + 9973 * chunk_idx) & 0x7FFFFFFF
        rng = np.random.default_rng(seed)

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

        # Pre-allocate slice buffers for fast indexing
        cat_total_sizes = {cat: self.actual_counts[cat].sum() for cat in self.actual_counts}
        idx_buffers = {cat: np.empty(cat_total_sizes[cat], dtype=np.int32) for cat in self.actual_counts}

        offsets = {cat: 0 for cat in self.actual_counts}
        decile_plans = []

        for key, cinfo in self.caas_by_key_counts.items():
            elig = self.eligible_by_key.get(key)
            if elig is None or elig.size == 0:
                continue

            # Build per-category (dest_start, dest_end) slices; skip keys with no positions.
            cat_slices = {}
            for cat in ['us', 'gs1', 'gs2', 'gs3', 'gs4']:
                n_cat = cinfo.get(f'n_{cat}', 0)
                if n_cat > 0:
                    cat_slices[cat] = (offsets[cat], offsets[cat] + n_cat)
                    offsets[cat] += n_cat
            if not cat_slices:
                continue
            decile_plans.append((elig, cat_slices, key))

        for r in range(chunk_size):
            uid = f"{os.getpid()}_{chunk_idx}_{r}" if self.export_individual else None

            for elig, cat_slices, key in decile_plans:
                # Each group draws its own positions from this decile's eligible
                # set, rather than partitioning one shared draw. That keeps every
                # group's marginal p-value calibrated against the same pool its
                # observed count came from. Cross-group dependence (a position can
                # be a CAAS under several schemes) is a property of the observed
                # data and is handled at combination time — see the module docstring.
                for cat, (dest_start, dest_end) in cat_slices.items():
                    n_cat = dest_end - dest_start
                    if len(elig) < n_cat:
                        idx_cat = rng.choice(elig, size=n_cat, replace=True)
                    else:
                        idx_cat = rng.choice(elig, size=n_cat, replace=False)
                    idx_buffers[cat][dest_start:dest_end] = idx_cat

                if self.export_individual:
                    caas_list = self.caas_by_key_lists.get(key, [])
                    if caas_list:
                        if len(elig) < len(caas_list):
                            idx_ind = rng.choice(elig, size=len(caas_list), replace=True)
                        else:
                            idx_ind = rng.choice(elig, size=len(caas_list), replace=False)
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

            for cat in self.actual_counts:
                cv = bincount(self.genes[idx_buffers[cat]], minlength=self.n_genes)
                S[cat][0] += cv
                S[cat][1] += cv.astype(np.float64) ** 2
                S[cat][2] += (cv >= self.actual_counts[cat]).astype(np.int64)

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

def _build_caas_payload(merged_df, randomization_type, decile_bins, pool_mask):
    caas_data = []
    for _, row in merged_df.loc[pool_mask].iterrows():
        if randomization_type == 'naive':
            key = 'global'
        elif randomization_type == 'cons_decile':
            k = int(np.digitize([row['cons_idx']], bins=decile_bins[:-1], right=False)[0])
            key = np.clip(k - 1, 0, len(decile_bins) - 2)
        else:
            continue
        caas_data.append({'position': row['position'], 'key': key})
    unique_keys = list({item['key'] for item in caas_data})
    return unique_keys, caas_data


# ---------------------------
# Permulation null (Tier 3E)
# ---------------------------
#
# Unlike 'naive'/'cons_decile', this type does not draw synthetic positions from
# an eligible pool. The randomised CAAS set for a gene at cycle i IS that gene's
# actually-detected positions in the CAAS permulation replay (caas_permulation.nf /
# gene_wrapper.py's run_permulation) at cycle i, subject to the same change_side
# filter the observed pool uses. This holds the phenotype-tree confound (the same
# thing the permulation null holds everywhere else in the pipeline) but is NOT
# conservation-decile matched like 'cons_decile' -- a different question, not a
# stricter version of it. The per-gene CAAS count is free per cycle (not fixed to
# the observed count), so this is an excess null like the rest of the permulation
# family, not an occupancy null.

_PERM_CATS = ['us', 'gs1', 'gs2', 'gs3', 'gs4']
_PERM_GROUP_OF_CAT = {'us': 'US', 'gs1': 'GS1', 'gs2': 'GS2', 'gs3': 'GS3', 'gs4': 'GS4'}
_PERM_CAT_OF_GROUP = {v: k for k, v in _PERM_GROUP_OF_CAT.items()}


def _iter_perm_detail_rows(detail_path):
    """Minimal, self-contained reader for perm_pos_detail: either a per-gene shard
    directory (perm_pos_detail/<Gene>.tsv.gz, current) or a single concatenated
    perm_pos_detail.tsv.gz (legacy). Mirrors CT_DISAMBIGUATION's gene_wrapper.py
    iter_detail_rows(), duplicated here rather than imported because
    CT_ACCUMULATION stages its own `local/` tree independently (ctacc_run.nf
    does `cp -R local/*`, not a cross-subworkflow import)."""
    import csv as _csv

    p = Path(detail_path)
    if p.is_dir():
        shards = sorted(p.glob("*.tsv.gz"))
        if not shards:
            raise FileNotFoundError(f"[permulation] no *.tsv.gz shards under {p}")
        for shard in shards:
            with gzip.open(shard, "rt", newline="") as f_in:
                yield from _csv.DictReader(f_in, delimiter="\t")
    else:
        with gzip.open(p, "rt", newline="") as f_in:
            yield from _csv.DictReader(f_in, delimiter="\t")


def _perm_row_passes_change_side(ct, cb, change_side):
    """Same semantics as the observed-side change_side filter in main(): 'top'/
    'bottom' each admit rows tagged with that direction OR 'both' (ct==1 alone, or
    cb==1 alone, already covers the 'both' case since a both-direction position has
    ct==1 AND cb==1); the default 'both' admits any row with directional signal."""
    if change_side == 'top':
        return ct == 1
    if change_side == 'bottom':
        return cb == 1
    return ct == 1 or cb == 1


def _read_authoritative_cycles(gene_cycle_scores_path):
    """The exact cycle enumeration, read from gene_cycle_scores.tsv rather than
    inferred from perm_pos_detail. gene_wrapper.py's _flush() writes one row per
    (gene, cycle) for EVERY cycle in cycle_tags unconditionally -- "a permuted
    labeling that produced no hit anywhere in this gene is a structural zero, not
    missing data" -- so every gene's rows alone already enumerate every cycle.
    Reading just the first gene's block is enough, but reading the whole file's
    distinct `cycle` column is cheap (this file is genes x cycles rows, not the
    full per-position detail) and does not depend on row ordering.
    """
    import csv as _csv
    cycles = set()
    with open(gene_cycle_scores_path, newline="") as f:
        reader = _csv.DictReader(f, delimiter="\t")
        for row in reader:
            cycles.add(row.get('cycle'))
    return cycles


def run_permulation_null(detail_path, gene_to_id, n_genes, actual_counts, change_side,
                          gene_cycle_scores_path=None):
    """Build the permulation null's (sum, sum_sq, count_above) accumulators by
    streaming perm_pos_detail once. Returns a `results` list with the same shape
    process_chunk() produces, so main()'s existing aggregation/output-writing tail
    (which just does np.sum([r[cat][...] for r in results]) over that list) needs
    no changes to consume it.

    Per (gene, cycle) category counts are accumulated in a dict rather than a
    dense genes x cycles x categories array: most (gene, cycle) pairs have zero
    detections (Pass A in gene_wrapper.py only ever emits a row for an actually
    detected position), so the dict is sized by the real detail-row volume, not
    by the full cross product.

    `gene_cycle_scores_path`, when given, provides the AUTHORITATIVE cycle count
    (see _read_authoritative_cycles). Without it, N falls back to the union of
    `cycle` values that happen to appear in perm_pos_detail, which silently
    understates N by one for any cycle with zero detections genome-wide across
    every gene and every scheme -- an edge case at production scale, but not a
    zero-probability one, and one the pipeline's own null-generation code
    explicitly guards against elsewhere (see the docstring above). Always pass
    gene_cycle_scores_path when it is available.
    """
    per_gene_cycle_counts = {}   # gene_id -> {cycle_str: np.ndarray(len(_PERM_CATS), int64)}
    all_cycles = set()
    n_gene_missing = 0
    n_group_unrecognised = 0

    for row in _iter_perm_detail_rows(detail_path):
        gid = gene_to_id.get(row.get('Gene'))
        if gid is None:
            n_gene_missing += 1
            continue
        cyc = row.get('cycle')
        all_cycles.add(cyc)

        try:
            ct = int(row.get('ct', 0) or 0)
            cb = int(row.get('cb', 0) or 0)
        except (TypeError, ValueError):
            ct, cb = 0, 0
        if not _perm_row_passes_change_side(ct, cb, change_side):
            continue

        cat = _PERM_CAT_OF_GROUP.get((row.get('caap_group') or '').strip().upper())
        if cat is None:
            n_group_unrecognised += 1
            continue

        cyc_map = per_gene_cycle_counts.setdefault(gid, {})
        vec = cyc_map.get(cyc)
        if vec is None:
            vec = np.zeros(len(_PERM_CATS), dtype=np.int64)
            cyc_map[cyc] = vec
        vec[_PERM_CATS.index(cat)] += 1

    if n_gene_missing:
        logging.warning(
            f"[permulation] {n_gene_missing} perm_pos_detail rows referenced genes "
            f"outside the current gene universe (dropped -- e.g. removed by post-processing "
            f"after the null was generated).")
    if n_group_unrecognised:
        logging.warning(
            f"[permulation] {n_group_unrecognised} perm_pos_detail rows had an "
            f"unrecognised caap_group (dropped).")

    if gene_cycle_scores_path:
        authoritative_cycles = _read_authoritative_cycles(gene_cycle_scores_path)
        n_missing_from_detail = len(authoritative_cycles - all_cycles)
        if n_missing_from_detail:
            logging.warning(
                f"[permulation] {n_missing_from_detail} cycle(s) had zero detections "
                f"genome-wide across every gene and scheme in perm_pos_detail -- caught "
                f"exactly because gene_cycle_scores.tsv was supplied; without it these "
                f"would have silently understated N.")
        n_total = len(authoritative_cycles)
    else:
        logging.warning(
            "[permulation] no --gene-cycle-scores supplied: N is inferred from the union "
            "of cycles seen in perm_pos_detail, which silently understates N by one for "
            "any cycle with zero detections genome-wide across every gene and scheme. "
            "Pass --gene-cycle-scores (gene_cycle_scores.tsv from the same "
            "CAAS_PERMULATION run) for an exact count.")
        n_total = len(all_cycles)

    if n_total == 0:
        raise RuntimeError(
            "[permulation] found zero cycles -- check --perm-pos-detail/--gene-cycle-scores")
    logging.info(f"[permulation] {n_total} cycles; "
                 f"{len(per_gene_cycle_counts)} genes have at least one detection")

    S = {cat: [np.zeros(n_genes, dtype=np.int64),
               np.zeros(n_genes, dtype=np.float64),
               np.zeros(n_genes, dtype=np.int64)]
         for cat in _PERM_CATS}

    # Genes with zero detections in every cycle: the null count is 0 in every
    # cycle, so count_above = n_total when the observed count is also <= 0
    # (0 >= obs holds for every cycle), else 0. sum/sum_sq stay at 0.
    for cat in _PERM_CATS:
        obs_vec = actual_counts[cat]
        S[cat][2][:] = np.where(obs_vec <= 0, n_total, 0)

    for gid, cyc_map in per_gene_cycle_counts.items():
        mat = np.array(list(cyc_map.values()), dtype=np.int64)  # (n_present_cycles, len(_PERM_CATS))
        n_present = mat.shape[0]
        n_absent = n_total - n_present
        for cat_idx, cat in enumerate(_PERM_CATS):
            col = mat[:, cat_idx]
            obs = actual_counts[cat][gid]
            S[cat][0][gid] = col.sum()
            S[cat][1][gid] = float((col.astype(np.float64) ** 2).sum())
            above_present = int((col >= obs).sum())
            above_absent = n_absent if obs <= 0 else 0
            S[cat][2][gid] = above_present + above_absent

    def _pack(Sv):
        return {'sum': Sv[0], 'sum_sq': Sv[1], 'count_above': Sv[2]}

    return [{'chunk_idx': 0, 'n_rands': n_total, **{cat: _pack(S[cat]) for cat in S}}]


def _write_empty_outputs_and_exit(args):
    """Write empty-but-valid outputs when no rows are available for randomization."""
    categories = ['us', 'gs1', 'gs2', 'gs3', 'gs4']

    base_dir = os.path.dirname(args.output_prefix) or '.'
    os.makedirs(base_dir, exist_ok=True)

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
        ])

        out_cat = f"{args.output_prefix}_{cat}_aggregated_results.csv"
        df_cat.to_csv(out_cat, index=False, compression='gzip' if args.compress else None)
        logging.info(f"Results ({cat}): {out_cat} [empty]")

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
    import os as _os
    if _os.path.getsize(args.caas_csv) == 0:
        logging.warning("CAAS input file is empty — no positions passed the filter; proceeding with all-null CAAS join")
        caas_raw = pd.DataFrame()
    else:
        caas_raw  = pd.read_csv(args.caas_csv, sep=None, engine='python')

    # Keep only needed columns to save memory
    required_cols = ['gene', 'msa_pos',
                     'change_side', 'tag', 'convergence_type', 'iscaap', 'caap_group',
                     'is_conserved_meta']

    # Normalise CAAS columns (handles global_meta_caas.tsv or original schema)
    # If caas_raw is empty, produce a minimal skeleton so the left-join keys exist
    if caas_raw.empty:
        caas_df = pd.DataFrame(columns=required_cols)
    else:
        caas_df = _remap_caas_df(caas_raw)

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
    cons_idx_arr  = merged_df['cons_idx'].values.astype(np.float32)

    # ── Eligible pool ────────────────────────────────────────────────────────
    # A row is eligible to be drawn by the null only if it satisfies BOTH
    # conditions an observed CAAS satisfies:
    #   1. not `masked` — the alignment column is not a gap;
    #   2. present in background.output — CAAStools actually tested it, in a gene
    #      that survived post-processing.
    # Condition 1 alone admits columns that were never testable, which the null
    # would then scatter CAAS into; see the module docstring for how that
    # miscalibrates each gene by its own tested/ungapped ratio.
    #
    # `--background-positions` is optional at the CLI so the module still runs
    # when only the global table is available; that path logs a warning because
    # the resulting pool is condition-1-only.
    masked = merged_df['masked'].values.astype(bool)
    if args.randomization_type == 'permulation':
        # No eligible-pool sampling for this type — the null's positions are the
        # ACTUAL detections replayed from perm_pos_detail (see run_permulation_null),
        # so --background-positions/--bg-caas are accepted but unused.
        logging.info(
            "[permulation] randomization type reads its null from --perm-pos-detail "
            "and does not draw from an eligible pool.")
    elif getattr(args, 'background_positions', None):
        universe_genes = None
        if getattr(args, 'bg_caas', None) and os.path.exists(args.bg_caas):
            universe_genes = load_universe_genes(args.bg_caas)
            logging.info(f"Universe (cleaned background) genes: {len(universe_genes)}")
        tested = load_tested_positions(args.background_positions, universe_genes)
        logging.info(f"Tested positions in eligible pool: {len(tested)}")

        gene_arr = merged_df['gene'].to_numpy()
        mpos_arr = merged_df['msa_pos'].to_numpy()
        in_tested = np.fromiter(
            ((g, int(p)) in tested for g, p in zip(gene_arr, mpos_arr)),
            dtype=bool, count=len(merged_df))
        masked = masked | ~in_tested

        n_elig = int((~masked).sum())
        logging.info(
            f"Eligible pool: {n_elig} positions "
            f"({n_elig / max(len(merged_df), 1):.1%} of {len(merged_df)} rows; "
            f"{int((~merged_df['masked'].values.astype(bool)).sum())} are ungapped, "
            f"of which the tested subset is kept)")

        # Every observed CAAS must live inside the pool. A CAAS outside it is one
        # the null can never reproduce, which both understates that gene's null
        # counts and pushes the per-decile draw toward the replace=True branch
        # below (when a decile holds fewer eligible positions than observed CAAS).
        caas_mask_chk = np.asarray(
            merged_df['iscaas'].replace({'TRUE': True, 'FALSE': False})
                .astype('boolean').fillna(False), dtype=bool)
        n_unreachable = int((caas_mask_chk & masked).sum())
        if n_unreachable:
            logging.warning(
                f"{n_unreachable} observed CAAS positions are OUTSIDE the eligible pool "
                f"and can never be drawn by the null — p-values for their genes will be "
                f"anti-conservative. Check that background.output matches this run's discovery.")
        else:
            logging.info("All observed CAAS positions are inside the eligible pool.")
    else:
        logging.warning(
            "No --background-positions supplied: the eligible pool is every ungapped "
            "alignment column, which includes positions CAAStools never tested and so "
            "cannot be the source of an observed CAAS. The pool is roughly 4x larger than "
            "the testable set, and each gene's null expectation is scaled by its own "
            "tested/ungapped ratio rather than by its genuine opportunity to carry a CAAS. "
            "Pass caastools background.output (with --bg-caas) to define the pool properly.")

    if args.randomization_type != 'permulation':
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

    # Base pool candidates by group + event side
    _gu = merged_df['caap_group'].fillna('').astype(str).str.strip().str.upper()
    pool_groups = {'US', 'GS1', 'GS2', 'GS3', 'GS4'}
    pool_group_mask = _gu.isin(pool_groups)
    if args.change_side in ('top', 'bottom'):
        # Include positions with change_side == target direction OR "both"
        caas_filter = pool_group_mask & merged_df['change_side'].isin([args.change_side, 'both'])
        logging.info(f"Direction filter '{args.change_side}': {caas_filter.sum()} positions retained")
    else:
        # Default: all positions that have any directional signal
        caas_filter = pool_group_mask & merged_df['change_side'].fillna('').ne('none')

    decile_bins = None
    if args.randomization_type == 'cons_decile':
        decile_bins = (
            _compute_bins_from_series(merged_df['cons_idx'])
            if args.decile_bins is None
            else np.array([float(x) for x in args.decile_bins.split(',')])
        )

    # Apply row validity and keep only per-group output pools.
    # Counts are kept strictly per group. There is deliberately no pooled
    # "all groups" count: a single physical position can satisfy several grouping
    # schemes, so summing the five counts would count that position once per
    # scheme. Cross-group evidence is instead combined at the p-value level
    # downstream (Cauchy Combination Test), which needs no such sum and tolerates
    # the positive dependence that the shared positions induce.
    pool_mask = caas_filter

    us_mask  = pool_mask & _gu.eq('US')
    gs1_mask = pool_mask & _gu.eq('GS1')
    gs2_mask = pool_mask & _gu.eq('GS2')
    gs3_mask = pool_mask & _gu.eq('GS3')
    gs4_mask = pool_mask & _gu.eq('GS4')

    logging.info(
        f"Category sizes: us={us_mask.sum()}, gs1={gs1_mask.sum()}, "
        f"gs2={gs2_mask.sum()}, gs3={gs3_mask.sum()}, gs4={gs4_mask.sum()}"
    )

    def _bc(mask):
        return np.bincount(np.asarray(genes_int[np.asarray(mask, dtype=bool)], dtype=np.int32), minlength=n_genes)

    # Each group gets an independent draw, which calibrates its marginal p-value.
    # Cross-group dependence is handled at combination time (CCT), not here.
    actual_counts_masks = dict([
        ('us',  us_mask),
        ('gs1', gs1_mask),
        ('gs2', gs2_mask),
        ('gs3', gs3_mask),
        ('gs4', gs4_mask),
    ])
    actual_counts = {cat: _bc(mask) for cat, mask in actual_counts_masks.items()}

    if args.randomization_type == 'permulation':
        if not getattr(args, 'perm_pos_detail', None):
            raise ValueError("--perm-pos-detail is required for --randomization-type permulation")
        logging.info(f"[permulation] reading null from {args.perm_pos_detail}")
        results = run_permulation_null(
            args.perm_pos_detail, gene_to_id, n_genes, actual_counts, args.change_side,
            gene_cycle_scores_path=getattr(args, 'gene_cycle_scores', None),
        )
    else:
        _, caas_data = _build_caas_payload(merged_df, args.randomization_type, decile_bins, pool_mask)

        def _key_of(row):
            if args.randomization_type == 'naive':
                return 'global'
            if decile_bins is None:
                raise ValueError("decile_bins is required for cons_decile randomization")
            k = int(np.digitize([row['cons_idx']], bins=decile_bins[:-1], right=False)[0])
            return np.clip(k - 1, 0, len(decile_bins) - 2)

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
                    caas_data if args.export_individual_rand else None,
                    position_to_tag if args.export_individual_rand else None,
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
        base_dir = os.path.dirname(args.output_prefix) or '.'
        os.makedirs(base_dir, exist_ok=True)
        out_cat = f"{args.output_prefix}_{cat}_aggregated_results.csv"
        df_cat.to_csv(out_cat, index=False,
                      compression='gzip' if args.compress else None)
        logging.info(f"Results ({cat}): {out_cat}")

    # Decile window info (diagnostic: how many eligible positions land in
    # each decile, globally and per gene — lets a run be sanity-checked
    # against the replace=True fallback in process_chunk without rerunning).
    if args.randomization_type == 'cons_decile':
        df_elig = merged_df[merged_df['masked'] == False].copy()
        bins = decile_bins if decile_bins is not None else _compute_bins_from_series(
            merged_df.loc[pool_mask, 'cons_idx']
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
    parser = argparse.ArgumentParser(description='CT_ACCUMULATION: permutation randomization (full pool + per-group)')
    parser.add_argument('--global_csv',   required=True)
    parser.add_argument('--caas_csv',     required=True)
    parser.add_argument('--background-positions', dest='background_positions', default=None,
                        help='caastools background.output (gene<TAB>tested msa positions). '
                             'Defines the eligible null pool; without it the pool falls back '
                             'to ungapped alignment columns, which is ~4x too large.')
    parser.add_argument('--bg-caas', dest='bg_caas', default=None,
                        help='cleaned_background_main.txt gene list; restricts the eligible '
                             'pool to post-processing-surviving genes.')
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--randomization-type', choices=['naive', 'cons_decile', 'permulation'], required=True)
    parser.add_argument('--perm-pos-detail', dest='perm_pos_detail', default=None,
                        help='perm_pos_detail/ shard dir (or legacy perm_pos_detail.tsv.gz) from '
                             'a prior CAAS_PERMULATION run. Required iff --randomization-type '
                             'permulation: the null-anchored accumulation type replays the SAME '
                             'permulation cycles CAAS uses elsewhere, rather than drawing from an '
                             'eligible pool.')
    parser.add_argument('--gene-cycle-scores', dest='gene_cycle_scores', default=None,
                        help='gene_cycle_scores.tsv from the SAME CAAS_PERMULATION run as '
                             '--perm-pos-detail. Optional but strongly recommended for '
                             '--randomization-type permulation: gives the exact cycle count '
                             '(N) instead of inferring it from which cycles happen to appear '
                             'in perm_pos_detail, which understates N for any cycle with zero '
                             'detections genome-wide.')
    parser.add_argument('--n-randomizations', type=int, default=100000)
    parser.add_argument('--workers',      type=int, default=None)
    parser.add_argument('--compress',     action='store_true')
    parser.add_argument('--export-individual-rand', action='store_true')
    parser.add_argument('--decile-bins',  type=str, default=None)
    parser.add_argument('--global-seed',  type=int, default=None)
    parser.add_argument('--precompute-masks', dest='precompute_masks', action='store_true')
    parser.add_argument('--no-precompute-masks', dest='precompute_masks', action='store_false')
    parser.set_defaults(precompute_masks=True)
    parser.add_argument('--change-side', dest='change_side', default='both',
                        choices=['top', 'bottom', 'both'],
                        help='Restrict the CAAS position pool to this phenotype direction. '
                             '"top" and "bottom" each include positions with change_side=="both". '
                             '"both" (default) retains all non-none positions (original behaviour).')
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')
    main(args)
