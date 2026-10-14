#!/usr/bin/env python3
# randomize.py — Per-gene empirical null for CAAS accumulation.
# PhyloPhere | subworkflows/CT_ACCUMULATION/local/src/randomization/

"""
Randomization: tests whether a gene accumulates more CAAS than a null predicts.

Writes <prefix>_aggregated_results.csv, plus conservation-decile parquet tables for 'cons_decile'.

Test design
-----------
Two randomization families, selected by --randomization-type:

'naive' / 'cons_decile': an occupancy null. For each group, the observed number of CAAS
positions is redrawn uniformly (without replacement, with replacement when the pool is
smaller than the draw) from the eligible position pool, and the per-gene counts are
recounted. 'naive' draws from the whole pool; 'cons_decile' keeps the number of CAAS per
conservation decile (deciles of the conservation index of the whole position table, applied
within the pool). Because the pool is position-level, gene length (strictly, each
gene's number of eligible positions) is controlled for by construction.

'permulation': an excess null (see run_permulation_null). The randomized CAAS set of a
gene at cycle i is that gene's detections in the filtered CAAS permulation null
(perm_pos_cycle_caas.tsv.gz, the table SCORING, POSENRICH and FCS also use), not a
synthetic draw from a pool, and the per-gene count is free in each cycle rather than
pinned to the observed count. N is the replay cycle roster when given. It holds the
phenotype-tree confound but is not conservation-decile matched.

In every case the empirical p-value per gene is
    p = (#{null count >= observed count} + 1) / (N + 1)
so a gene is significant when it accumulates more CAAS than its null predicts.

Eligible pool
-------------
The pool is the set of positions an observed CAAS could actually have come from: the
positions CAAStools tested (background.output), restricted to the genes that survive
post-processing (cleaned_background_main.txt). This mirrors build_background() of
posenrich_enrich.py. A position that failed the gap, missingness or pattern filters can
never be a CAAS, and neither can any position of a gene that post-processing removed.

When the filtered permulation null (--null-cycle-caas) is given, 'naive' and 'cons_decile'
draw from the positional universe U shared with SCORING and POSENRICH: the tested positions
that score caas_score > 0 in at least one null cycle (any side), plus the observed CAAS
positions. A position the null never scores cannot produce a CAAS under the permulation
labelings, so it is not offered to the draws. Without the null table the pool is the tested
positions as described above, and the log says which pool was used. Keys (deciles, or the
whole pool for 'naive') that hold fewer pool positions than observed CAAS are drawn with
replacement and reported in the log and in the decile tables.

The pool must be proportional to each gene's genuine opportunity to carry a CAAS, because
a gene's null expectation is proportional to its share of the pool and the per-gene
p-values are only comparable under that condition. Widening the pool to every ungapped
alignment column breaks the proportionality gene by gene: poorly covered genes receive a
too-high expectation (conservative) and well covered genes a too-low one
(anti-conservative), which concentrates false positives in the best-covered genes.

Each gene's empirical p-value is derived from the null distribution of its CAAS count and adjusted directly by
Benjamini-Hochberg (BH FDR) downstream in accum_gene_lists.nf and scoring_compute.R.

No significance, convergence or divergence category families, FDR or gene lists are
produced here.

Called by:  CT_ACCUMULATION Nextflow process (ctacc_run.nf CT_ACCUMULATION_RANDOMIZE → main.py --tool randomize)
Inputs:     --global_csv (<prefix>_global.csv of concatenate.py), --caas_csv (filtered_discovery.tsv),
            --background-positions and --bg-caas (eligible pool), --null-cycle-caas (universe U of
            'naive' and 'cons_decile'; the null of 'permulation') and --cycle-roster (permulation only)
Outputs:    <output-prefix>_aggregated_results.csv (GeneID, Gene, ActualCount, RandsAbove,
            RandsBelow, PValueEmpirical, MeanSim, SD, NumRands);
            for cons_decile also <output-prefix>_cons_decile_global.parquet (decile, NumPositions,
            NumPool, NumObserved, Replace), <output-prefix>_cons_decile_per_gene.parquet (ungapped
            positions per gene and decile) and <output-prefix>_cons_decile_pool_per_gene.parquet
            (draw-pool positions per gene and decile)
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import os

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from collections import defaultdict
import pyarrow as pa
import pyarrow.parquet as pq
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import random
import multiprocessing.shared_memory as shm


# ── Eligible position pool ────────────────────────────────────────────────────


def load_universe_genes(path):
    """Genes that survive post-processing (cleaned_background_main.txt: one gene per line)."""
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
    positions`, one line per gene. The positions are in the same coordinate space
    as filtered_discovery's `Position`, i.e. `msa_pos` after _remap_caas_df, which
    is also the key the global/CAAS merge joins on. Do not match on the global
    table's `position` column: it is a globally unique running row index, not an
    alignment coordinate.

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


# ── Positional universe U (occupancy nulls) ───────────────────────────────────
#
# With the filtered permulation null available, the occupancy nulls draw from the same positional universe as
# SCORING and POSENRICH:
#   U = (null-positive positions inside the tested pool) union (observed CAAS positions)
# where null-positive means caas_score > 0 in at least one cycle and side of perm_pos_cycle_caas.tsv.gz and the
# tested pool is background.output restricted to the genes of the cleaned background. An observed position stays
# in U whether or not the null scored it, so every observed count is reachable by the draws.


def load_null_positive_positions(null_path):
    """Distinct (Gene, Position) with caas_score > 0 in any cycle and side of the filtered null table.

    The table is read in chunks and reduced to distinct pairs chunk by chunk, so memory follows the number of
    distinct null-positive positions, not the number of (position, cycle) rows. A file without content gives an
    empty table.
    """
    cols = ['Gene', 'Position', 'caas_score']
    empty = pd.DataFrame({'Gene': pd.Series(dtype=str), 'Position': pd.Series(dtype=np.int64)})
    if os.path.getsize(null_path) == 0:
        return empty
    header = pd.read_csv(null_path, sep='\t', nrows=0).columns
    missing = [c for c in cols if c not in header]
    if missing:
        raise ValueError(f"{null_path} lacks column(s) {missing}: regenerate the CAAS permulation null (CAAS_CORE_MERGE).")
    parts = [empty]
    for ch in pd.read_csv(null_path, sep='\t', usecols=cols, chunksize=2_000_000,
                          dtype={'Gene': str, 'Position': str, 'caas_score': str}):
        score = pd.to_numeric(ch['caas_score'], errors='coerce')
        ch = ch.loc[score > 0, ['Gene', 'Position']]
        ch['Position'] = pd.to_numeric(ch['Position'], errors='coerce')
        parts.append(ch.dropna().astype({'Position': np.int64}).drop_duplicates())
    return pd.concat(parts, ignore_index=True).drop_duplicates()


def _pair_keys(gene_ids, positions):
    """One int64 key per (gene id, position): gene_id * 2^32 + position."""
    return np.asarray(gene_ids, dtype=np.int64) * (1 << 32) + np.asarray(positions, dtype=np.int64)


def in_null_positive(gene_ids, msa_pos, null_pos, gene_to_id):
    """Row-level flag: the (gene, msa_pos) of each row is a null-positive position.

    gene_ids are the integer gene ids of the rows (gene_to_id maps gene name to id); null_pos is the output of
    load_null_positive_positions. Genes of the null outside gene_to_id cannot match any row and are ignored.
    """
    ng = null_pos['Gene'].map(gene_to_id)
    ok = ng.notna().to_numpy()
    null_keys = _pair_keys(ng[ok].to_numpy(), null_pos['Position'].to_numpy()[ok])
    return np.isin(_pair_keys(gene_ids, msa_pos), null_keys)


def eligible_rows(positions, masked):
    """Indices of the draw pool: one row per distinct position that is not masked.

    `positions` is the globally unique position id, so the repeats of a position in the merged table (one per
    side or hypothesis) count once; a position that is a CAAS on both sides would otherwise sit twice in
    the pool and the null would re-land preferentially on genes that are already CAAS-dense.
    """
    _, first_idx = np.unique(positions, return_index=True)
    distinct = np.zeros(len(positions), dtype=bool)
    distinct[first_idx] = True
    return np.flatnonzero(~np.asarray(masked, dtype=bool) & distinct)


def decile_keys(cons_idx, decile_bins):
    """0-based conservation decile of each value; the extremes are clipped into the first and last decile."""
    k = np.digitize(cons_idx, bins=decile_bins[:-1], right=False)
    return np.clip(k - 1, 0, len(decile_bins) - 2)


def restrict_pool_to_universe(masked, gene_ids, msa_pos, observed, null_pos, gene_to_id):
    """Masked flag of the universe U: (unmasked and null-positive) or observed.

    `masked` is True where the row is outside the tested pool (gap, untested or removed gene), `observed` flags the
    rows of observed CAAS positions (either side). Returns the new masked array (True = outside U).
    """
    in_null = in_null_positive(gene_ids, msa_pos, null_pos, gene_to_id)
    return ~((~np.asarray(masked, dtype=bool) & in_null) | np.asarray(observed, dtype=bool))


# ── Column remapping ──────────────────────────────────────────────────────────


def _remap_caas_df(df):
    """Normalize a filtered_discovery.tsv DataFrame to the internal column schema.

    Renames Gene and Position to gene and msa_pos, fills the optional columns that are
    absent (tag from tag_support, side from change_side or 'both', is_conserved_meta as
    False). Other columns pass through unchanged.
    """
    # Structural keys only; the other columns already carry their lowercase names
    df = df.rename(columns={'Gene': 'gene', 'Position': 'msa_pos', 'position': 'msa_pos'})

    _bool = lambda x: str(x).strip().lower() in {'true', 't', '1', 'yes', 'y'}

    # Conserved-state column to bool (optional)
    if 'is_conserved_meta' in df.columns:
        df['is_conserved_meta'] = df['is_conserved_meta'].map(_bool)
    else:
        df['is_conserved_meta'] = False

    # tag: fallback to tag_support if tag not present
    if 'tag' not in df.columns:
        if 'tag_support' in df.columns:
            df['tag'] = df['tag_support']
        else:
            df['tag'] = ''

    # side: fallback to change_side or 'both' if not present
    if 'side' not in df.columns:
        if 'change_side' in df.columns:
            df['side'] = df['change_side']
        else:
            df['side'] = 'both'

    logging.info("CAAS CSV: remapped filtered_discovery.tsv columns to internal schema")
    return df


# ── Worker class ──────────────────────────────────────────────────────────────


class RandomizationWorker:
    """Per-process state of the occupancy null ('naive' and 'cons_decile').

    Maps the shared-memory arrays of the position table once per worker process and
    draws chunks of randomizations from them (process_chunk). Each draw redraws the
    observed CAAS count of every group, key (global, or conservation decile) by key.
    """

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
        actual_counts,          # dict: group name -> np.ndarray(n_genes, int64), observed CAAS per gene
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
        # One draw slot per physical position, not one per side row (see eligible_rows).
        all_eligible = eligible_rows(self.positions, self.masked)
        if self.randomization_type == 'naive':
            self.eligible_by_key['global'] = all_eligible
            print(f"Total eligible positions for 'naive' randomization: {len(all_eligible)}")
            return
        if self.randomization_type == 'cons_decile':
            assert self.decile_bins is not None
            dec = decile_keys(self.cons_idx, self.decile_bins)
            in_pool = np.zeros(self.n_rows, dtype=bool)
            in_pool[all_eligible] = True
            for d in range(len(self.decile_bins) - 1):
                self.eligible_by_key[d] = self.row_indices[(dec == d) & in_pool]
            return
        raise ValueError(f"Unknown randomization_type: {self.randomization_type}")

    def process_chunk(self, chunk_size, chunk_idx):
        """Run chunk_size randomizations and return the per-gene accumulators.

        The seed derives from global_seed and chunk_idx, so a chunk is reproducible
        and chunks do not share a stream. Returns {'chunk_idx', 'n_rands', <group>:
        {'sum', 'sum_sq', 'count_above'}} with one value per gene: the sum and sum of
        squares of the null counts, and the number of draws whose count reaches the
        observed count.
        """
        logging.info(f"Starting chunk {chunk_idx} with {chunk_size} randomizations in PID {os.getpid()}")
        seed = None
        if self.global_seed is not None:
            seed = ((self.global_seed or 0) + 9973 * chunk_idx) & 0x7FFFFFFF
        rng = np.random.default_rng(seed)

        # Accumulators per group: [sum, sum_sq, count_above], updated in place.
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

        # One draw buffer per group, filled slice by slice (one slice per key)
        cat_total_sizes = {cat: self.actual_counts[cat].sum() for cat in self.actual_counts}
        idx_buffers = {cat: np.empty(cat_total_sizes[cat], dtype=np.int32) for cat in self.actual_counts}

        offsets = {cat: 0 for cat in self.actual_counts}
        decile_plans = []

        for key, cinfo in self.caas_by_key_counts.items():
            elig = self.eligible_by_key.get(key)
            if elig is None or elig.size == 0:
                continue

            # (start, end) slice of each group's buffer for this key; keys without CAAS are skipped
            cat_slices = {}
            for cat in self.actual_counts:
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
                # Each category draws its positions from this key's eligible set.
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


# ── Worker helpers (module level, for pickling) ───────────────────────────────


def init_worker(
    randomization_type, decile_bins, export_individual, output_dir, global_seed,
    precompute_masks, n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
    caas_data, position_to_tag,
    actual_counts,              # dict: group name -> np.ndarray
):
    """ProcessPoolExecutor initializer: build the process-wide RandomizationWorker."""
    global worker
    worker = RandomizationWorker(
        randomization_type, decile_bins, export_individual, output_dir, global_seed,
        precompute_masks, n_rows, n_genes, shm_names, shapes, dtypes, extra_key_sizes,
        caas_data, position_to_tag,
        actual_counts,
    )

def process_wrapper(ch):
    """Run one (chunk_size, chunk_idx) task on the process-wide worker."""
    size, idx = ch
    return worker.process_chunk(size, idx)


# ── Utility functions ─────────────────────────────────────────────────────────


def _compute_bins_from_series(series):
    """Decile edges (0, 10, ..., 100th percentile) of the non-missing values."""
    return np.percentile(series.dropna(), np.arange(0, 101, 10))


def _build_caas_payload(merged_df, randomization_type, decile_bins, pool_mask):
    """Observed CAAS rows (those in pool_mask) with their draw key: 'global', or the decile.

    Returns (unique keys, [{'position', 'key'}]). Other randomization types yield nothing.
    """
    caas_data = []
    for _, row in merged_df.loc[pool_mask].iterrows():
        if randomization_type == 'naive':
            key = 'global'
        elif randomization_type == 'cons_decile':
            key = int(decile_keys(np.float32(row['cons_idx']), decile_bins))
        else:
            continue
        caas_data.append({'position': row['position'], 'key': key})
    unique_keys = list({item['key'] for item in caas_data})
    return unique_keys, caas_data


# ── Permulation null ──────────────────────────────────────────────────────────
#
# Unlike 'naive' and 'cons_decile', this type draws no synthetic positions from an
# eligible pool. The randomized CAAS set of a gene at cycle i is the set of positions
# the CAAS permulation detects in that gene at cycle i. It is read from
# perm_pos_cycle_caas.tsv.gz, the same filtered null table SCORING (p.emp, p.emp_fact,
# p.adj_bh), POSENRICH and FCS use: cluster trains and removed (cycle, gene)
# units are already out of it, as they are out of filtered_discovery.tsv on the observed
# side, so both sides count the same kind of position. The null keeps the phenotype-tree
# confound but is not conservation-decile matched: it answers a different question, not a
# stricter version of the same one. The per-gene CAAS count is free in each cycle (not
# fixed to the observed count), so this is an excess null, not an occupancy null.
#
# One counting rule for both sides (see _side_mask): a row is a hit when it is a
# detection (caas_score > 0 in the null, the rows of filtered_discovery on the
# observed side) whose side is the selected one; with direction 'both' a position on
# both sides is two rows, on either side of the comparison.

_PERM_CATS = ['caas']


def _side_mask(side, change_side):
    """Boolean mask of the rows whose `side` the direction selector keeps.

    'top' or 'bottom' keep that side only; 'both' keeps top and bottom rows (a position
    detected on both sides is two rows). Rows with side 'none' or missing never count.
    The observed pool (main) and the null (null_counts_by_gene_cycle) both go through
    this function.
    """
    s = pd.Series(side).astype(str).str.strip().str.lower()
    if change_side in ('top', 'bottom'):
        return s.eq(change_side)
    return s.isin(('top', 'bottom'))


def _base_cycle(tag):
    """'b_5~H3' -> 'b_5': the cycle of a per-hypothesis replay tag."""
    return str(tag).split('~', 1)[0]


def null_counts_by_gene_cycle(null_path, change_side):
    """Hits of the permulation null per (Gene, cycle), plus the set of cycles present in the file.

    Reads perm_pos_cycle_caas.tsv.gz (Gene, Position, side, cycle, caas_score, ...) and keeps the
    rows with caas_score > 0 on the selected side. A (Gene, Position, side, cycle) key counts once,
    as the observed table is deduplicated per (gene, position, side). The cycles present are
    taken before the filter, so a cycle whose rows were all dropped is still a present cycle.

    Returns (counts, cycles_present): counts is a Series indexed by (Gene, cycle) with the number
    of hits, cycles_present a set of base cycle labels.
    """
    cols = ['Gene', 'Position', 'side', 'cycle', 'caas_score']
    header = pd.read_csv(null_path, sep='\t', nrows=0).columns
    missing = [c for c in cols if c not in header]
    if missing:
        raise ValueError(f"[permulation] {null_path} lacks column(s) {missing}: regenerate the CAAS permulation null "
                         f"(CAAS_CORE_MERGE).")
    df = pd.read_csv(null_path, sep='\t', usecols=cols,
                     dtype={'Gene': 'category', 'side': 'category', 'cycle': 'category'})
    if (df['cycle'].astype(str).str.contains('~', regex=False)).any():
        logging.warning("[permulation] the null carries '<base>~H*' hypothesis-replay tags; they are collapsed to "
                        "the base cycle (domain-pooled, not independent draws).")
    df['cycle'] = df['cycle'].astype(str).map(_base_cycle)
    cycles_present = set(df['cycle'].unique())
    score = pd.to_numeric(df['caas_score'], errors='coerce')
    df = df[(score > 0) & _side_mask(df['side'], change_side).to_numpy()]
    df = df.drop_duplicates(['Gene', 'Position', 'side', 'cycle'])
    counts = df.groupby(['Gene', 'cycle'], observed=True).size()
    return counts, cycles_present


def read_cycle_roster(roster_path):
    """Cycle labels of cycle_roster.txt (one per line), as a set; an empty set without a path."""
    if not roster_path:
        return set()
    with open(roster_path) as fh:
        return {_base_cycle(line.strip()) for line in fh if line.strip()}


def resolve_n_cycles(roster, cycles_present):
    """N of the permulation null: the replay roster when it lists every cycle present, else the cycles present.

    The roster counts the cycles that detected nothing; the cycles present in the null file cannot,
    so falling back to them can only raise the p-values. The choice is logged either way.
    """
    if roster:
        stray = cycles_present - roster
        if not stray:
            logging.info(f"[permulation] N = {len(roster)} cycles (replay roster; {len(cycles_present)} of them "
                         f"present in the null file)")
            return len(roster)
        logging.warning(f"[permulation] the roster does not list {len(stray)} cycle(s) of the null file "
                        f"(e.g. {sorted(stray)[:3]}): ignoring it, N = cycles present ({len(cycles_present)}).")
    else:
        logging.warning(f"[permulation] no cycle roster supplied: N = cycles present in the null file "
                        f"({len(cycles_present)}), which understates N by the cycles with no detection.")
    return len(cycles_present)


def run_permulation_null(null_counts, n_total, gene_to_id, n_genes, actual_counts):
    """Build the permulation null's (sum, sum_sq, count_above) accumulators from the per (gene, cycle) hit counts.

    Returns a `results` list shaped like the output of process_chunk(), so the
    aggregation and output tail of main() consumes it unchanged.

    `null_counts` is the Series of null_counts_by_gene_cycle (hits per (Gene, cycle), only the
    pairs with at least one hit); `n_total` is N. A (gene, cycle) pair without a hit counts 0 in
    that cycle, so a gene's null count in the N - n_present remaining cycles is 0.
    """
    per_gene_counts = {}   # gene_id -> int64 array of the positive counts of the cycles with a hit
    n_gene_missing = 0
    for gene, col in null_counts.groupby(level='Gene', observed=True):
        gid = gene_to_id.get(gene)
        if gid is None:
            n_gene_missing += 1
            continue
        per_gene_counts[gid] = col.to_numpy(dtype=np.int64)
    if n_gene_missing:
        logging.warning(
            f"[permulation] {n_gene_missing} genes of the null are outside the current gene universe (dropped, "
            f"e.g. removed by post-processing after the null was generated).")

    if n_total == 0:
        raise RuntimeError("[permulation] found zero cycles: check --null-cycle-caas / --cycle-roster")
    logging.info(f"[permulation] {n_total} cycles; {len(per_gene_counts)} genes have at least one detection")

    S = {cat: [np.zeros(n_genes, dtype=np.int64),
               np.zeros(n_genes, dtype=np.float64),
               np.zeros(n_genes, dtype=np.int64)]
         for cat in _PERM_CATS}

    # Genes with no detection in any cycle: the null count is 0 in every cycle, so
    # count_above = n_total when the observed count is also <= 0 (0 >= obs holds in every
    # cycle), else 0. sum and sum_sq stay at 0.
    for cat in _PERM_CATS:
        S[cat][2][:] = np.where(actual_counts[cat] <= 0, n_total, 0)

    for gid, col in per_gene_counts.items():
        n_absent = n_total - col.size
        for cat in _PERM_CATS:
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
    """Write the per-group result files with a header and no rows (and empty decile tables).

    Used when the global/CAAS join is empty, so downstream steps still find their inputs.
    """
    categories = ['caas']

    base_dir = os.path.dirname(args.output_prefix) or '.'
    os.makedirs(base_dir, exist_ok=True)

    for cat in categories:
        df_cat = pd.DataFrame(columns=[
            'GeneID',
            'Gene',
            'ActualCount',
            'RandsAbove',
            'RandsBelow',
            'PValueEmpirical',
            'MeanSim',
            'SD',
            'NumRands',
        ])

        out_cat = f"{args.output_prefix}_aggregated_results.csv"
        df_cat.to_csv(out_cat, index=False, compression='gzip' if args.compress else None)
        logging.info(f"Results ({cat}): {out_cat} [empty]")

    if args.randomization_type == 'cons_decile':
        pd.DataFrame(columns=['decile', 'NumPositions', 'NumPool', 'NumObserved', 'Replace']).to_parquet(
            f"{args.output_prefix}_cons_decile_global.parquet", index=False
        )
        for name in ('per_gene', 'pool_per_gene'):
            pd.DataFrame(columns=['gene']).to_parquet(
                f"{args.output_prefix}_cons_decile_{name}.parquet", index=False
            )

    logging.info("Randomization analysis complete (empty input).")


# ── Main ──────────────────────────────────────────────────────────────────────


def main(args):
    """Run the randomization test selected by args.randomization_type and write the results.

    Reads from args: global_csv, caas_csv, randomization_type, output_prefix, change_side,
    background_positions and bg_caas (eligible pool), null_cycle_caas (universe U of the
    occupancy nulls, null of permulation) and cycle_roster (permulation), n_randomizations, workers, global_seed, decile_bins,
    compress, export_individual_rand, precompute_masks.
    """
    logging.info("Loading data")
    # global_csv holds the positional data (cons_idx) and the per-column flags (masked, iscaas)
    global_df = pd.read_csv(args.global_csv)
    import os as _os
    if _os.path.getsize(args.caas_csv) == 0:
        logging.warning("CAAS input file is empty — no positions passed the filter; proceeding with all-null CAAS join")
        caas_raw = pd.DataFrame()
    else:
        caas_raw  = pd.read_csv(args.caas_csv, sep=None, engine='python')

    # Keep only the columns that are used, to save memory
    required_cols = ['gene', 'msa_pos',
                     'side', 'tag', 'convergence_type',
                     'is_conserved_meta']

    # Normalize the CAAS columns; an empty input becomes a skeleton so the left-join keys exist
    if caas_raw.empty:
        caas_df = pd.DataFrame(columns=required_cols)
    else:
        caas_df = _remap_caas_df(caas_raw)

    available = [c for c in required_cols if c in caas_df.columns]
    caas_df = caas_df[available]
    logging.info(f"CAAS df: {len(caas_df)} rows, columns: {available}")

    # One slot per (gene, position, side): a position detected on both
    # sides is two rows (side top and bottom), and any other repeat of the key
    # (for example one row per hypothesis) is a duplicate of the same physical unit.
    if not caas_df.empty and {'gene', 'msa_pos'}.issubset(caas_df.columns):
        n_before = len(caas_df)
        _keys = ['gene', 'msa_pos']
        if 'side' in caas_df.columns:
            _keys = _keys + ['side']
        caas_df = caas_df.drop_duplicates(_keys, keep='first')
        logging.info(
            f"CAAS df after (gene, position, side) dedup: {len(caas_df)} rows "
            f"(collapsed {n_before - len(caas_df)} hypothesis-replay duplicates)")

    # One row per global position
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
    # A row is eligible to be drawn by the null only if it satisfies both
    # conditions an observed CAAS satisfies:
    #   1. not `masked`: the alignment column is not a gap;
    #   2. present in background.output: CAAStools tested it, in a gene that
    #      survived post-processing.
    # Condition 1 alone admits columns that were never testable, where the null
    # would then scatter CAAS; the module docstring explains how that miscalibrates
    # each gene by its own tested/ungapped ratio.
    #
    # --background-positions is optional so the module still runs when only the
    # global table is available; that path logs a warning because the pool is then
    # defined by condition 1 alone.
    masked = merged_df['masked'].values.astype(bool)
    if args.randomization_type == 'permulation':
        # No eligible-pool sampling for this type: the null's positions are the
        # hits of the filtered permulation null (see run_permulation_null), so
        # --background-positions and --bg-caas are accepted but unused.
        logging.info(
            "[permulation] randomization type reads its null from --null-cycle-caas "
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

        # Count distinct positions (the draw pool is deduplicated by position in
        # _prepare_keyed_eligibles), not rows.
        _pos_arr = merged_df['position'].to_numpy()
        n_elig = int(np.unique(_pos_arr[~masked]).size)
        n_ungapped = int(np.unique(
            _pos_arr[~merged_df['masked'].values.astype(bool)]).size)
        logging.info(
            f"Eligible pool: {n_elig} distinct positions "
            f"({n_ungapped} ungapped, of which the tested subset is kept)")

    else:
        logging.warning(
            "No --background-positions supplied: the eligible pool is every ungapped "
            "alignment column, which includes positions CAAStools never tested and so "
            "cannot be the source of an observed CAAS. The pool is roughly 4x larger than "
            "the testable set, and each gene's null expectation is scaled by its own "
            "tested/ungapped ratio rather than by its genuine opportunity to carry a CAAS. "
            "Pass caastools background.output (with --bg-caas) to define the pool properly.")

    # Observed CAAS rows (top or bottom side): the rows of the position table matched by a detection
    # of filtered_discovery; a side 'none' row is not a detection
    caas_bool = np.asarray(
        merged_df['iscaas'].replace({'TRUE': True, 'FALSE': False}).astype('boolean').fillna(False), dtype=bool)
    observed_rows = merged_df['side'].notna().to_numpy() & _side_mask(merged_df['side'], 'both').to_numpy()

    pool_label = None
    if args.randomization_type != 'permulation':
        n_outside_tested = int((observed_rows & masked).sum())
        null_path = getattr(args, 'null_cycle_caas', None)
        null_pos = load_null_positive_positions(null_path) if null_path else None
        if null_pos is not None and len(null_pos) == 0:
            logging.warning(f"{null_path} has no position with caas_score > 0 (a null without cycles, for example "
                            f"--caas_full_perms 0): the universe U would be the observed positions alone, so the "
                            f"previous pool is used instead.")
            null_pos = None
        if null_pos is not None:
            masked = restrict_pool_to_universe(masked, genes_int, merged_df['msa_pos'].to_numpy(), observed_rows,
                                               null_pos, gene_to_id)
            pool_label = "U (tested pool restricted to null-positive positions, plus observed positions)"
            logging.info(f"{len(null_pos)} distinct null-positive positions in {null_path}; "
                         f"{n_outside_tested} observed rows lie outside the tested pool and stay in U")
        else:
            pool_label = ("tested pool (no null table given: positions able to produce a CAAS are not restricted)"
                          if getattr(args, 'background_positions', None) else "ungapped alignment columns")
        n_elig = int(eligible_rows(positions, masked).size)
        logging.info(f"Draw pool of '{args.randomization_type}': {pool_label}; {n_elig} distinct positions")
        # Every observed CAAS must lie inside the pool, or the null can never reproduce it
        n_unreachable = int((observed_rows & masked).sum())
        if n_unreachable:
            logging.warning(
                f"{n_unreachable} observed CAAS rows are OUTSIDE the draw pool and can never be drawn by the null: "
                f"p-values for their genes will be anti-conservative. Check that background.output matches this "
                f"run's discovery.")
        else:
            logging.info("All observed CAAS positions are inside the draw pool.")

    if args.randomization_type != 'permulation':
        genes_int_arr = genes_int.astype(np.int32)

        # Shared memory: the position table is mapped once and read by every worker
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

    # CAAS rows, filtered by side with the rule the permulation null uses
    # (_side_mask): one row per side, a position on both sides is two rows.
    caas_filter = _side_mask(merged_df['side'].fillna(''), args.change_side)
    logging.info(f"Direction filter '{args.change_side}': {caas_filter.sum()} observed CAAS rows retained")
    if not caas_df.empty and not caas_df['side'].astype(str).str.lower().isin(['top', 'bottom']).any():
        logging.warning("the CAAS table has no row with side top or bottom: no observed CAAS counts under any direction.")

    decile_bins = None
    if args.randomization_type == 'cons_decile':
        decile_bins = (
            _compute_bins_from_series(merged_df['cons_idx'])
            if args.decile_bins is None
            else np.array([float(x) for x in args.decile_bins.split(',')])
        )

    pool_mask = caas_filter
    caas_mask = pool_mask

    logging.info(f"Category size: caas={caas_mask.sum()}")
    if args.randomization_type != 'permulation':
        # The rows counted as observed belong to the pool by construction of U; a violation is a bug
        n_counted_outside = int((caas_mask.to_numpy(dtype=bool) & masked).sum())
        if n_counted_outside and pool_label.startswith('U'):
            raise AssertionError(f"{n_counted_outside} observed rows used as the observed count lie outside U")

    def _bc(mask):
        return np.bincount(np.asarray(genes_int[np.asarray(mask, dtype=bool)], dtype=np.int32), minlength=n_genes)

    actual_counts_masks = {'caas': caas_mask}
    actual_counts = {cat: _bc(mask) for cat, mask in actual_counts_masks.items()}

    if args.randomization_type == 'permulation':
        if not getattr(args, 'null_cycle_caas', None):
            raise ValueError("--null-cycle-caas is required for --randomization-type permulation")
        logging.info(f"[permulation] reading null from {args.null_cycle_caas}")
        null_counts, cycles_present = null_counts_by_gene_cycle(args.null_cycle_caas, args.change_side)
        n_total = resolve_n_cycles(read_cycle_roster(getattr(args, 'cycle_roster', None)), cycles_present)
        results = run_permulation_null(null_counts, n_total, gene_to_id, n_genes, actual_counts)
    else:
        _, caas_data = _build_caas_payload(merged_df, args.randomization_type, decile_bins, pool_mask)

        def _key_of(row):
            if args.randomization_type == 'naive':
                return 'global'
            if decile_bins is None:
                raise ValueError("decile_bins is required for cons_decile randomization")
            return int(decile_keys(np.float32(row['cons_idx']), decile_bins))

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

        # Keys whose pool holds fewer positions than the observed count: process_chunk draws them with replacement
        pool_idx = eligible_rows(positions, masked)
        if args.randomization_type == 'naive':
            pool_size_by_key = {'global': int(pool_idx.size)}
        else:
            pool_dec = decile_keys(cons_idx_arr[pool_idx], decile_bins)
            pool_size_by_key = {int(d): int((pool_dec == d).sum()) for d in range(len(decile_bins) - 1)}
        short_keys = {k: (pool_size_by_key.get(k, 0), v['n_caas']) for k, v in extra_key_sizes.items()
                      if pool_size_by_key.get(k, 0) < v.get('n_caas', 0)}
        if short_keys:
            logging.warning(
                f"{len(short_keys)} of {len(extra_key_sizes)} draw keys hold fewer pool positions than observed CAAS "
                f"rows and are drawn with replacement (key: pool, observed): {dict(sorted(short_keys.items(), key=str))}")
        else:
            logging.info("Every draw key holds at least as many pool positions as observed CAAS rows (no replacement).")

        position_to_tag = dict(zip(merged_df['position'].to_numpy(),
                                   merged_df['tag'].astype(str).fillna('').to_numpy() if 'tag' in merged_df.columns else np.array([''] * len(merged_df))))

        # Parallel plan: the randomizations are split into one chunk per worker
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
            'ActualCount':     act,
            'RandsAbove':      cnt_ge_cat,
            'RandsBelow':      total_n_rands - cnt_ge_cat,
            'PValueEmpirical': p_emp,
            'MeanSim':         mean_cat,
            'SD':              sd_cat,
            'NumRands':        total_n_rands,
        })
        base_dir = os.path.dirname(args.output_prefix) or '.'
        os.makedirs(base_dir, exist_ok=True)
        out_cat = f"{args.output_prefix}_aggregated_results.csv"
        df_cat.to_csv(out_cat, index=False,
                      compression='gzip' if args.compress else None)
        logging.info(f"Results ({cat}): {out_cat}")

    # Decile occupancy (diagnostic), globally and per gene, to check a run against the replace=True fallback of
    # process_chunk. Per decile: NumPositions (ungapped positions, as before), NumPool (distinct positions of the draw
    # pool: U when the null table was given), NumObserved (observed CAAS rows redrawn in that decile) and Replace
    # (NumPool < NumObserved). The per-gene pool occupancy is written next to the per-gene ungapped table.
    if args.randomization_type == 'cons_decile':
        # One row per physical position (as in eligible_rows): merged_df repeats each position once per side.
        df_elig = (merged_df[merged_df['masked'] == False]
                   .drop_duplicates('position').copy())
        bins = decile_bins if decile_bins is not None else _compute_bins_from_series(
            merged_df.loc[pool_mask, 'cons_idx']
        )
        df_elig['decile'] = decile_keys(df_elig['cons_idx'].to_numpy(np.float32), bins) + 1
        global_decile     = df_elig.groupby('decile').size().reset_index(name='NumPositions')
        per_gene_decile   = df_elig.groupby(['gene', 'decile']).size().reset_index(name='NumPositions')
        per_gene_pivot    = per_gene_decile.pivot(
            index='gene', columns='decile', values='NumPositions'
        ).fillna(0).reset_index()

        df_pool = merged_df.iloc[pool_idx][['gene', 'cons_idx']].copy()
        df_pool['decile'] = decile_keys(df_pool['cons_idx'].to_numpy(np.float32), bins) + 1
        pool_global = df_pool.groupby('decile').size().reset_index(name='NumPool')
        obs_global = pd.DataFrame(
            [(int(k) + 1, int(v.get('n_caas', 0))) for k, v in extra_key_sizes.items()],
            columns=['decile', 'NumObserved'])
        global_decile = (global_decile.merge(pool_global, on='decile', how='outer')
                         .merge(obs_global, on='decile', how='outer').fillna(0).sort_values('decile'))
        for c in ('NumPositions', 'NumPool', 'NumObserved'):
            global_decile[c] = global_decile[c].astype(np.int64)
        global_decile['Replace'] = global_decile['NumPool'] < global_decile['NumObserved']
        pool_per_gene = df_pool.groupby(['gene', 'decile']).size().reset_index(name='NumPool').pivot(
            index='gene', columns='decile', values='NumPool').fillna(0).reset_index()

        # Parquet needs string column names: the per-gene tables have one column per decile
        for tbl in (per_gene_pivot, pool_per_gene):
            tbl.columns = [str(c) for c in tbl.columns]
        global_decile.to_parquet(f"{args.output_prefix}_cons_decile_global.parquet",   index=False)
        per_gene_pivot.to_parquet(f"{args.output_prefix}_cons_decile_per_gene.parquet", index=False)
        pool_per_gene.to_parquet(f"{args.output_prefix}_cons_decile_pool_per_gene.parquet", index=False)

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
    parser.add_argument('--null-cycle-caas', dest='null_cycle_caas', default=None,
                        help='perm_pos_cycle_caas.tsv.gz of a CAAS_CORE_MERGE run (Gene, Position, side, '
                             'cycle, caas_score): the filtered null that SCORING also uses. Required for '
                             '--randomization-type permulation; for naive and cons_decile it restricts the '
                             'draw pool to the universe U (null-positive tested positions plus observed '
                             'positions).')
    parser.add_argument('--cycle-roster', dest='cycle_roster', default=None,
                        help='cycle_roster.txt of the same run (one replayed cycle per line): N of the '
                             'permulation null, which counts the cycles that detected nothing. Without '
                             'it N is the number of cycles present in --null-cycle-caas, with a warning.')
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
                             'The discovery table is per-side (a "both" position is two '
                             'rows, side=top and side=bottom), so "top"/"bottom" match '
                             '`side` exactly. "both" (default) retains all non-none rows.')
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')
    main(args)
