#!/usr/bin/env python3
"""Blockwise counters for direct position p: nd and inclusive k_emp.

U is the null-positive positional universe plus observed positions, as before.
The legacy null_classes.parquet/null_fit.tsv staging names and CLI arguments are
retained for Nextflow compatibility; a single inert class is written. No class
pool, conditional tail or factored p is computed. Historical c_k columns are ignored.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

# core.scores lives with the disambiguation core; staged in a work dir the script finds the copy the process places
# next to it, run from the repository it falls back to CT_DISAMBIGUATION/local.
_here = Path(__file__).resolve()
for _cand in (_here.parent, _here.parents[3] / "CT_DISAMBIGUATION" / "local"):
    if (_cand / "src" / "core" / "scores.py").exists():
        sys.path.insert(0, str(_cand))
        break
from src.core.scores import SCORE_RULE, TIE_TOL, cycle_order_key, position_score  # noqa: E402

RATE_MIN = 0.5                  # share of null-positive positions that must lie in the tested set to restrict the null to it
POS_BITS = 32                   # position keys: gene id << 32 | position


# ── Observed positions ────────────────────────────────────────────────────────


def observed_best(postproc_file: Path) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """The observed positive positions: (genes, positions, best-side score), one row per (gene, position).

    The score of a position side is its asr_path_score when positive (core.scores.position_score); only the sides top and
    bottom count, and the position keeps its largest side, the pooled statistic of scoring_compute.R.
    """
    import csv

    best: Dict[Tuple[str, int], float] = {}
    seen = set()
    with open(postproc_file, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["Gene"], int(row["Position"]), row["side"])
            if key in seen:
                raise ValueError(f"duplicate row for {key}: a position side is scored once")
            seen.add(key)
            if row["side"] not in ("top", "bottom"):
                continue
            try:
                score = position_score(float(row["asr_path_score"]))
            except (TypeError, ValueError):
                score = None
            if score is not None:
                k = (key[0], key[1])
                best[k] = max(best.get(k, score), score)
    keys = sorted(best)
    return [g for g, _ in keys], np.array([p for _, p in keys], dtype=np.int64), np.array([best[k] for k in keys], dtype=np.float64)


class GeneCodes:
    """A gene vocabulary: position keys are ``code << 32 | position`` int64."""

    def __init__(self, genes):
        self.genes = sorted(set(genes))
        self._arr = pa.array(self.genes, pa.string())

    def codes(self, column) -> np.ndarray:
        col = column.combine_chunks() if isinstance(column, pa.ChunkedArray) else column
        idx = pc.index_in(col, value_set=self._arr).to_numpy(zero_copy_only=False)
        if np.isnan(idx.astype(float)).any():
            raise ValueError("a gene is missing from the vocabulary")
        return idx.astype(np.int64)

    def keys(self, genes, positions) -> np.ndarray:
        g = self.codes(pa.array(list(genes), pa.string()) if not isinstance(genes, (pa.Array, pa.ChunkedArray)) else genes)
        return (g << POS_BITS) + np.asarray(positions, dtype=np.int64)


def tested_membership(background_file: Path, gene_idx: np.ndarray, gene_names: List[str], positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """For each (gene, position) row: whether CAAStools tested it, and whether its gene has a line in the background file.

    The rows name their gene by an index into `gene_names`. background.output is one line per gene,
    ``gene<TAB>comma-separated positions``; only the genes of `gene_names` are parsed.
    """
    wanted = set(gene_names)
    tested: Dict[str, np.ndarray] = {}
    with open(background_file) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            gene, tab, rest = line.partition("\t")
            if not tab or gene not in wanted:
                continue
            tested[gene] = np.array([int(p) for p in rest.split(",") if p.strip().lstrip("-").isdigit()], dtype=np.int64)
    in_tested = np.zeros(len(positions), dtype=bool)
    covered = np.zeros(len(positions), dtype=bool)
    order = np.argsort(gene_idx, kind="stable")
    bounds = np.searchsorted(gene_idx[order], np.arange(len(gene_names) + 1), side="left")
    for g, name in enumerate(gene_names):
        if name in tested and bounds[g + 1] > bounds[g]:
            rows = order[bounds[g]:bounds[g + 1]]
            covered[rows] = True
            in_tested[rows] = np.isin(positions[rows], tested[name])
    return in_tested, covered


def build_universe(null: pa.Table, obs_genes: List[str], obs_pos: np.ndarray, tested_file: Optional[Path],
                   universe_file: Optional[Path], roster: List[str]):
    """The positions of U with their class, the fit and the base metrics of the null.

    Returns (classes table (gene, position, nd, cls), breaks, class sizes {cls: detections}, metrics rows, note).
    """
    encoded = null["gene"].combine_chunks().dictionary_encode()
    gene_names = encoded.dictionary.to_pylist()
    gene_idx = encoded.indices.to_numpy(zero_copy_only=False).astype(np.int64)
    positions = null["position"].to_numpy(zero_copy_only=False).astype(np.int64)
    nd = null["nd"].to_numpy(zero_copy_only=False).astype(np.int64)
    first = null["first_rank"].to_numpy(zero_copy_only=False).astype(np.int64)
    n_raw = len(positions)
    vocab = GeneCodes(gene_names + list(obs_genes))
    name_code = vocab.codes(pa.array(gene_names, pa.string())) if gene_names else np.empty(0, dtype=np.int64)
    null_keys = (name_code[gene_idx] << POS_BITS) + positions if n_raw else np.empty(0, dtype=np.int64)
    obs_keys = np.unique(vocab.keys(obs_genes, obs_pos)) if len(obs_genes) else np.empty(0, dtype=np.int64)

    note = "no tested positions given: every null-positive position"
    keep = np.ones(n_raw, dtype=bool)
    if tested_file is not None and Path(tested_file).is_file():
        in_tested, covered = tested_membership(tested_file, gene_idx, gene_names, positions)
        n_cov = int(covered.sum())
        rate = in_tested[covered].sum() / n_cov if n_cov > 0 else 0.0
        if rate >= RATE_MIN:
            if universe_file is not None and Path(universe_file).is_file():
                u = {g.strip() for g in Path(universe_file).read_text().splitlines() if g.strip() and g.strip() != "Gene"}
                in_univ = np.array([g in u for g in gene_names], dtype=bool)[gene_idx] if n_raw else np.ones(0, dtype=bool)
                note = "null-positive positions inside the tested set and the cleaned background genes"
            else:
                in_univ = np.ones(n_raw, dtype=bool)
                note = "null-positive positions inside the tested set"
            keep = (in_tested & in_univ) | np.isin(null_keys, obs_keys)
        else:
            note = (f"tested positions ignored (only {100 * rate:.1f}% of the null-positive positions lie in them): "
                    "every null-positive position")

    order = sorted(roster, key=cycle_order_key)
    n_cycles = len(order)
    # seq(N / 10, N, length.out = 10) as R builds it: from + i * by, never above N
    by = (n_cycles - n_cycles / 10) / 9
    ks = sorted({max(1, int(round(min(n_cycles / 10 + i * by, n_cycles)))) for i in range(10)}) if n_cycles else []
    first_u = first[keep] + 1
    curve = [(k, int((first_u <= k).sum())) for k in ks]

    nd_u = nd[keep]
    breaks = np.empty(0)  # inert legacy fit; not a statistical model
    cls = np.ones(len(nd_u), dtype=np.int64)
    sizes = {int(k): int(nd_u[cls == k].sum()) for k in np.unique(cls)}
    classes = pa.table({"gene": null["gene"].filter(pa.array(keep)), "position": pa.array(positions[keep].astype(np.int32)),
                        "nd": pa.array(nd_u.astype(np.int32)), "cls": pa.array(cls.astype(np.int32))})
    metrics = {"n_null_raw": n_raw, "n_null_u": int(keep.sum()), "curve": curve, "rate_note": note}
    return classes, breaks, sizes, metrics


def gene_list_identity(path: Path) -> Tuple[int, str]:
    """(number of genes, digest) of a gene list file: one gene per line, optional 'Gene' header.

    The digest is the SHA-256 of the sorted unique names joined by newlines, the one posenrich_universe.gene_list_digest
    computes, so consumers compare lists by content and not by file layout.
    """
    import hashlib

    genes = {g.strip() for g in Path(path).read_text().splitlines() if g.strip() and g.strip() != "Gene"}
    return len(genes), hashlib.sha256("\n".join(sorted(genes)).encode()).hexdigest()


def write_universe(out: Path, classes, breaks, sizes, metrics, n_cycles: int, n_present: int, gene_list=None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    pq.write_table(classes, out / "null_classes.parquet", compression="zstd")
    with open(out / "null_fit.tsv", "w") as fh:
        fh.write("cls\tupper\tsize\n")
        n_cls = len(breaks) + 1
        for k in range(1, n_cls + 1):
            upper = repr(float(breaks[k - 1])) if k <= len(breaks) else "Inf"
            fh.write(f"{k}\t{upper}\t{sizes.get(k, 0)}\n")
    with open(out / "null_universe_base.tsv", "w") as fh:
        fh.write("metric\tcycles\tvalue\tnote\n")
        fh.write(f"N_cycles\tNA\t{n_cycles}\tcycle_roster.txt\n")
        fh.write(f"N_cycles_present_in_null\tNA\t{n_present}\t\n")
        fh.write(f"null_positive_positions_all\tNA\t{metrics['n_null_raw']}\t\n")
        fh.write(f"null_positive_positions_in_universe\tNA\t{metrics['n_null_u']}\t{metrics['rate_note']}\n")
        if gene_list is not None:  # value = number of genes, note = digest of the gene list the universe received
            fh.write(f"gene_list_sha256\tNA\t{gene_list[0]}\t{gene_list[1]}\n")
        for k, v in metrics["curve"]:
            fh.write(f"null_positive_positions_cumulative\t{k}\t{v}\tdistinct null-positive positions of U seen in the first k cycles\n")


# ── Counting one block ────────────────────────────────────────────────────────


def count_block(positions_file: Path, classes: pa.Table, breaks: np.ndarray, obs_genes: List[str], obs_pos: np.ndarray,
                obs_score: np.ndarray) -> pa.Table:
    """Direct exceedance counts from one block; legacy class arguments are ignored."""
    table = pq.read_table(positions_file, columns=["gene", "position", "all"])
    vocab = GeneCodes(set(pc.unique(table["gene"]).to_pylist()) | set(pc.unique(classes["gene"]).to_pylist()) | set(obs_genes))
    keys = vocab.keys(table["gene"], table["position"].to_numpy(zero_copy_only=False))
    value = table["all"].to_numpy(zero_copy_only=False).astype(np.float64)
    obs_keys = vocab.keys(obs_genes, obs_pos) if len(obs_genes) else np.empty(0, dtype=np.int64)
    thr = obs_score - TIE_TOL

    # k_emp: the cycles of the block whose best-side score of the position reaches the observed one (the rows of an
    # observed position are always kept: the universe never removes them)
    k_emp = np.zeros(len(obs_keys), dtype=np.int64)
    if len(keys) and len(obs_keys):
        order = np.argsort(keys, kind="stable")
        sk, sv = keys[order], value[order]
        lo = np.searchsorted(sk, obs_keys, side="left")
        hi = np.searchsorted(sk, obs_keys, side="right")
        for i in np.flatnonzero(hi > lo).tolist():
            k_emp[i] = int((sv[lo[i]:hi[i]] >= thr[i]).sum())

    cols = {"gene": pa.array(obs_genes, pa.string()), "position": pa.array(obs_pos.astype(np.int32)),
            "k_emp": pa.array(k_emp)}
    return pa.table(cols)


# ── Reduce ────────────────────────────────────────────────────────────────────


def reduce_counts(count_files: List[Path], classes: pa.Table, obs_genes: List[str], obs_pos: np.ndarray,
                  obs_score: np.ndarray, n_cls: int) -> Dict[str, np.ndarray]:
    """Sum the per-block counts and attach the nd of each observed position: the row of null_counters.tsv."""
    n = len(obs_genes)
    k_emp = np.zeros(n, dtype=np.int64)
    for f in count_files:
        t = pq.read_table(f)
        if t.num_rows != n:
            raise ValueError(f"{f} counts {t.num_rows} positions, {n} observed")
        k_emp += t["k_emp"].to_numpy(zero_copy_only=False)
    vocab = GeneCodes(set(pc.unique(classes["gene"]).to_pylist()) | set(obs_genes))
    nd = np.zeros(n, dtype=np.int64)
    if classes.num_rows and n:
        ckeys = vocab.keys(classes["gene"], classes["position"].to_numpy(zero_copy_only=False))
        corder = np.argsort(ckeys)
        ckeys = ckeys[corder]
        cnd = classes["nd"].to_numpy(zero_copy_only=False).astype(np.int64)[corder]
        okeys = vocab.keys(obs_genes, obs_pos)
        at = np.searchsorted(ckeys, okeys)
        at[at >= len(ckeys)] = 0
        hit = ckeys[at] == okeys
        nd = np.where(hit, cnd[at], 0)
    return {"k_emp": k_emp, "nd": nd}


def write_counters(path: Path, obs_genes, obs_pos, red, n_cls: int) -> None:
    with open(path, "w") as fh:
        fh.write("Gene\tPosition\tk_emp\tnd\n")
        for i, (g, p) in enumerate(zip(obs_genes, obs_pos.tolist())):
            fh.write(f"{g}\t{p}\t{int(red['k_emp'][i])}\t{int(red['nd'][i])}\n")



def read_fit(path: Path) -> Tuple[np.ndarray, Dict[int, int]]:
    breaks, sizes = [], {}
    for line in path.read_text().splitlines()[1:]:
        k, upper, size = line.split("\t")
        if upper != "Inf":
            breaks.append(float(upper))
        sizes[int(k)] = int(size)
    return np.asarray(breaks, dtype=np.float64), sizes


# ── CLI ───────────────────────────────────────────────────────────────────────


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    u = sub.add_parser("universe", help="nd -> U -> classes")
    u.add_argument("--null-positions", required=True, type=Path)
    u.add_argument("--postproc", required=True, type=Path, help="filtered_discovery.tsv (the observed positions)")
    u.add_argument("--tested", default=None, type=Path, help="caastools background.output")
    u.add_argument("--universe", default=None, type=Path, help="cleaned_background_main.txt")
    u.add_argument("--roster", required=True, type=Path, help="cycle_roster.txt")
    u.add_argument("--counts-json", default=None, type=Path, help="counts.json of reduce_null_blocks.py (cycles with a row)")
    u.add_argument("--output-dir", required=True, type=Path)
    c = sub.add_parser("count", help="direct k_emp of one block")
    c.add_argument("--positions", required=True, type=Path)
    c.add_argument("--classes", required=True, type=Path)
    c.add_argument("--fit", required=True, type=Path)
    c.add_argument("--postproc", required=True, type=Path)
    c.add_argument("--output", required=True, type=Path)
    r = sub.add_parser("reduce", help="sum the counts of the blocks")
    r.add_argument("--counts", required=True, nargs="+", type=Path)
    r.add_argument("--classes", required=True, type=Path)
    r.add_argument("--fit", required=True, type=Path)
    r.add_argument("--postproc", required=True, type=Path)
    r.add_argument("--output", required=True, type=Path)
    a = ap.parse_args(argv)

    obs_genes, obs_pos, obs_score = observed_best(a.postproc)
    if a.cmd == "universe":
        roster = [x for x in a.roster.read_text().splitlines() if x.strip()]
        if len(roster) != len(set(roster)):
            raise ValueError("Duplicate null cycle IDs in roster")
        null = pq.read_table(a.null_positions)
        classes, breaks, sizes, metrics = build_universe(null, obs_genes, obs_pos, a.tested, a.universe, roster)
        present = len(roster)
        if a.counts_json is not None and a.counts_json.is_file():
            counts = json.loads(a.counts_json.read_text())
            if counts["cycles"] != len(roster):
                raise ValueError("Null cycle count differs from roster")
            present = int(counts["cycles_present"])
            if counts.get("score_rule") != SCORE_RULE:
                print(f"the CAAS permulation null was scored with the rule '{counts.get('score_rule')}' and the observed scores with "
                      f"'{SCORE_RULE}': p.emp and p.adj_bh would compare different statistics. Rebuild the null (CAAS_CORE_MERGE).",
                      file=sys.stderr)
                return 1
        gene_list = gene_list_identity(a.universe) if a.universe is not None and Path(a.universe).is_file() else None
        write_universe(a.output_dir, classes, breaks, sizes, metrics, len(roster), present, gene_list)
        print(f"[null] universe: {metrics['n_null_u']} of {metrics['n_null_raw']} null-positive positions kept "
              f"({metrics['rate_note']}); {len(breaks) + 1} classes", file=sys.stderr)
        return 0
    breaks, sizes = read_fit(a.fit)
    classes = pq.read_table(a.classes)
    if a.cmd == "count":
        pq.write_table(count_block(a.positions, classes, breaks, obs_genes, obs_pos, obs_score), a.output, compression="zstd")
        return 0
    red = reduce_counts(a.counts, classes, obs_genes, obs_pos, obs_score, len(breaks) + 1)
    write_counters(a.output, obs_genes, obs_pos, red, len(breaks) + 1)
    print(f"[null] counters of {len(obs_genes)} observed positions from {len(a.counts)} blocks", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
