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
import os
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
pa.set_cpu_count(max(1,int(os.environ.get('OMP_NUM_THREADS','1'))))
pa.set_io_thread_count(max(1,int(os.environ.get('OMP_NUM_THREADS','1'))))

# core.scores lives with the disambiguation core; staged in a work dir the script finds the copy the process places
# next to it, run from the repository it falls back to CT_DISAMBIGUATION/local.
_here = Path(__file__).resolve()
for _cand in (_here.parent, _here.parents[3] / "CT_DISAMBIGUATION" / "local"):
    if (_cand / "src" / "core" / "scores.py").exists():
        sys.path.insert(0, str(_cand))
        break
from src.core.scores import SCORE_RULE, TIE_TOL, cycle_order_key, position_score  # noqa: E402

sys.path.insert(0, str(_here.parent))
from scoring_input import read_rows
POS_BITS = 32                   # position keys: gene id << 32 | position


# ── Observed positions ────────────────────────────────────────────────────────


def observed_best(postproc_file: Path) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """The observed positive positions: (genes, positions, best-side score), one row per (gene, position).

    The score of a position side is its asr_path_score when positive (core.scores.position_score); only the sides top and
    bottom count, and the position keeps its largest side, the pooled statistic of scoring_compute.R.
    """
    best: Dict[Tuple[str, int], float] = {}
    for (g,p,side), value in read_rows(postproc_file).items():
        if side in ('top','bottom'):
            best[(g,p)] = max(best.get((g,p),value),value)
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


def read_tested(path):
    tested = {}
    for line in Path(path).read_text().splitlines():
        g, tab, text = line.partition('\t')
        if not tab or not g or g.strip()!=g or g in tested:
            raise ValueError('Invalid/duplicate tested gene')
        values = text.split(',') if text not in ('','NULL') else []
        if any(not p.isdigit() or int(p)>=2**31 for p in values):
            raise ValueError('Invalid zero-based tested position')
        ps = [int(p) for p in values]
        if len(set(ps))!=len(ps):
            raise ValueError('Duplicate tested position')
        tested[g] = np.asarray(ps,dtype=np.int64)
    return tested


def read_universe(path):
    genes = Path(path).read_text().splitlines()
    if genes and genes[0] in ('Gene','gene'):
        genes = genes[1:]
    if any(not g or g.strip()!=g for g in genes) or len(set(genes))!=len(genes):
        raise ValueError('Invalid/duplicate universe gene')
    return set(genes)


def tested_membership(background_file, gene_idx, gene_names, positions):
    tested = read_tested(background_file)
    inside = np.zeros(len(positions),dtype=bool)
    covered = np.zeros(len(positions),dtype=bool)
    order = np.argsort(gene_idx,kind='stable')
    bounds = np.searchsorted(gene_idx[order],np.arange(len(gene_names)+1))
    for g,name in enumerate(gene_names):
        rows = order[bounds[g]:bounds[g+1]]
        if name in tested:
            covered[rows] = True
            inside[rows] = np.isin(positions[rows],tested[name])
    return inside, covered


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

    if tested_file is None or universe_file is None:
        raise ValueError('Calibrated positions require explicit tested positions and cleaned background; no detection-selected BH family')
    tested = read_tested(tested_file)
    universe = read_universe(universe_file)
    for g,p in zip(obs_genes,obs_pos):
        if g not in universe or g not in tested or p not in tested[g]:
            raise ValueError(f'Observed position outside declared tested universe: {(g,int(p))}')
    in_tested, _ = tested_membership(tested_file,gene_idx,gene_names,positions)
    in_univ = np.asarray([g in universe for g in gene_names],dtype=bool)[gene_idx]
    # Null-only genes/positions outside the declared family are excluded explicitly.
    keep = in_tested & in_univ
    family = sum(len(tested.get(g,())) for g in universe)
    note = 'fixed tested positions intersected with cleaned background genes'

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
    metrics = {"n_null_raw": n_raw, "n_null_u": int(keep.sum()), "curve": curve, "rate_note": note, "family": family}
    return classes, breaks, sizes, metrics


def gene_list_identity(path: Path) -> Tuple[int, str]:
    """(number of genes, digest) of a gene list file: one gene per line, optional 'Gene' header.

    The digest is the SHA-256 of the sorted unique names joined by newlines, the one posenrich_universe.gene_list_digest
    computes, so consumers compare lists by content and not by file layout.
    """
    import hashlib

    genes = read_universe(path)
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
        fh.write(f"tested_positions_family\tNA\t{metrics['family']}\tfixed tested family; undetected observed positions have p=1\n")
        fh.write(f"null_positive_positions_all\tNA\t{metrics['n_null_raw']}\t\n")
        fh.write(f"null_positive_positions_in_universe\tNA\t{metrics['n_null_u']}\t{metrics['rate_note']}\n")
        if gene_list is not None:  # value = number of genes, note = digest of the gene list the universe received
            fh.write(f"gene_list_sha256\tNA\t{gene_list[0]}\t{gene_list[1]}\n")
        for k, v in metrics["curve"]:
            fh.write(f"null_positive_positions_cumulative\t{k}\t{v}\tdistinct null-positive positions of U seen in the first k cycles\n")


def file_digest(path):
    import hashlib
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for chunk in iter(lambda:fh.read(1048576),b''):
            h.update(chunk)
    return h.hexdigest()


def validated_positions(path):
    table = pq.read_table(path, columns=['gene','cycle','position','all'])
    if any(table[n].null_count for n in table.column_names) or not pa.types.is_integer(table['position'].type):
        raise ValueError('Invalid null position schema or missing values')
    for name in ('gene','cycle'):
        if not (pa.types.is_string(table[name].type) or pa.types.is_dictionary(table[name].type)):
            raise ValueError('Invalid null position IDs')
        values = pc.unique(table[name]).to_pylist()
        if any(not isinstance(v,str) or not v or v.strip()!=v for v in values):
            raise ValueError('Invalid null position IDs')
    p = table['position'].to_numpy(zero_copy_only=False)
    v = table['all'].to_numpy(zero_copy_only=False)
    if np.any(p<0) or np.any(p>=2**31) or np.any(~np.isfinite(v)) or np.any(v<=0) or np.any(v>1):
        raise ValueError('Invalid null position/score')
    unique = table.group_by(['gene','cycle','position']).aggregate([('all','count')])
    if unique.num_rows!=table.num_rows:
        raise ValueError('Duplicate gene/cycle/position in null')
    return table


def universe_digest(gene, positions):
    """Identity of a gene's positional universe: SHA-256 of its canonical background.output line (ascending positions,
    NULL when empty). The null tiles record it with src/core/universe.py:digest, which defines the same string."""
    import hashlib
    ordered = sorted(int(p) for p in positions)
    return hashlib.sha256(f"{gene}\t{','.join(map(str, ordered)) if ordered else 'NULL'}".encode()).hexdigest()


def validate_null(blocks, roster, null, counts, output=None, observed_contract=None, tested=None):
    import re
    if not roster or len(set(roster))!=len(roster) or any(not re.fullmatch(r'b_[1-9][0-9]*',c) for c in roster):
        raise ValueError('Invalid/empty null cycle roster')
    rank = {c:i for i,c in enumerate(sorted(roster,key=cycle_order_key))}
    acc, seen, present, total = {}, set(), set(), 0
    hashes = []
    for i,block in enumerate(blocks,1):
        cy = (block/'cycles.txt').read_text().splitlines()
        if not cy or len(set(cy))!=len(cy) or seen.intersection(cy) or not set(cy)<=set(roster):
            raise ValueError('Missing/mixed/overlapping block cycles')
        seen.update(cy)
        meta = json.loads((block/'counts.json').read_text())
        if meta.get('score_rule')!=SCORE_RULE or meta.get('cycles')!=len(cy):
            raise ValueError('Block count/rule mismatch')
        if observed_contract is not None:
            manifest = block/'scoring_artifacts.json'
            if not manifest.is_file():
                raise ValueError('Missing null artifact checksums; regenerate canonical blocks')
            saved = json.loads(manifest.read_text())
            required = {'matrices.rds','positions.parquet','cycles.txt','counts.json'}
            if saved.get('version')!=1 or not isinstance(saved.get('files'),dict) or not required<=set(saved['files']):
                raise ValueError('Invalid null artifact manifest')
            import hashlib
            for name,checksum in saved['files'].items():
                if Path(name).name!=name:
                    raise ValueError('Invalid null artifact path')
                h=hashlib.md5()
                with (block/name).open('rb') as fh:
                    for chunk in iter(lambda:fh.read(1048576),b''):h.update(chunk)
                if h.hexdigest()!=checksum:
                    raise ValueError(f'Null artifact checksum mismatch: {name}')
            prov = meta.get('scoring_provenance')
            if not isinstance(prov,dict) or prov.get('version')!=1:
                raise ValueError('Unverifiable legacy null scientific/filter provenance; regenerate observed and null artifacts')
            for key in ('scientific','cluster','gene_filter','strict'):
                if prov.get(key)!=observed_contract.get(key) or prov.get(key) is None:
                    raise ValueError(f'Observed/null {key} mismatch')
            for key in ('maps','lengths','length_sources','alignments'):
                for gene in set(prov.get(key,{})) & set(observed_contract.get(key,{})):
                    if prov[key][gene]!=observed_contract[key][gene]:
                        raise ValueError(f'Observed/null {key} mismatch for {gene}')
            # The null was scored over the positions of the positional universe: it is valid only with that universe.
            used = prov.get('universe')
            if not isinstance(used, dict) or not used:
                raise ValueError('Null without positional-universe provenance; regenerate the null with the positional universe')
            if tested is None:
                raise ValueError('The tested positions (background.output) are required to check the null positional universe')
            for gene, recorded in used.items():
                if gene not in tested or universe_digest(gene, tested[gene]) != recorded:
                    raise ValueError(f'Null positional universe differs from the tested positions (background.output) for {gene}')
        t = validated_positions(block/'positions.parquet')
        if meta.get('scored_positions')!=t.num_rows:
            raise ValueError('Block scored-position count mismatch')
        total += t.num_rows
        cp = meta.get('cycles_present')
        if not isinstance(cp,list) or len(cp)!=len(set(cp)) or not set(cp)<=set(cy):
            raise ValueError('Invalid block cycles_present')
        present.update(cp)
        if not set(pc.unique(t['cycle']).to_pylist())<=set(cp):
            raise ValueError('Position cycle outside block roster/present cycles')
        maxima = t.group_by(['gene','cycle']).aggregate([('all','max')]).sort_by([('gene','ascending'),('cycle','ascending')])
        if output is not None:
            Path(output).mkdir(parents=True,exist_ok=True)
            with (Path(output)/f'block_{i:06d}.tsv').open('w') as fh:
                fh.write('gene\tcycle\tmax\n')
                for g,c,v in zip(*(maxima[n].to_pylist() for n in ('gene','cycle','all_max'))):
                    fh.write(f'{g}\t{c}\t{v!r}\n')
        ranks = pc.index_in(t['cycle'],value_set=pa.array(sorted(roster,key=cycle_order_key),pa.string()))
        detections = t.append_column('rank',ranks).group_by(['gene','position']).aggregate([('rank','min'),('all','count')])
        for g,p,first,nd in zip(*(detections[n].to_pylist() for n in ('gene','position','rank_min','all_count'))):
            old,earliest = acc.get((g,p),(0,len(roster)))
            acc[g,p] = old+nd,min(earliest,first)
        for name in ('matrices.rds','positions.parquet','cycles.txt','counts.json'):
            hashes.append((block.name,name,file_digest(block/name)))
    if seen!=set(roster) or counts.get('cycles')!=len(roster) or counts.get('cycles_present')!=len(present) or counts.get('scored_positions')!=total:
        raise ValueError('Null aggregate counts/roster disagree with blocks')
    actual = {}
    for g,p,nd,first in zip(*(null[n].to_pylist() for n in ('gene','position','nd','first_rank'))):
        if (g,p) in actual:
            raise ValueError('Duplicate aggregate null positions')
        actual[g,p] = nd,first
    if actual!=acc:
        raise ValueError('Null positions disagree with block content')
    return hashes


# ── Counting one block ────────────────────────────────────────────────────────


def count_block(positions_file: Path, classes: pa.Table, breaks: np.ndarray, obs_genes: List[str], obs_pos: np.ndarray,
                obs_score: np.ndarray) -> pa.Table:
    """Direct exceedance counts from one block; legacy class arguments are ignored."""
    table = validated_positions(positions_file)
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
        if t['gene'].to_pylist()!=list(obs_genes) or t['position'].to_pylist()!=list(obs_pos):
            raise ValueError('Counter keys/order differ from observed positions')
        vals = t['k_emp'].to_numpy(zero_copy_only=False)
        if np.any(~np.isfinite(vals)) or np.any(vals<0) or np.any(vals!=np.floor(vals)):
            raise ValueError('Invalid exceedance counts')
        k_emp += vals.astype(np.int64)
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
    if np.any(k_emp>nd):
        raise ValueError('Exceedance counts exceed detections; duplicate/mixed blocks')
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
    u.add_argument("--counts-json", required=True, type=Path, help="counts.json of reduce_null_blocks.py (cycles with a row)")
    u.add_argument("--observed-contract", required=True,type=Path)
    u.add_argument("--blocks", nargs="+", required=True, type=Path)
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

    if a.cmd == 'universe':
        if not a.counts_json.is_file():
            raise ValueError('counts.json with the CAAS score rule is required; regenerate observed and null artifacts together')
        saved_rule = json.loads(a.counts_json.read_text()).get('score_rule')
        if saved_rule != SCORE_RULE:
            raise ValueError(f'Retired or unavailable CAAS null score rule {saved_rule!r}; expected {SCORE_RULE}. Regenerate observed and null artifacts together')
    obs_genes, obs_pos, obs_score = observed_best(a.postproc)
    if a.cmd == "universe":
        roster = [x for x in a.roster.read_text().splitlines() if x.strip()]
        if len(roster) != len(set(roster)):
            raise ValueError("Duplicate null cycle IDs in roster")
        null = pq.read_table(a.null_positions)
        if a.tested is None or a.universe is None:
            raise ValueError('Calibrated positions require explicit tested positions and cleaned background')
        contract = json.loads(a.observed_contract.read_text())
        if contract.get('version')!=1 or contract.get('observed_sha256')!=file_digest(a.postproc) or not contract.get('scientific'):
            raise ValueError('Observed scoring provenance missing or checksum mismatch; regenerate CT_POSTPROC')
        if set(contract.get('excluded_genes',[])) & read_universe(a.universe):
            raise ValueError('Cleaned background contains explicitly excluded observed genes; provide an evaluable background')
        hashes = validate_null(a.blocks, roster, null, json.loads(a.counts_json.read_text()),a.output_dir,contract,read_tested(a.tested))
        (a.output_dir/'validated_inputs.json').write_text(json.dumps({'version':1,'blocks':hashes,'null_positions_sha256':file_digest(a.null_positions),'counts_sha256':file_digest(a.counts_json),'roster_sha256':file_digest(a.roster),'observed_sha256':file_digest(a.postproc),'tested_sha256':file_digest(a.tested) if a.tested else None,'background_sha256':file_digest(a.universe) if a.universe else None}))
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
