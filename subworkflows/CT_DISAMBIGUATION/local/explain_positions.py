#!/usr/bin/env python3
# explain_positions.py — Evidence of the best-scored positions: what each domain of each hypothesis saw.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/

"""
CAAS_EVIDENCE: Evidence of the N best positions of a run: what each domain of each hypothesis saw.

Picks the N best positions of position_scores.tsv (`core.evidence.select_top_positions`), takes their rows from
discovery.tab, scores them again through `core.observed` with the unpooled rows kept and writes
`evidence_top<N>.tsv` (one row per entry and domain, `core.evidence.EVIDENCE_COLUMNS`) and `top_positions.tsv`
(gene, position, CAAS_score, p.emp in rank order). The scoring is the one of observed_b0_main.py, so the numbers
are those of the master before the hypotheses of a position are pooled. Only the schemes a position was scored with
are explained (`scheme_set` of position_scores.tsv, union over its sides; discovery.tab can hold more). The PSS
weights act only in the pooling, so they are not an input. The tables are written with whatever could be explained;
a chosen gene with no alignment, ASR or discovery rows, or a chosen position with no discovery row, makes the run
exit 1 after writing them, so that missing evidence never goes unnoticed.

Called by:  CAAS_EVIDENCE (subworkflows/CT_DISAMBIGUATION/ct_evidence.nf) → explain_positions.py
Inputs:     --discovery        discovery.tab of the run
            --position-scores  position_scores.tsv of the run
            --design           observed design (trait file, or directory of traitfile_H*.tab)
            --alignment-dir, --tree, --asr-cache-dir   alignments, species tree and ASR cache
            --top              number of positions to explain (0 writes empty tables)
Outputs:    <output-dir>/evidence_top<N>.tsv and <output-dir>/top_positions.tsv
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import csv
import logging
import multiprocessing as mp
import sys
import time
from pathlib import Path

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "src"))

from observed_b0_main import read_b0_rows
from src.core.driver import load_gene_context
from src.core.evidence import EVIDENCE_COLUMNS, evidence_rows, select_top_positions
from src.core.labelings import read_trait_pairs
from src.core.observed import analyze_observed, observed_entries
from src.data.loaders import load_ensembl_genes
from src.utils.concurrency import init_worker, plan_concurrency
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)

TOP_COLUMNS = ["gene", "position", "CAAS_score", "p.emp"]  # columns of top_positions.tsv


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--alignment-dir", required=True, help="Directory with alignment files")
    p.add_argument("--tree", required=True, help="Phylogenetic tree (Newick)")
    p.add_argument("--discovery", required=True, help="discovery.tab of the run")
    p.add_argument("--position-scores", required=True, help="position_scores.tsv of the run")
    p.add_argument("--top", type=int, required=True, help="Number of positions to explain (0 writes empty tables)")
    p.add_argument("--design", required=True, help="Observed design: the trait file, or the directory of traitfile_H*.tab")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--asr-model", default="lg")
    p.add_argument("--asr-cache-dir", required=True, help="ASR cache dir")
    p.add_argument("--posterior-threshold", type=float, default=0.0)
    p.add_argument("--taxid-mapping", default=None)
    p.add_argument("--ensembl-genes-file", default=None)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--log-file", type=Path, default=None)
    return p.parse_args()


# ── Evidence ──────────────────────────────────────────────────────────────────


def _explain_gene(job):
    """Pool worker: unpooled evidence rows of one gene (job is the tuple built in main), or (gene, None) when it cannot be scored."""
    (gene, rows, alignment_dir, tree, taxid, model, cache, threshold, ensembl, trait_pairs) = job
    try:
        ctx = load_gene_context(gene, alignment_dir, tree, taxid, model, cache, threshold, ensembl)
        if ctx is None:
            return gene, None
        _, diag = analyze_observed(ctx, gene, observed_entries(gene, rows), trait_pairs, None, threshold, keep_unpooled=True)
        return gene, evidence_rows(diag["unpooled"])
    except Exception as exc:  # noqa: BLE001 - one gene must not end the run; the summary reports it
        logger.error(f"[evidence] {gene} failed: {exc}", exc_info=True)
        return gene, None


def scored_schemes(score_rows):
    """{(gene, position): schemes} a position was scored with: the union over its sides of `scheme_set` ('GS1+US').
    A position whose scheme_set is missing or NA is not in the mapping (no restriction)."""
    out = {}
    for r in score_rows:
        names = {s for s in str(r.get("scheme_set") or "").split("+") if s and s != "NA"}
        if names:
            out.setdefault((str(r["Gene"]), str(r["Position"])), set()).update(names)
    return out


def _write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        w.writerows([[r[c] for c in columns] for r in rows])


def main():
    """Select the top positions, explain them gene by gene and write the evidence and top-position tables."""
    args = parse_arguments()
    configure_logging(verbose=args.verbose, log_file=args.log_file)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.position_scores, newline="") as fh:
        score_rows = list(csv.DictReader(fh, delimiter="\t"))
    chosen = select_top_positions(score_rows, args.top)
    _write_tsv(out_dir / "top_positions.tsv", TOP_COLUMNS,
               [{"gene": g, "position": p, "CAAS_score": repr(i["score"]),
                 "p.emp": "NA" if i["p_emp"] is None else repr(i["p_emp"])} for g, p, i in chosen])
    evidence_path = out_dir / f"evidence_top{args.top}.tsv"
    if not chosen:
        _write_tsv(evidence_path, EVIDENCE_COLUMNS, [])
        return

    wanted = {}
    for gene, pos, _ in chosen:
        wanted.setdefault(gene, set()).add(pos)
    by_gene = read_b0_rows(args.discovery)
    rows_of = {g: [r for r in by_gene.get(g, []) if str(r["position"]) in pos] for g, pos in wanted.items()}
    absent_positions = []
    for gene, pos in sorted(wanted.items()):
        absent = sorted(pos - {str(r["position"]) for r in rows_of[gene]}, key=int)
        if absent:
            absent_positions.append(f"{gene}:{','.join(absent)}")
            logger.error(f"[evidence] {gene}: no discovery rows at position(s) {','.join(absent)}")

    trait_pairs = read_trait_pairs(Path(args.design))
    ensembl = (load_ensembl_genes(Path(args.ensembl_genes_file)) or set()) if args.ensembl_genes_file else None
    jobs = [(g, rows_of[g], args.alignment_dir, args.tree, args.taxid_mapping, args.asr_model, args.asr_cache_dir,
             args.posterior_threshold, ensembl, trait_pairs) for g in sorted(rows_of) if rows_of[g]]
    workers, _ = plan_concurrency(args.workers, 1, logger)
    t0 = time.time()
    by_gene_evidence, skipped = {}, []
    with mp.Pool(processes=workers, initializer=init_worker, initargs=(1, None)) as pool:
        for gene, rows in pool.imap_unordered(_explain_gene, jobs, chunksize=1):
            if rows is None:
                skipped.append(gene)
            else:
                by_gene_evidence[gene] = rows
    skipped += sorted(g for g in wanted if not rows_of[g])
    logger.info(f"[evidence] {len(by_gene_evidence)} genes explained in {time.time() - t0:.1f}s; {len(skipped)} left out"
                + (f" (no alignment, ASR or discovery rows): {sorted(skipped)[:10]}" if skipped else ""))

    # rank order, not arrival order: the rows of a position keep the order evidence_rows gave them
    rank = {(g, p): k for k, (g, p, _) in enumerate(chosen)}
    schemes = scored_schemes(score_rows)
    flat = [r for rows in by_gene_evidence.values() for r in rows
            if (r["gene"], r["msa_pos"]) not in schemes or r["caap_group"] in schemes[(r["gene"], r["msa_pos"])]]
    flat.sort(key=lambda r: rank[(r["gene"], r["msa_pos"])])
    _write_tsv(evidence_path, EVIDENCE_COLUMNS, flat)
    if skipped or absent_positions:
        logger.error(f"[evidence] incomplete: genes left out {sorted(skipped)}; positions without discovery rows {absent_positions}")
        sys.exit(1)


if __name__ == "__main__":
    main()
