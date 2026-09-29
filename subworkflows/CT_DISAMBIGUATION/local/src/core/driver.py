"""Per-gene replay of labelings: the code the observed labeling (b_0) and the permuted ones share.

One gene's ASR is loaded once (`load_gene_context`), every labeling is scored against it through
`analyze_gene_disambiguation` (`score_labelings`), and the hypotheses of a cycle are collapsed into one
domain-pooled record per position and scheme (`pool_labelings`). Reading discovery files, splitting a
gene into chunks and writing tables stay with the callers.
"""
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from src.asr.asr_single import load_alignment_and_mappings, load_and_match_tree, run_asr_pipeline
from src.convergence.disambiguate_single import PositionAxes, analyze_gene_disambiguation
from src.convergence.fop_pool import base_cycle, pool_domains
from src.core.labelings import hyp_id, trait_pairs_from
from src.core.pooling import pooled_sides
from src.phylo.tree_utils import build_tree_node_mapping, extract_tip_labels
from src.utils.io_utils import find_gene_alignment

logger = logging.getLogger(__name__)


def load_gene_context(
    gene: str,
    alignment_dir: str,
    tree_file: str,
    taxid_mapping_path: Optional[str],
    asr_model: str,
    asr_cache_dir: str,
    posterior_threshold: float,
    ensembl_genes: Optional[Set[str]] = None,
) -> Optional[Dict[str, Any]]:
    """Load alignment + tree + precomputed ASR posteriors for one gene ONCE.

    Mirrors the precomputed-ASR load path of process_single_gene (the
    phenotype-invariant part), including the PAML tree rebuild for node alignment.
    Returns a context dict, or None if alignment/ASR unavailable.
    """
    alignment_path = find_gene_alignment(Path(alignment_dir), gene, ensembl_genes)
    if not alignment_path:
        return None

    alignment_data = load_alignment_and_mappings(
        alignment_path,
        Path(taxid_mapping_path) if taxid_mapping_path else None,
        gene_name=gene,
    )
    tree_data = load_and_match_tree(
        Path(tree_file), alignment_data,
        Path(taxid_mapping_path) if taxid_mapping_path else None,
    )

    from src.asr.asr_single import SingleGeneASRConfig, run_asr_pipeline

    asr_config = SingleGeneASRConfig(
        alignment_path=alignment_path,
        tree_path=Path(tree_file),
        taxid_path=Path(taxid_mapping_path) if taxid_mapping_path else None,
        model=asr_model,
        posterior_threshold=posterior_threshold,
        output_dir=Path(asr_cache_dir),
    )
    # run_asr_pipeline loads the cached ASR when present (rst + rst1) and COMPUTES
    # it into asr_cache_dir on a miss. This makes the permulation replay robust to
    # genes that appear only under a null labeling — never in the observed run, so
    # never cached by CT_DISAMBIGUATION_RUN — which in asr_mode=compute would
    # otherwise be silently dropped from the null. Observed-significant genes are
    # already cached by the time this runs (the asr_ready gate in
    # caas_permulation.nf), so this only computes the null-only tail.
    try:
        node_posteriors = run_asr_pipeline(
            gene, asr_config, skip_if_exists=True,
            alignment_data=alignment_data, tree_data=tree_data,
        )
    except Exception as exc:  # noqa: BLE001 — codeml / parse failure for this one gene
        logger.warning(f"[perms] ASR unavailable for {gene} ({exc}) — excluded from the null")
        return None
    if not node_posteriors:
        return None

    rst_file = getattr(node_posteriors, "rst_file", None)
    paml_tree_file = getattr(node_posteriors, "tree_file", None)
    if rst_file and paml_tree_file and Path(paml_tree_file).exists():
        try:
            ordered_nodes, id_mapping = build_tree_node_mapping(
                tree_file=Path(paml_tree_file), rst_file=Path(rst_file)
            )
            tree_data.nodes = ordered_nodes
            tree_data.root = ordered_nodes[-1]
            tree_data.node_mapping = id_mapping

            def _tip_taxid(label: str) -> str:
                return label.split("_")[-1] if "_" in label else label

            tree_data.tip_set = {
                _tip_taxid(lbl) for lbl in extract_tip_labels(tree_data.root)
            }
        except Exception as e:
            logger.warning(f"[perms] {gene}: could not rebuild tree_data from PAML tree: {e}")

    return {
        "alignment_data": alignment_data,
        "tree_data": tree_data,
        "node_posteriors": node_posteriors,
    }


def score_labelings(
    ctx: Dict[str, Any],
    gene: str,
    tags: List[str],
    labelings: Dict[str, Tuple[List[str], List[str]]],
    entries_by_tag: Dict[str, List[Any]],
    posterior_threshold: float,
) -> List[Tuple[str, List[Any]]]:
    """Score every labeling of `tags` that has discovery entries: [(tag, per-position records)].

    A labeling whose scoring raises is skipped (logged at debug), as the null has always done.
    """
    alignment_data = ctx["alignment_data"]
    tree_data = ctx["tree_data"]
    node_posteriors = ctx["node_posteriors"]
    full_posteriors = getattr(node_posteriors, "posteriors_node", None)

    axes_only = os.environ.get("CAAS_PERMS_AXES_ONLY", "1") not in ("0", "false", "False")
    per_site_dist_cache: Dict[int, Any] = {} if axes_only else None

    results = []
    for cyc in tags:
        labeling = labelings.get(cyc)
        if not labeling:
            continue
        caas_entries = entries_by_tag.get(cyc, [])
        if not caas_entries:
            continue
        fg, bg = labeling
        trait_pairs = trait_pairs_from(fg, bg)
        try:
            biochem_results, _ = analyze_gene_disambiguation(
                gene=gene,
                alignment_data=alignment_data,
                tree_data=tree_data,
                caas_positions=[],
                caas_entries=caas_entries,
                caas_metadata_path=Path("dummy_path"),
                trait_pairs=trait_pairs,
                taxid_mapping=alignment_data.species_to_taxid,
                posterior_data=full_posteriors,
                posterior_threshold=posterior_threshold,
                diagnostics_dir=None,
                asr_mode="precomputed",
                axes_only=axes_only,
                per_site_dist_cache=per_site_dist_cache,
            )
            if biochem_results:
                results.append((cyc, biochem_results))
        except Exception as e:
            logger.debug(f"[perms] {gene} cycle {cyc} failed: {e}")
            continue
    return results


def _expand_pooled(pos, grp, hyp_label, pooled):
    """pool_domains return -> <= 2 per-side PositionAxes rows (`side` authoritative, `asr_path_score` =
    that side's core_s). No participating domain on either side -> one `side="none"` row."""
    out = [
        PositionAxes(
            position=pos, caap_group=grp,
            asr_path_score=sd["asr_path_score"],
            side=sd["side"], hypothesis=hyp_label,
            derived_agreement=sd["derived_agreement"],
            domain_scores=sd["domain_scores"],
        )
        for sd in pooled_sides(pooled)
    ]
    if not out:
        out.append(PositionAxes(
            position=pos, caap_group=grp, asr_path_score=0.0,
            side="none", hypothesis=hyp_label,
        ))
    return out


def pool_labelings(
    results: List[Tuple[str, List[Any]]],
    fop_pairs: Optional[Dict[str, Dict[Tuple[str, int], float]]],
) -> List[Tuple[str, List[Any]]]:
    """Collapse scored labelings into one per-side, domain-pooled record set per cycle.

    Every axes-only record carries `.sides`, the raw compute_domain_scores return
    {"top": {...}, "bottom": {...}, "domain_meta": {...}} for one (position, scheme, hypothesis).
    fop_pool.pool_domains reduces M >= 1 such records to one score per phenotype side (M == 1 is
    the plain PSS-weighted mean over the K fixed Voronoi domains), the same statistic the observed
    path emits through disambiguate_single._emit_pooled_side_rows.

    With `fop_pairs` the "<base>~H<m>" hypotheses of a base cycle are pooled together, weighted by
    their PSS; without it each labeling is its own single-hypothesis cycle.
    """
    if fop_pairs is not None:
        # (base_cyc, pos, scheme) -> [ {hyp, sides} ]; one pool_domains call each.
        by_pos: Dict[Tuple[str, int, str], List[Dict[str, Any]]] = {}
        for cyc, results_list in results:
            base = base_cycle(cyc)
            hyp = hyp_id(cyc.split("~", 1)[1]) if "~" in cyc else "H1"
            for r in results_list:
                pos = getattr(r, "position", None)
                if pos is None:
                    continue
                grp = getattr(r, "caap_group", "US")
                by_pos.setdefault((base, pos, grp), []).append(
                    {"hyp": hyp, "sides": getattr(r, "sides", None) or {}}
                )

        pooled_by_cycle: Dict[str, List[Any]] = {}
        for (base, pos, grp), hyp_recs in by_pos.items():
            # pool_domains wants only {(hyp, domain) -> pss} of this cycle (missing -> equal weight).
            pss_map = fop_pairs.get(base, {}) or None
            pooled = pool_domains(hyp_recs, pss_map)
            pooled_by_cycle.setdefault(base, []).extend(_expand_pooled(pos, grp, None, pooled))
        return list(pooled_by_cycle.items())

    # Non-FOP: one hypothesis per record. pool_domains still runs (M == 1) so the observed and
    # the null go through the exact same reducer.
    expanded = []
    for cyc, recs in results:
        out_recs = []
        for r in recs:
            pos = getattr(r, "position", None)
            if pos is None:
                continue
            hyp = getattr(r, "hypothesis", None)
            pooled = pool_domains(
                [{"hyp": hyp or "H1", "sides": getattr(r, "sides", None) or {}}], None)
            out_recs.extend(_expand_pooled(pos, getattr(r, "caap_group", "US"), hyp, pooled))
        expanded.append((cyc, out_recs))
    return expanded
