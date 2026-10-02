"""The observed labeling (b_0) as full records: discovery rows -> CAAS entries -> pooled results -> master rows.

The null scores its labelings as thin axes records; the observed labeling needs everything the master carries
(which hypotheses took part, the support tallies, the per-domain residues), so it goes through
`analyze_gene_disambiguation` in full mode, ONE call per gene over all of that gene's hypothesis rows: the call
pools the hypotheses of a position and scheme itself, with the same reducer (`fop_pool.pool_domains`) the null's
`pool_labelings` uses.

The rows come from the discovery of b_0: one row per (position, hypothesis, caap_group), with `trait` naming the
hypothesis ('H3', 'traitfile_H3.tab' or 'b_0~H3'). The order of the entries decides the order of the output rows
within a position and breaks the ties of the modal residues (the first residue seen wins), so a producer should
present them in the order discovery.tab has (position, then hypothesis, then scheme). Nothing else depends on it.
"""
from typing import Any, Dict, Iterable, List, Optional, Tuple

from src.convergence.disambiguate_single import analyze_gene_disambiguation
from src.core.meta import caas_id
from src.data.loaders import _parse_conserved_pair, as_bool, normalize_amino_list
from src.data.models import CAASPosition


def observed_entries(gene: str, rows: Iterable[Dict[str, Any]]) -> List[CAASPosition]:
    """CAAS entries of one gene from its discovery rows (mappings with position, trait, caap_group, caas,
    amino_encoded, pattern and, when present, is_conserved_meta and conserved_pair), in the order given.
    The tag of an entry is the content id of its row (`core.meta.caas_id`)."""
    entries: List[CAASPosition] = []
    for row in rows:
        pos0 = int(row["position"])
        caas = str(row.get("caas") or "")
        parts = caas.split("/") if caas else []
        group = str(row.get("caap_group") or "US")
        amino = str(row.get("amino_encoded") or "")
        trait = str(row.get("trait") or "")
        entries.append(CAASPosition(
            position=pos0,
            position_one_based=pos0 + 1,
            tag=caas_id(gene, pos0, trait, group, caas, amino, row.get("pattern", "")),
            caas=caas,
            trait1_aa=normalize_amino_list(list(parts[0])) if len(parts) == 2 else [],
            trait0_aa=normalize_amino_list(list(parts[1])) if len(parts) == 2 else [],
            caap_group=group,
            amino_encoded=amino,
            is_conserved_meta=as_bool(row.get("is_conserved_meta")),
            conserved_pair=_parse_conserved_pair(str(row.get("conserved_pair") or "")),
            trait=trait,
        ))
    return entries


def score_observed(
    ctx: Dict[str, Any],
    gene: str,
    entries: List[CAASPosition],
    trait_pairs: Dict[int, List[Tuple[str, str]]],
    hyp_pairs_pss: Optional[Dict[Tuple[str, int], float]],
    posterior_threshold: float,
) -> List[Any]:
    """Pooled full-mode ConvergenceResults of one gene. `ctx` is core.driver.load_gene_context's return;
    `trait_pairs` is core.labelings.read_trait_pairs of the design and `hyp_pairs_pss` its observed_pss."""
    node_posteriors = ctx["node_posteriors"]
    results, _ = analyze_gene_disambiguation(
        gene=gene,
        alignment_data=ctx["alignment_data"],
        tree_data=ctx["tree_data"],
        caas_positions=sorted({e.position for e in entries}),
        caas_entries=entries,
        trait_pairs=trait_pairs,
        taxid_mapping=ctx["alignment_data"].species_to_taxid,
        posterior_data=node_posteriors.posteriors_node if node_posteriors else None,
        posterior_threshold=posterior_threshold,
        hyp_pairs_pss=hyp_pairs_pss,
    )
    return results


def observed_master_rows(gene: str, results: Iterable[Any], master_fields: List[str]) -> List[Tuple[str, Any, Dict[str, str]]]:
    """(gene, msa_pos, row) triples for core.master.write_master_csv."""
    # imported here: gene_wrapper imports the core, so a top-level import would be circular
    from src.core.master import master_row
    from src.utils.gene_wrapper import convert_convergence_result_to_dict

    rows = []
    for r in results:
        record = convert_convergence_result_to_dict(r, multi_hypothesis=None)
        rows.append((gene, getattr(r, "position", None), master_row(record, master_fields)))
    return rows
