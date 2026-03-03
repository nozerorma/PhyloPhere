#!/usr/bin/env python3

from pathlib import Path
import logging
import pandas as pd
from typing import Any, Dict, List, Optional, Tuple
from collections import defaultdict
import json

import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend
import matplotlib.pyplot as plt

from src.plots.plot_utils import build_result_from_row, find_tree_file
from src.asr.tree_parser import build_node_mapping

logger = logging.getLogger(__name__)

HAS_MATPLOTLIB = True
HAS_SEABORN = True


def plot_random_gene_trees(
    df: pd.DataFrame,
    output_dir: Path,
    asr_root: Path,
    node_dumps_root: Path,
    tip_details_root: Path,
    n: int = 5,
):
    """Render gene tree plots for N random positions using existing ASR outputs."""
    if df.empty:
        logger.warning("No rows available for gene tree plotting")
        return

    # Check for required columns
    if "gene" not in df.columns:
        logger.warning("No 'gene' column in dataframe; skipping gene tree plots")
        return
    if "position" not in df.columns and "msa_pos" not in df.columns:
        logger.warning(
            "No 'position' or 'msa_pos' column in dataframe; skipping gene tree plots"
        )
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Cache per-gene result entries assembled from CSV rows
    gene_results_cache: Dict[str, list] = {}
    gene_posteriors_cache: Dict[str, Optional[Dict]] = {}

    def _load_node_posteriors(
        gene: str, focus_pos: Optional[int] = None, explicit_path: Optional[Path] = None
    ) -> Optional[Dict[str, Any]]:
        """Load per-node posterior states from JSONL for a gene.

        JSONL format expected: one JSON object per line: {"node_id": int, "positions": {"1": {"A":0.99}}}
        Optional metadata line: {"__metadata__": {"node_id_map": {...}}}
        """
        dump_path = Path(node_dumps_root) / f"{gene.lower()}_posteriors.jsonl"
        # Allow explicit per-gene dump path override (e.g., from DB alignment_extras)
        if explicit_path:
            dump_path = explicit_path
        if not dump_path.exists():
            return None
        # tip_details may be JSON or JSONL; try JSON first then JSONL
        tip_jsonl_path = Path(tip_details_root) / f"{gene.lower()}_tip_details.jsonl"
        if not tip_jsonl_path.exists():
            return None

        import json

        per_node: Dict = {}
        try:
            with open(dump_path, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    try:
                        obj = json.loads(line)
                    except Exception:
                        continue
                    # skip optional metadata
                    if obj.get("__metadata__"):
                        continue
                    nid = obj.get("node_id")
                    positions = obj.get("positions") or {}
                    try:
                        nid = int(nid)
                    except Exception:
                        continue
                    for pos_s, aa_map in positions.items():
                        try:
                            pos = int(pos_s)
                        except Exception:
                            continue
                        if focus_pos is not None and pos != focus_pos:
                            continue
                        # pick modal amino acid for this node/pos from aa_map
                        if not isinstance(aa_map, dict):
                            continue
                        best_aa = None
                        best_prob = -1.0
                        for aa, prob in aa_map.items():
                            try:
                                prob = float(prob)
                            except Exception:
                                continue
                            if prob > best_prob:
                                best_prob = prob
                                best_aa = aa
                        if best_aa is not None:
                            per_node[nid] = (best_aa, best_prob)
                            per_node[str(nid)] = (best_aa, best_prob)
            if per_node:
                logger.info(
                    f"Loaded poster node JSONL for {gene}: {len(per_node)} nodes (including string keys)"
                )
                return {"per_node": per_node}
            return None
        except Exception as err:
            logger.debug("Failed to load node posteriors for %s: %s", gene, err)
            return None

    saved = 0
    # Shuffle all rows to improve chances, stop after n saved
    for _, row in df.sample(frac=1, random_state=42).iterrows():
        if saved >= n:
            break

        gene = str(row["gene"])

        # Try both position and msa_pos columns
        pos = None
        if "position" in row and pd.notna(row["position"]):
            pos = int(row["position"])
        elif "msa_pos" in row and pd.notna(row["msa_pos"]):
            pos = int(row["msa_pos"])

        if pos is None:
            continue

        if gene not in gene_results_cache:
            results: list = []
            # Build minimal candidates from aggregated CSV rows
            for _, gene_row in df[df["gene"] == gene].iterrows():
                built = build_result_from_row(gene_row)
                if built:
                    results.append(built)

            gene_results_cache[gene] = results

        results = gene_results_cache.get(gene) or []
        if not results:
            continue

        tree_file = find_tree_file(gene, asr_root)
        if not tree_file:
            continue

        # Load node posteriors from JSONL and attach to results
        # NOTE: JSONL contains 1-based PAML positions, so convert from CAAS 0-based
        if gene not in gene_posteriors_cache:
            # Try default node_dumps_root first
            post = _load_node_posteriors(
                gene, focus_pos=pos + 1 if pos is not None else None
            )
            # If missing, try to find explicit path in results (e.g., prior DB export added posterior_dump_path)
            if not post:
                dump_path = None
                for r in results:
                    if isinstance(r, dict) and r.get("posterior_dump_jsonl"):
                        posterior_dump = r.get("posterior_dump_jsonl")
                        if posterior_dump is not None:
                            dump_path = Path(posterior_dump)
                        break
                if dump_path:
                    post = _load_node_posteriors(
                        gene,
                        focus_pos=pos + 1 if pos is not None else None,
                        explicit_path=dump_path,
                    )
            gene_posteriors_cache[gene] = post

        node_posteriors = gene_posteriors_cache[gene]
        if node_posteriors:
            for result in results:
                result["node_posteriors"] = node_posteriors

        # Load tip details from JSON and attach to results
        # NOTE: JSONL contains 0-based CAAS positions
        tip_jsonl_path = tip_details_root / f"{gene.lower()}_tip_details.jsonl"
        if tip_jsonl_path.exists():
            try:
                tip_details = []
                with open(tip_jsonl_path, "r", encoding="utf-8") as f:
                    for line in f:
                        if not line.strip():
                            continue
                        try:
                            obj = json.loads(line)
                            if obj.get("position") == pos:
                                tip_details.append(obj)
                        except Exception:
                            continue

                # Optional light pre-flattening to help downstream code
                flattened_pairs: List[Dict[str, Any]] = []
                if isinstance(tip_details, list):
                    for rec in tip_details:
                        if isinstance(rec, dict) and "pair_details" in rec:
                            flattened_pairs.extend(rec.get("pair_details") or [])
                elif isinstance(tip_details, dict):
                    flattened_pairs = tip_details.get("pair_details") or []

                for result in results:
                    result["tip_details"] = tip_details
                    if flattened_pairs:
                        result["pair_details"] = flattened_pairs

            except Exception as e:
                logger.debug(f"Failed to load tip details for {gene}: {e}")

        plot_path = create_gene_tree_state_plot(
            results=results,
            output_dir=output_dir,
            gene_name=gene,
            tree_file=tree_file,
            focus_position=pos,
        )
        if plot_path:
            saved += 1
            logger.info("✓ Saved gene tree plot for %s:%s -> %s", gene, pos, plot_path)

    if saved == 0:
        logger.warning(
            "Gene tree plotting skipped: no eligible genes with ASR details and tree files"
        )


def create_gene_tree_state_plot(
    results: List[Dict[str, Any]],
    output_dir: Path,
    gene_name: str,
    tree_file: Path,
    focus_position: Optional[int] = None,
) -> Optional[Path]:
    """
    Plot annotated gene tree for a specific position (or first available) using ASR states.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("Matplotlib not available, skipping gene tree plot")
        return None

    candidate = None
    if focus_position is not None:
        # First try to find a result with tip_details (JSON-loaded results are preferred)
        for r in results:
            pos0 = r.get("position")
            pos1 = r.get("position_one_based")
            if (
                (
                    (pos0 is not None and pos0 == focus_position)
                    or (
                        pos1 is not None
                        and pos0 is not None
                        and pos1 == focus_position + 1
                    )
                )
                and r.get("node_mapping")
                and r.get("node_state_details")
                and r.get("tip_details")
            ):
                # Prefer results with tip_details (from JSON)
                candidate = r
                break

        if candidate is None:
            logger.info(
                f"No result entry matches focus_position={focus_position}; skipping gene tree plot"
            )
            return None
    else:
        logger.debug("No focus_position provided, selecting first available result")
        # Prefer JSON-based result (has tip_details)
        for r in results:
            if (
                r.get("node_mapping")
                and r.get("node_state_details")
                and r.get("tip_details")
            ):
                candidate = r
                break

        # Fallback to any result with node data
        if candidate is None:
            for r in results:
                if r.get("node_mapping") and r.get("node_state_details"):
                    candidate = r
                    break

    if not candidate or not tree_file or not tree_file.exists():
        return None

    # Use RST file for PAML node IDs (same as used in analysis)
    rst_file = tree_file.parent / "rst"
    if not rst_file.exists():
        logger.warning(f"RST file not found for gene tree plot: {rst_file}")
        return None

    try:
        # Build tree with PAML node IDs from RST (consistent with analysis)
        ordered_nodes, id_mapping = build_node_mapping(rst_file=rst_file)
        root = ordered_nodes[-1]  # Last node in postorder is root
    except Exception as e:
        logger.warning(f"Could not build tree from RST for gene tree plot: {e}")
        return None

    # Use node object as key to handle leaves (which have node_id=None)
    positions_by_node: Dict[object, Tuple[float, float]] = {}
    leaf_counter = [0]

    def layout(node, depth=0):
        if node.is_leaf():
            y = leaf_counter[0]
            leaf_counter[0] += 1
        else:
            child_ys = [layout(child, depth + 1) for child in node.children]
            y = sum(child_ys) / len(child_ys) if child_ys else leaf_counter[0]
        positions_by_node[node] = (depth, y)
        return y

    layout(root, depth=0)

    scale_factor = 1.2
    for node, (x, y) in list(positions_by_node.items()):
        positions_by_node[node] = (x * scale_factor, y)

    # Build positions dict by node_id (only for nodes with IDs)
    positions: Dict[int, Tuple[float, float]] = {
        getattr(node, "node_id"): pos
        for node, pos in positions_by_node.items()
        if hasattr(node, "node_id") and getattr(node, "node_id") is not None
    }

    leaf_count = max(leaf_counter[0], 1)
    fig_height = min(15, max(5.5, leaf_count * 0.25))
    fig_width = min(14, max(8, leaf_count * 0.22 + 6))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    def draw_edges(node):
        x1, y1 = positions_by_node[node]
        for child in node.children:
            x2, y2 = positions_by_node[child]
            ax.plot([x1, x2], [y1, y2], color="#444444", linewidth=1)
            draw_edges(child)

    draw_edges(root)

    node_state_details = candidate.get("node_state_details", {})
    node_posteriors = candidate.get("node_posteriors") or {}
    per_node_posteriors = (
        node_posteriors.get("per_node") or node_posteriors.get("nodes") or {}
    )
    node_mapping = candidate.get("node_mapping", {})
    # Build mapping from node_id -> role; handle scalar and list (focal_nodes) values
    node_roles_by_id = {}
    focal_nodes = []
    for role, node_id in node_mapping.items():
        if isinstance(node_id, int):
            node_roles_by_id[node_id] = role
        elif isinstance(node_id, list):
            # Capture focal_nodes so we can annotate and attach states below
            if role == "focal_nodes":
                focal_nodes = node_id

    role_colors = {
        "root": "#2f1b9e",
        "mrca_contrast": "#be66b2",
        "focal_nodes": "#7570b3",
    }

    # Prepare tip residue lookup tables from both pairs (JSON) and tip_details (legacy)
    tip_state_map: Dict[str, str] = {}
    tip_group_map: Dict[str, str] = {}
    taxid_to_species: Dict[str, str] = {}  # Reverse lookup: taxid -> species name

    def _store_tip(
        key: Optional[str],
        residue: Optional[str],
        group_label: Optional[str],
        species: Optional[str] = None,
        taxid: Optional[str] = None,
    ):
        if not key:
            return
        normalized_key = str(key)
        normalized_key_lc = normalized_key.lower()
        if residue:
            tip_state_map.setdefault(normalized_key, residue)
            tip_state_map.setdefault(normalized_key_lc, residue)
        if group_label:
            tip_group_map.setdefault(normalized_key, group_label)
            tip_group_map.setdefault(normalized_key_lc, group_label)
        # Build taxid -> species reverse lookup
        if taxid and species:
            taxid_to_species.setdefault(str(taxid), species)
            taxid_to_species.setdefault(str(taxid).lower(), species)

    # --- Tip details normalization (handles wrapper-list JSON shape) ---
    tip_details_raw = candidate.get("tip_details") or []
    logger.debug(f"Processing tip_details data: {tip_details_raw}")

    pair_details = candidate.get("pair_details")

    # If pair_details not directly available, try to extract from tip_details_raw
    if not pair_details:
        if isinstance(tip_details_raw, dict):
            pair_details = tip_details_raw.get("pair_details") or []
        elif isinstance(tip_details_raw, list):
            # Case A: wrapper list with nested pair_details
            if (
                tip_details_raw
                and isinstance(tip_details_raw[0], dict)
                and "pair_details" in tip_details_raw[0]
            ):
                pair_details = []
                for rec in tip_details_raw:
                    pair_details.extend(rec.get("pair_details") or [])
            else:
                # Case B: already a list of pair dicts
                pair_details = tip_details_raw
        else:
            pair_details = []
    else:
        # Ensure it's a list if someone accidentally stored a dict
        if isinstance(pair_details, dict):
            pair_details = [pair_details]

    # Defensive: guarantee list
    if pair_details is None:
        pair_details = []

    # Ingest tip residues from normalized pair_details
    for pair in pair_details:
        for group_name, group_label in (
            ("top_tip_residues", "TOP"),
            ("bottom_tip_residues", "BOTTOM"),
        ):
            for rec in pair.get(group_name) or []:
                logger.debug(f"Processing tip record: {rec}")
                residue = rec.get("residue")
                species = rec.get("species")
                taxid = rec.get("taxid")
                label = rec.get("label")
                for key in filter(None, [taxid, species, label]):
                    _store_tip(key, residue, group_label, species=species, taxid=taxid)

    # Prepare per-node state summaries (supports dynamic focal_N lists)
    state_by_node: Dict[int, Tuple[str, Optional[float]]] = {}

    # Scalar roles (root, mrca_contrast, etc.)
    for role, node_id in node_mapping.items():
        if isinstance(node_id, int):
            state = node_state_details.get(role)
            if state:
                prob = node_state_details.get(f"{role}_prob")
                state_by_node[node_id] = (state, prob)

    # Dynamic focal nodes: align focal_nodes with focal_states/focal_probs lists
    focal_states = node_state_details.get("focal_states") or []
    focal_probs = node_state_details.get("focal_probs") or []
    for idx, node_id in enumerate(focal_nodes, start=1):
        if node_id is None:
            continue
        node_roles_by_id[node_id] = f"focal_{idx}"
        state = focal_states[idx - 1] if idx - 1 < len(focal_states) else None
        prob = focal_probs[idx - 1] if idx - 1 < len(focal_probs) else None
        if state:
            state_by_node[node_id] = (state, prob)

    # Fallback from pair_details for MRCA modal AA at the specific pair node
    for pair in pair_details:
        node_id = pair.get("node_id")
        modal = pair.get("mrca_modal_aa")
        if node_id and modal and node_id not in state_by_node:
            state_by_node[node_id] = (modal, None)

    # Ingest per_node_posteriors: handles both legacy dict {'state': aa, 'prob': prob}
    # and normalized tuple (aa, prob) formats
    if per_node_posteriors:
        for node_key, info in per_node_posteriors.items():
            try:
                node_id = int(node_key)
            except (TypeError, ValueError):
                continue

            # Handle both tuple format (aa, prob) and dict format {'state': aa, 'prob': prob}
            if isinstance(info, tuple) and len(info) >= 1:
                # Tuple format: (aa, prob)
                state = info[0]
                prob = info[1] if len(info) > 1 else None
            elif isinstance(info, dict):
                # Dict format: {'state': aa, 'prob': prob}
                state = info.get("aa") or info.get("modal_aa") or info.get("state")
                prob = info.get("prob")
            else:
                continue

            if state and node_id not in state_by_node:
                state_by_node[node_id] = (state, prob)

    pair_annotations = defaultdict(list)
    for pair in pair_details:
        node_id = pair.get("node_id")
        if node_id:
            pair_annotations[node_id].append(pair)

    # Annotate all nodes with node IDs and available states
    # Only prominently display nodes with roles or specific states
    node_label_offset = 0.45
    highlighted_count = 0
    for node_id, (x, y) in positions.items():
        role = node_roles_by_id.get(node_id)

        # Only draw and label nodes with roles or state information
        state_entry = state_by_node.get(node_id)
        if not role and not state_entry:
            # Draw small gray dot for unmarked internal nodes
            ax.scatter([x], [y], color="#d0d0d0", s=15, zorder=3, alpha=0.5)
            continue

        highlighted_count += 1
        color = role_colors.get(role, "#4c566a") if role else "#4c566a"
        marker_size = 120 if role else 50
        ax.scatter(
            [x],
            [y],
            color=color,
            s=marker_size,
            zorder=5,
            edgecolor="white",
            linewidth=1.0,
        )

        # Use PAML node ID directly (already set from build_node_mapping)
        label_lines = []
        if role:
            # Friendlier labels for MRCA/focal roles
            display_role = {"mrca_contrast": "MRCA_ALL"}.get(role, role)
            if role.startswith("focal_"):
                try:
                    idx = int(role.split("_")[1])
                    display_role = f"MRCA_{idx}"
                except Exception:
                    display_role = role
            label_lines.append(display_role)
        label_lines.append(f"Node {node_id}")

        if state_entry:
            state, prob = state_entry
            if prob is not None:
                label_lines.append(f"{state} ({prob:.2f})")
            else:
                label_lines.append(state)
        else:
            label_lines.append("?")

        label_text = "\n".join(label_lines)
        ax.text(
            x + node_label_offset,
            y + (0.2 if role else 0.05),
            label_text,
            va="center",
            ha="left",
            fontsize=8 if role else 7,
            fontweight="bold" if role else "normal",
            color=color,
            bbox=dict(
                boxstyle="round,pad=0.3",
                fc="white",
                ec=color,
                alpha=0.9,
                linewidth=1.5 if role else 1.0,
            ),
        )

    # Annotate tip states (residues derived from trait groups)
    max_x = max(coord[0] for coord in positions_by_node.values())

    def annotate_tips(node):
        if node.is_leaf():
            raw_label = node.name or str(node.node_id)
            label = str(raw_label)
            residue = None
            group = None

            # Extract taxid from PAML format: {node_id}_{taxid}
            taxid_from_label = None
            if "_" in label:
                parts = label.split("_")
                if len(parts) >= 2:
                    # Try last part as taxid
                    taxid_from_label = parts[-1]

            lookup_keys = [
                label,
                label.lower(),
                taxid_from_label,  # Try extracted taxid
                label.replace("_", " ") if label else None,
                label.replace(" ", "_") if label else None,
                label.lower().replace("_", " ") if label else None,
                label.lower().replace(" ", "_") if label else None,
            ]

            matched_key = None
            for key in lookup_keys:
                if key:
                    normalized_key = str(key)
                    if normalized_key in tip_state_map:
                        residue = tip_state_map[normalized_key]
                        matched_key = normalized_key
                    if normalized_key in tip_group_map:
                        group = tip_group_map.get(normalized_key, group)
                    if residue is not None or group is not None:
                        break

            if residue is None and group is None:
                return  # Skip tips outside the species of interest

            # Determine display label: prefer species name over taxid
            display_label = label
            if (
                matched_key
                and matched_key == taxid_from_label
                and taxid_from_label in taxid_to_species
            ):
                # We matched by taxid - use species name for display
                display_label = taxid_to_species[taxid_from_label]

            # Trait coloring: red for TOP (trait 1), green for BOTTOM (trait 0)
            color = (
                "#e74c3c"
                if group == "TOP"
                else "#27ae60" if group == "BOTTOM" else "#95a5a6"
            )
            x, y = positions_by_node[node]  # Use positions_by_node for leaves
            ax.scatter(
                [x],
                [y],
                color=color,
                s=80,
                zorder=4,
                edgecolor="white",
                linewidth=1.2,
                marker="o",
            )
            # Show species name and residue with trait color
            ax.text(
                max_x + 0.6,
                y,
                f"{display_label}: {residue or '?'}",
                va="center",
                ha="left",
                fontsize=7,
                fontweight="bold",
                color=color,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    fc="white",
                    ec=color,
                    alpha=0.9,
                    linewidth=1.5,
                ),
            )
        else:
            for child in node.children:
                annotate_tips(child)

    annotate_tips(root)

    position_label = candidate.get("position")
    focus_label = focus_position
    if position_label is None and focus_label is not None:
        position_label = focus_label
    suffix = f"_pos{focus_label}" if focus_label is not None else ""
    plot_path = output_dir / f"{gene_name.lower()}_gene_tree_states{suffix}.png"
    title = f"{gene_name} tree – position {position_label if position_label is not None else '?'} (0-based)"
    if focus_label is not None:
        title = f"{title} (focus {focus_label})"
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.axis("off")

    fig.tight_layout()
    fig.savefig(plot_path, bbox_inches="tight", dpi=300)
    if plt is not None:
        plt.close(fig)

    logger.info("✓ Saved gene tree state plot")
    return plot_path
