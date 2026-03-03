from pathlib import Path
import os
import json
import logging
from typing import Any, Dict, Optional, Tuple
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

def load_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    # Try tab-separated first
    try:
        df = pd.read_csv(path, sep="\t")
        # Verify it worked by checking if multiple columns exist
        if len(df.columns) > 1:
            return df
    except Exception:
        pass
    # Fall back to comma-separated
    return pd.read_csv(path)


def build_result_from_row(row: pd.Series) -> Optional[Dict[str, Any]]:
    """Construct a minimal ASR result entry from an aggregated CSV row."""
    # Accept either 'msa_pos' (0-based) or 'position' (1-based); normalize both
    pos0 = None
    pos1 = None
    if "msa_pos" in row and not pd.isna(row["msa_pos"]):
        try:
            pos0 = int(row["msa_pos"])
            pos1 = pos0 + 1
        except Exception:
            return None
    elif "position" in row and not pd.isna(row["position"]):
        try:
            pos1 = int(row["position"])
            pos0 = pos1 - 1
        except Exception:
            return None
    else:
        return None

    node_mapping: Dict[str, Any] = {}
    all_mrca_node = row.get("all_mrca_node") if hasattr(row, "get") else None
    if all_mrca_node is not None and not pd.isna(all_mrca_node):
        try:
            node_mapping["mrca_contrast"] = int(all_mrca_node)
        except Exception:
            pass

    focal_nodes: list[int] = []
    idx = 1
    while f"mrca_{idx}_node" in row:
        node_val = row.get(f"mrca_{idx}_node")
        if node_val is not None and not pd.isna(node_val):
            try:
                node_id = int(node_val)
                focal_nodes.append(node_id)
                node_mapping[f"focal_{idx}"] = node_id
            except Exception:
                pass
        idx += 1
    if focal_nodes:
        node_mapping["focal_nodes"] = focal_nodes

    node_state_details: Dict[str, Any] = {}
    all_mrca_state = row.get("all_mrca_state") if hasattr(row, "get") else None
    if all_mrca_state is not None and not pd.isna(all_mrca_state):
        node_state_details["mrca_contrast"] = str(all_mrca_state)
        prob = row.get("all_mrca_posterior")
        if prob is not None and not pd.isna(prob):
            try:
                node_state_details["mrca_contrast_prob"] = float(prob)
            except Exception:
                pass

    focal_states: list[str] = []
    focal_probs: list[float] = []
    for idx, node_id in enumerate(focal_nodes, start=1):
        state = row.get(f"mrca_{idx}_state") if hasattr(row, "get") else None
        prob = row.get(f"mrca_{idx}_posterior") if hasattr(row, "get") else None
        if state is not None and not pd.isna(state):
            focal_states.append(str(state))
        if prob is not None and not pd.isna(prob):
            try:
                focal_probs.append(float(prob))
            except Exception:
                pass

    if focal_states:
        node_state_details["focal_states"] = focal_states
    if focal_probs:
        node_state_details["focal_probs"] = focal_probs

    if not node_mapping or not node_state_details:
        return None

    return {
        # `position` is canonical 0-based index used across plotting code
        "position": pos0,
        # `position_one_based` is PAML-compatible 1-based index
        "position_one_based": pos1,
        "node_mapping": node_mapping,
        "node_state_details": node_state_details,
        "pair_details": [],  # Tip-level residues unavailable in bulk outputs
    }


def find_tree_file(gene: str, asr_root: Optional[Path]) -> Optional[Path]:
    """Locate the PAML tree file for a gene (tree_paml.nwk).

    If `asr_root` is None or standard candidates fail, this function will attempt
    broader autodetection by checking common `results` locations and scanning the
    repository for `tree_paml.nwk` files that appear associated with `gene`.
    """
    gene_caps = gene.upper()
    candidates = []

    # Common layouts observed in different pipelines:
    # 1) <asr_root>/asr_GENE/tree_paml.nwk
    # 2) <asr_root>/GENE/asr_GENE/tree_paml.nwk
    # 3) <asr_root>/GENE/tree_paml.nwk
    # 4) <asr_root>/asr_GENE/tree_paml.nwk (already covered)
    if asr_root is not None:
        candidates.append(asr_root / f"asr_{gene_caps}" / "tree_paml.nwk")
        candidates.append(asr_root / gene_caps / f"asr_{gene_caps}" / "tree_paml.nwk")
        candidates.append(asr_root / gene_caps / "tree_paml.nwk")

        # Also check lower/upper variations and direct children named with gene prefix
        # (e.g., asr_root / 'asr_BRCA2')
        candidates.append(asr_root / f"asr_{gene_caps}" / "tree_paml.nwk")
        if asr_root is not None:
            candidates.append(
                asr_root / gene_caps.lower() / f"asr_{gene_caps}" / "tree_paml.nwk"
            )

    for candidate in candidates:
        try:
            if candidate.exists():
                return candidate
        except Exception:
            continue

    # Fallback: scan immediate children of asr_root for directories containing the gene name
    try:
        if asr_root and asr_root.exists() and asr_root.is_dir():
            for child in asr_root.iterdir():
                name = child.name.lower()
                if gene.lower() in name or f"asr_{gene.lower()}" in name:
                    cand = child / "tree_paml.nwk"
                    if cand.exists():
                        return cand
    except Exception:
        pass

    # Check ASR cache env var first (precomputed cache)
    common_roots = []
    asr_cache_env = os.environ.get("ASR_CACHE_DIR")
    if asr_cache_env:
        try:
            common_roots.append(Path(asr_cache_env))
        except Exception:
            pass

    # Broader autodetection: check common results locations relative to cwd
    common_roots.extend(
        [
            Path.cwd() / "results" / "asr",
            Path.cwd() / "results_toy" / "asr",
            Path.cwd() / "results",
        ]
    )
    for root in common_roots:
        try:
            if not root.exists():
                continue
            for child in root.rglob("tree_paml.nwk"):
                # prefer files whose path contains the gene name or asr_ prefix
                parts = "/".join(p.name.lower() for p in child.parts)
                if gene.lower() in parts or f"asr_{gene.lower()}" in parts:
                    return child
        except Exception:
            continue

    # Last resort: scan the whole repo/workdir for any tree file that mentions the gene
    try:
        for candidate in Path.cwd().rglob("tree_paml.nwk"):
            parts = "/".join(p.lower() for p in candidate.parts)
            if gene.lower() in parts or f"asr_{gene.lower()}" in parts:
                return candidate
    except Exception:
        pass

    return None