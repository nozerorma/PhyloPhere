from pathlib import Path
import os
import json
import logging
from typing import Any, Dict, Optional, Tuple
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def chr_sort_key(val: str) -> Tuple[int, str]:
    s = val.lower()
    if s.startswith("chr"):
        s = s[3:]
    special = {"x": 23, "xx": 23, "y": 24, "yy": 24, "m": 25, "mt": 25, "mito": 25}
    if s in special:
        return (special[s], s)
    if s.isdigit():
        return (int(s), "")
    return (1000, s)


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


def safe_save(fig, path: Path):
    """
    Save matplotlib figure with parent directory creation.

    Args:
        fig: Matplotlib Figure object
        path: Output file path (.png assumed)
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    logger.info("✓ Saved plot: %s", path)


def load_gene_results(gene: str, json_root: Path) -> Optional[list]:
    """Load per-gene convergence summary JSON produced by single-gene pipeline."""
    candidates = [
        json_root / f"{gene.lower()}_convergence_summary.json",
        json_root / f"{gene.upper()}_convergence_summary.json",
    ]
    for path in candidates:
        if not path.exists():
            continue
        try:
            with path.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
            results = []
            if isinstance(payload, dict):
                if "results" in payload:
                    results = payload.get("results") or []
                elif "positions" in payload:
                    results = payload.get("positions") or []
            elif isinstance(payload, list):
                results = payload

            # Convert pairs array to pair_details format for compatibility with plots.py
            for result in results:
                if "pairs" in result and "pair_details" not in result:
                    # Convert pairs array to pair_details list format
                    result["pair_details"] = result.get("pairs", [])

            return results if results else []
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("Failed to load JSON summary for %s: %s", gene, exc)
    return None


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
            candidates.append(asr_root / gene_caps.lower() / f"asr_{gene_caps}" / "tree_paml.nwk")

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
    common_roots.extend([
        Path.cwd() / "results" / "asr",
        Path.cwd() / "results_toy" / "asr",
        Path.cwd() / "results",
    ])
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


def compute_mean_evidence_score(
    group: pd.DataFrame,
    p_floor: float = 1e-300,
) -> float:
    """
    Composite evidence score for visuals:
        E = mean( -log10(p_hyper), -log10(p_boot) )
    where each term is computed as the mean -log10(p) within the group
    for that method.

    - Uses first-available synonym column for each method to avoid double counting.
    - Floors p-values to avoid inf/-inf.
    - Returns NaN if no usable p-value columns exist.
    """

    if group is None or len(group) == 0:
        return float("nan")

    def first_existing(cols):
        for c in cols:
            if c in group.columns:
                return c
        return None

    def mean_log_evidence(col):
        s = pd.to_numeric(group[col], errors="coerce").astype(float)
        s = s.dropna()
        if s.empty:
            return None
        s = s.clip(lower=p_floor)
        return float((-np.log10(s)).mean())

    # Hypergeometric-like p-value column synonyms
    hyper_col = first_existing(("pvalue", "p_value"))

    # Bootstrap p-value column synonyms
    boot_col = first_existing(("pvalue.boot", "pvalue_boot"))

    evidences = []
    if hyper_col:
        e = mean_log_evidence(hyper_col)
        if e is not None and not np.isnan(e):
            evidences.append(e)

    if boot_col:
        e = mean_log_evidence(boot_col)
        if e is not None and not np.isnan(e):
            evidences.append(e)

    if not evidences:
        return float("nan")

    return float(np.mean(evidences))



def get_pattern_abbreviation(pattern: str) -> str:
    """Map pattern type to abbreviated label."""
    abbrev_map = {
        "convergent": "CON",
        "divergent": "DIV",
        "parallel": "PAR_ALL",
        "parallel_convergence": "PAR_CON",
        "parallel_divergence": "PAR_DIV",
        "parallel_mixed": "PAR_MIX",
        "parallel_codivergent": "PAR_CDV",
        "codivergent_top_converges": "CDT",
        "codivergent_bottom_converges": "CDB",
        "codivergent_top": "CDT",
        "codivergent_bottom": "CDB",
        "codivergent": "CDV",
        "no_change": "N/C",
        "ambiguous": "AMB",
    }
    return abbrev_map.get(pattern.lower(), pattern.upper()[:3])



def get_pattern_colors(pattern_type: str) -> str:
    """
    Get color for pattern type.

    Color scheme:
    - Red family: codivergent_top (top converges)
    - Green family: codivergent_bottom (bottom converges)
    - Purple family: parallel variants (convergence, divergence, mixed, codivergent)
    - Gray: other patterns
    """
    color_map = {
        "convergent": "#9467bd",  # Purple (simple convergent)
        "divergent": "#8c564b",  # Brown (simple divergent)
        "codivergent_top": "#e74c3c",  # Red (top converges)
        "codivergent_bottom": "#27ae60",  # Green (bottom converges)
        "codivergent": "#afa86b",  # Olive
        "parallel_convergence": "#da6fd6",  # Light purple
        "parallel_divergence": "#d88372",  # Light brown
        "parallel_mixed": "#70b4db",  # Sky blue
        "parallel_codivergent": "#e2d985",  # Light olive
        "no_change": "#727272",  # Light gray
        "ambiguous": "#cccccc",  # Medium gray
    }
    return color_map.get(pattern_type.lower(), "#cccccc")
