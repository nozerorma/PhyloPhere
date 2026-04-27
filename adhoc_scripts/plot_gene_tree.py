#!/usr/bin/env python3
"""
Simple wrapper to generate gene tree plots using the disambiguation module.
Directly uses the actual disambiguation plotter functions.
"""

import sys
import json
import argparse
import logging
from pathlib import Path

# Setup paths
WORKSPACE_ROOT = Path(__file__).parent.parent
SUBWORKFLOW_ROOT = WORKSPACE_ROOT / "subworkflows" / "CT_DISAMBIGUATION" / "local"

# Add to path
sys.path.insert(0, str(SUBWORKFLOW_ROOT))

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Generate gene tree diagrams using the disambiguation module",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python plot_gene_tree.py MLPH 200
  python plot_gene_tree.py MLPH 365 --output /tmp/trees
        """
    )
    
    parser.add_argument("gene", help="Gene name (e.g., MLPH)")
    parser.add_argument("position", type=int, help="Position in alignment (0-based)")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output directory for plot (default: adhoc_scripts/output)"
    )
    parser.add_argument(
        "--debug-root",
        type=Path,
        default=WORKSPACE_ROOT / "debug" / "mlph",
        help="Debug root directory"
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=None,
        help="Tree file (auto-detected if not provided)"
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=20,
        help="Figure width in inches (default: 20)"
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=15,
        help="Figure height in inches (default: 15)"
    )
    
    args = parser.parse_args()
    
    # Import here after path setup
    import matplotlib.pyplot as plt
    from src.plots.gene_trees_bulk import create_gene_tree_state_plot
    
    # Monkey-patch plt.subplots to use custom figure sizes
    _orig_subplots = plt.subplots
    
    def _patched_subplots(*pargs, **kwargs):
        if 'figsize' in kwargs:
            # Replace calculated size with user settings
            kwargs['figsize'] = (args.fig_width, args.fig_height)
        return _orig_subplots(*pargs, **kwargs)
    
    plt.subplots = _patched_subplots
    
    output_dir = args.output or (WORKSPACE_ROOT / "adhoc_scripts" / "output")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load convergence JSONL
    jsonl_path = args.debug_root / "json_summaries" / f"{args.gene.lower()}_convergence_positions.jsonl"
    
    if not jsonl_path.exists():
        logger.error(f"JSONL not found: {jsonl_path}")
        return 1
    
    results = []
    with open(jsonl_path) as f:
        for line in f:
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                if record.get("position") == args.position:
                    results.append(record)
            except json.JSONDecodeError:
                continue
    
    if not results:
        logger.error(f"No records found for {args.gene}:{args.position}")
        return 1
    
    logger.info(f"Found {len(results)} record(s)")
    
    # Load node posteriors
    posteriors_path = args.debug_root / "diagnostics" / "node_dumps" / f"{args.gene.lower()}_posteriors.jsonl"
    per_node = {}
    if posteriors_path.exists():
        with open(posteriors_path) as f:
            for line in f:
                if not line.strip():
                    continue
                try:
                    obj = json.loads(line)
                    if obj.get("__metadata__"):
                        continue
                    nid = obj.get("node_id")
                    positions = obj.get("positions") or {}
                    if nid is not None:
                        nid = int(nid)
                        for pos_s, aa_map in positions.items():
                            try:
                                pos = int(pos_s)
                            except:
                                continue
                            if pos == args.position + 1:  # posteriors are 1-based, input is 0-based
                                if isinstance(aa_map, dict):
                                    best_aa = None
                                    best_prob = -1.0
                                    for aa, prob in aa_map.items():
                                        try:
                                            prob = float(prob)
                                        except:
                                            continue
                                        if prob > best_prob:
                                            best_prob = prob
                                            best_aa = aa
                                    if best_aa:
                                        per_node[nid] = (best_aa, best_prob)
                                        per_node[str(nid)] = (best_aa, best_prob)
                except json.JSONDecodeError:
                    continue
        
        if per_node:
            logger.info(f"Loaded posteriors for {len(per_node)} nodes")
            for result in results:
                result["node_posteriors"] = {"per_node": per_node}
    
    # Load tip details
    tip_path = args.debug_root / "diagnostics" / "tip_details" / f"{args.gene.lower()}_tip_details.jsonl"
    if tip_path.exists():
        with open(tip_path) as f:
            for line in f:
                if not line.strip():
                    continue
                try:
                    obj = json.loads(line)
                    if obj.get("position") == args.position:
                        for result in results:
                            result["tip_details"] = obj.get("pair_details", [])
                except json.JSONDecodeError:
                    continue
    
    # Auto-detect tree
    tree_file = args.tree
    if not tree_file:
        candidates = [
            WORKSPACE_ROOT / "Data" / "5.Phylogeny" / "versioned-trees" / "science.abn7829_data_s4.nex.fixed.tree",
            WORKSPACE_ROOT / "Data" / "5.Phylogeny" / "versioned-trees" / "science.abn7829_data_s4.nex.tree",
            WORKSPACE_ROOT / "Data" / "5.Phylogeny" / "science.abn7829_data_s4.nex.pruned.tree",
        ]
        for p in candidates:
            if p.exists():
                tree_file = p
                break
    
    if not tree_file or not Path(tree_file).exists():
        logger.error(f"Tree file not found: {tree_file}")
        return 1
    
    logger.info(f"Using tree: {tree_file}")
    
    # Create temporary RST symlink so plotter can find it
    rst_file = args.debug_root / "diagnostics" / "asr" / f"asr_{args.gene}" / "rst"
    rst_dest = Path(tree_file).parent / "rst"
    
    symlink_created = False
    try:
        if rst_file.exists() and not rst_dest.exists():
            logger.info(f"Creating RST symlink for plotter")
            rst_dest.symlink_to(rst_file)
            symlink_created = True
        
        # Call the actual disambiguation plotter
        plot_path = create_gene_tree_state_plot(
            results=results,
            output_dir=output_dir,
            gene_name=args.gene,
            tree_file=Path(tree_file),
            focus_position=args.position,
        )
        
        if not plot_path:
            logger.error("Failed to generate plot")
            return 1
        
        logger.info(f"✓ Plot saved: {plot_path}")
        print(f"\n✓ Success! Plot saved to:\n  {plot_path}\n")
        return 0
    
    finally:
        if symlink_created and rst_dest.is_symlink():
            rst_dest.unlink()
            logger.info("Cleaned up RST symlink")
        plt.subplots = _orig_subplots


if __name__ == "__main__":
    sys.exit(main())
