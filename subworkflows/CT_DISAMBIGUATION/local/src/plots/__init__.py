"""
Plotting Module
================

Output generation and summary statistics.

Components:
    - plots: Visualization creation for analysis results
"""

# Plotting functions
from .gene_trees_bulk import (
    create_gene_tree_state_plot,
)
from .bulk_plots import (
    generate_bulk_plots,
)


# Export all available functions
__all__ = [
    # Plotting functions
    "create_gene_tree_state_plot",
    "generate_bulk_plots",
]
