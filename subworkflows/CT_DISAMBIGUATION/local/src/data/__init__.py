"""
Data Layer Package

ReST-style Package Documentation
=================================

Purpose
--------
Provides unified data models and loaders for CAAS (Convergent Amino Acid Sites) analysis.
Handles metadata parsing, species/trait definitions, position indexing (zero/one-based),
taxid mapping, and validation with comprehensive logging.

Export Categories
-------------------

**Models**:
  - CAASPosition: CAAS site metadata (tag, position, amino acids, significance, p-value)
  - BiochemResults: Comprehensive convergence/biochemistry results object
  - ContrastDefinition: Species pair definition with taxid mapping

**Loaders**:
  - read_caas_metadata_table: Load CAAS metadata for specified positions into CAASPosition objects
  - list_gene_caas_positions: Raw DataFrame loading with flexible separator detection
  - get_caas_position_info: Fetch single CAAS position metadata by gene/position
  - parse_trait_pairs: Unified tab-separated trait parser (deterministic ordering)

Integration Points
-------------------
- **Upstream**: Consumed by single_gene_pipeline, aggregation_main, convergence analysis
- **Downstream**: Exports CAASPosition/BiochemResults to convergence & biochem modules
- **Dependencies**: pandas (DataFrame), pathlib (Path), csv (parsing), logging

Author: Refactored per biochem recipe | Date: 2025-12
"""

from .models import CAASPosition, BiochemResults, ContrastDefinition
from .loaders import (
    read_caas_metadata_table,
    build_caas_positions_map,
    list_gene_caas_positions,
    get_caas_position_info,
    parse_trait_pairs,
    load_ensembl_genes
)

__all__ = [
    # Models
    "CAASPosition",
    "BiochemResults",
    "ContrastDefinition",
    # Loaders
    "read_caas_metadata_table",
    "build_caas_positions_map",
    "list_gene_caas_positions",
    "get_caas_position_info",
    "parse_trait_pairs",
    "load_ensembl_genes"
]
