"""PhyloPhere validation harness — tier-agnostic building blocks.

Nothing here imports the pipeline runtime. Metrics can be computed from an
existing ``runs/`` output tree with only stdlib + numpy.
"""

from .phenotype import PhenotypeSpec, DatasetSpec, load_dataset_spec
from .trait_matrix import emit_trait_matrix

__all__ = [
    "PhenotypeSpec",
    "DatasetSpec",
    "load_dataset_spec",
    "emit_trait_matrix",
]
