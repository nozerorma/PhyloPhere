# amino.py — Normalization of amino-acid code lists.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Amino-acid helpers shared by the disambiguation layers.

The module sits outside the ``data`` and ``convergence`` packages so that both can
depend on it without a circular import.

Imported by: src/core/observed.py, src/utils/gene_wrapper.py
"""

from typing import List, Optional


STANDARD_AAS = frozenset('ACDEFGHIKLMNPQRSTVWY')


def clean_amino(value):
    aa = str(value).strip().upper() if value is not None else ''
    return aa if aa in STANDARD_AAS else None


def normalize_amino_list(values: List[Optional[str]]) -> List[str]:
    return list(dict.fromkeys(aa for value in values or [] if (aa := clean_amino(value)) is not None))
