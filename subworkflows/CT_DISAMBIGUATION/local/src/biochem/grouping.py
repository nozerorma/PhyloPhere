# grouping.py — Amino acid grouping schemes (US, GS1 to GS4) used to encode residues.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/biochem/

"""
Amino-acid grouping schemes used for convergence typing: US, GS1, GS2, GS3, GS4.

The tables are identical to the ones of the discovery step in
``subworkflows/CT/local/modules/caas_id.py``, so a residue is encoded the same way when a CAAS is
found and when it is scored. Each scheme maps a residue to a group label; the comment above a
table lists its groups.

Imported by: src/convergence/path_scores.py, src/convergence/disambiguate_single.py
Inputs: none (constant tables)
Outputs: `get_grouping_scheme(aa, scheme)`, the group label of a residue
"""

from typing import Dict, Optional

# US: identity mapping (each residue is its own group)
US: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "Y": "Y",
}

# GS1: CV // AGPS // NDQE // RHK // ILMFWY // T
GS1: Dict[str, str] = {
    "C": "t",
    "V": "t",
    "A": "n",
    "G": "n",
    "P": "n",
    "S": "n",
    "N": "p",
    "D": "p",
    "Q": "p",
    "E": "p",
    "R": "b",
    "H": "b",
    "K": "b",
    "I": "h",
    "L": "h",
    "M": "h",
    "F": "h",
    "W": "h",
    "Y": "h",
    "T": "o",
}

# GS2: C // AGV // DE // NQHW // RK // ILFP // YMTS
GS2: Dict[str, str] = {
    "C": "c",
    "A": "s",
    "G": "s",
    "V": "s",
    "D": "a",
    "E": "a",
    "N": "n",
    "Q": "n",
    "H": "n",
    "W": "n",
    "R": "b",
    "K": "b",
    "I": "h",
    "L": "h",
    "F": "h",
    "P": "h",
    "Y": "x",
    "M": "x",
    "T": "x",
    "S": "x",
}

# GS3: C // AGPST // NDQE // RHK // ILMV // FWY
GS3: Dict[str, str] = {
    "C": "c",
    "A": "n",
    "G": "n",
    "P": "n",
    "S": "n",
    "T": "n",
    "N": "s",
    "D": "s",
    "Q": "s",
    "E": "s",
    "R": "b",
    "H": "b",
    "K": "b",
    "I": "l",
    "L": "l",
    "M": "l",
    "V": "l",
    "F": "g",
    "W": "g",
    "Y": "g",
}

# GS4: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
GS4: Dict[str, str] = {
    "C": "c",
    "A": "h",
    "I": "h",
    "L": "h",
    "V": "h",
    "S": "o",
    "T": "o",
    "N": "p",
    "Q": "p",
    "D": "a",
    "E": "a",
    "R": "b",
    "H": "b",
    "G": "g",
    "P": "r",
    "K": "k",
    "M": "m",
    "F": "f",
    "W": "y",
    "Y": "y",
}

SCHEMES: Dict[str, Dict[str, str]] = {
    "US": US,
    "GS1": GS1,
    "GS2": GS2,
    "GS3": GS3,
    "GS4": GS4,
}


def get_grouping_scheme(aa: str, scheme: str) -> Optional[str]:
    """Return the group label for an amino acid under a GS scheme."""
    aa_u = (aa or "").strip().upper()
    table = SCHEMES.get((scheme or "").strip().upper())
    if not table or not aa_u:
        return None
    return table.get(aa_u)
