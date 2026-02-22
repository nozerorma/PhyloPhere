"""Amino-acid grouping schemes used for convergence typing.

Definitions are aligned with the paired CAAP implementation in
``subworkflows/CT/local/modules/caap_id.py`` and include:
US, GS0, GS1, GS2, GS3, GS4.
"""

from typing import Dict, Optional

# US: identity mapping (classical CAAS)
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

# GS0: random aggregation control
# Partition: GWDC // PM // K // IQLS // EATVYF // NHR
GS0: Dict[str, str] = {
    "G": "a",
    "W": "a",
    "D": "a",
    "C": "a",
    "P": "b",
    "M": "b",
    "K": "c",
    "I": "d",
    "Q": "d",
    "L": "d",
    "S": "d",
    "E": "e",
    "A": "e",
    "T": "e",
    "V": "e",
    "Y": "e",
    "F": "e",
    "N": "f",
    "H": "f",
    "R": "f",
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
    "GS0": GS0,
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
