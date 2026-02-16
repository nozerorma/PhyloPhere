"""Amino-acid grouping schemes used for convergence typing.

This module intentionally keeps only the minimal functionality required by the
active pipeline: group definitions (GS1-GS4) and amino-acid to group lookup.
"""

from typing import Dict, Optional

# GS1: CV // AGPS // NDQE // RHK // ILMFWY // T
GS1: Dict[str, str] = {
    "C": "CV",
    "V": "CV",
    "A": "AGPS",
    "G": "AGPS",
    "P": "AGPS",
    "S": "AGPS",
    "N": "NDQE",
    "D": "NDQE",
    "Q": "NDQE",
    "E": "NDQE",
    "R": "RHK",
    "H": "RHK",
    "K": "RHK",
    "I": "ILMFWY",
    "L": "ILMFWY",
    "M": "ILMFWY",
    "F": "ILMFWY",
    "W": "ILMFWY",
    "Y": "ILMFWY",
    "T": "T",
}

# GS2: C // AGV // DE // NQHW // RK // ILFP // YMTS
GS2: Dict[str, str] = {
    "C": "C",
    "A": "AGV",
    "G": "AGV",
    "V": "AGV",
    "D": "DE",
    "E": "DE",
    "N": "NQHW",
    "Q": "NQHW",
    "H": "NQHW",
    "W": "NQHW",
    "R": "RK",
    "K": "RK",
    "I": "ILFP",
    "L": "ILFP",
    "F": "ILFP",
    "P": "ILFP",
    "Y": "YMTS",
    "M": "YMTS",
    "T": "YMTS",
    "S": "YMTS",
}

# GS3: C // AGPST // NDQE // RHK // ILMV // FWY
GS3: Dict[str, str] = {
    "C": "C",
    "A": "AGPST",
    "G": "AGPST",
    "P": "AGPST",
    "S": "AGPST",
    "T": "AGPST",
    "N": "NDQE",
    "D": "NDQE",
    "Q": "NDQE",
    "E": "NDQE",
    "R": "RHK",
    "H": "RHK",
    "K": "RHK",
    "I": "ILMV",
    "L": "ILMV",
    "M": "ILMV",
    "V": "ILMV",
    "F": "FWY",
    "W": "FWY",
    "Y": "FWY",
}

# GS4: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
GS4: Dict[str, str] = {
    "C": "C",
    "A": "AILV",
    "I": "AILV",
    "L": "AILV",
    "V": "AILV",
    "S": "ST",
    "T": "ST",
    "N": "NQ",
    "Q": "NQ",
    "D": "DE",
    "E": "DE",
    "R": "RH",
    "H": "RH",
    "G": "G",
    "P": "P",
    "K": "K",
    "M": "M",
    "F": "F",
    "W": "WY",
    "Y": "WY",
}

SCHEMES: Dict[str, Dict[str, str]] = {
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
