"""CAAS ids: a function of the content of one discovery.tab row.

A CAAS is one row of discovery.tab: a gene, a 0-based alignment position, a hypothesis, a caap_group (scheme)
and the substitution found there (caas, amino_encoded, pattern). Its id hashes exactly those fields, so the same
CAAS gets the same id in every run, every file and every process, and it does not depend on which other rows exist
or on their order.

The id is `CAAS_` plus the first `ID_HEX_CHARS` hexadecimal characters of the SHA-256 of the fields. With 16
characters (64 bits) a table of 4.4 million rows has a collision probability of about 5e-7; `assign_ids` still
raises if two different rows share an id.
"""
import hashlib
from typing import Iterable, List, Tuple

from src.core.labelings import hyp_id

ID_PREFIX = "CAAS_"
ID_HEX_CHARS = 16
_SEP = "\x1f"  # unit separator: no field can contain it, so ("AB", "C") and ("A", "BC") hash differently


def caas_id(gene, position, hypothesis, caap_group, caas, amino_encoded, pattern, hex_chars: int = ID_HEX_CHARS) -> str:
    """Id of one CAAS row. `hypothesis` is anything `hyp_id` accepts ('H3', 'traitfile_H3.tab', 'b_0~H3')."""
    fields = (str(gene), str(int(position)), hyp_id(hypothesis), str(caap_group), str(caas), str(amino_encoded), str(pattern))
    digest = hashlib.sha256(_SEP.join(fields).encode("utf-8")).hexdigest()
    return ID_PREFIX + digest[:hex_chars].upper()


def row_id(gene, row) -> str:
    """Id of a discovery.tab row given as a mapping (position, trait, caap_group, caas, amino_encoded, pattern);
    an empty or missing field counts as the empty string, and a missing caap_group as 'US'."""
    pattern = row.get("pattern")
    return caas_id(gene, row["position"], row.get("trait") or "", row.get("caap_group") or "US",
                   row.get("caas") or "", row.get("amino_encoded") or "", "" if pattern is None else pattern)


def assign_ids(rows: Iterable[Tuple], hex_chars: int = ID_HEX_CHARS) -> List[str]:
    """Ids of (gene, position, hypothesis, caap_group, caas, amino_encoded, pattern) rows, in the order given.

    Identical rows share their id; two different rows with the same id raise ValueError.
    """
    ids: List[str] = []
    seen = {}
    for row in rows:
        key = (str(row[0]), int(row[1]), hyp_id(row[2])) + tuple(str(v) for v in row[3:7])
        cid = caas_id(*row[:7], hex_chars=hex_chars)
        if seen.setdefault(cid, key) != key:
            raise ValueError(f"CAAS id collision: {seen[cid]} and {key} both give {cid}")
        ids.append(cid)
    return ids
