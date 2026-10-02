"""The master CSV of the observed labeling: rows are serialized once, in the worker, and written once.

caas_convergence_master.csv is the table the rest of the pipeline reads. It used to be produced by
re-reading every result from the aggregation SQLite database; the rows are now built from the records the
workers already hold and written here, in the order the database export used (gene, msa_pos, then the order
a gene's rows were produced). The database keeps only the decoration outputs.
"""
import csv
import gzip
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


def serialize_value(val: Any) -> str:
    """CSV cell of a record value: lists comma-joined, None empty, everything else str()."""
    if isinstance(val, (list, tuple)):
        return ",".join(str(v) for v in val)
    if val is None:
        return ""
    return str(val)


def master_row(record: Dict[str, Any], fields: Sequence[str]) -> Dict[str, str]:
    return {k: serialize_value(record.get(k)) for k in fields}


def write_master_csv(
    rows: Iterable[Tuple[str, Any, Dict[str, str]]], path: Path, fields: Sequence[str]
) -> int:
    """Write (gene, msa_pos, row) triples sorted by gene then msa_pos (stable: production order within
    a position is kept). Returns the number of rows written; an empty input writes just the header.
    A path ending in .gz is written compressed."""
    ordered = sorted(rows, key=lambda t: (t[0], -1 if t[1] is None else t[1]))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for _gene, _pos, row in ordered:
            writer.writerow(row)
    return len(ordered)
