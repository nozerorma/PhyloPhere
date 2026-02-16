#!/usr/bin/env python3
"""
SQLite helpers for disambiguation storage
=========================================

Small SQLite helpers to store alignment metadata and per-position per-hypothesis
result dictionaries (as JSON). This module is intentionally lightweight and
does not perform heavy biological normalization; consolidations and merging
should occur upstream (worker) or downstream (export writers).

Schema
------
The module maintains two tables: `gene_alignment` and `results`.

Usage Example
-------------
::

    from src.utils.disambiguation_db import init_db, get_connection, insert_gene_alignment

    db_path = Path('aggregation.sqlite3')
    init_db(db_path)
    conn = get_connection(db_path)
    insert_gene_alignment(conn, 'BRCA2', alignment_obj)

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date
----
2025-12-09
"""

from __future__ import annotations

import json
import logging
import sqlite3
import pickle
from pathlib import Path
from typing import Any, Dict, Iterator, Optional, Tuple

logger = logging.getLogger(__name__)


def _sanitize_for_json(obj: Any) -> Any:
    """Recursively convert common non-JSON types to JSON-safe values."""
    from pathlib import Path as _Path

    if obj is None:
        return None
    if isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, bytes):
        try:
            return obj.decode("utf-8")
        except Exception:
            return str(obj)
    if isinstance(obj, _Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_sanitize_for_json(v) for v in obj]
    try:
        return str(obj)
    except Exception:
        return None


def _load_posteriors_from_jsonl(jsonl_path: Path) -> Dict[int, Dict[int, Dict[str, float]]]:
    """Load posteriors from a JSON Lines (JSONL) file exported by export_posteriors_to_jsonl.

    Each line is expected to be a JSON object like::

        {"node_id": 234, "positions": {"1": {"A": 0.95, "R": 0.02}, ...}}

    :param jsonl_path: Path to the JSONL file to read.
    :type jsonl_path: Path
    :returns: Mapping from node_id -> position -> {AA: posterior}
    :rtype: Dict[int, Dict[int, Dict[str, float]]]
    :example: ::

        loaded = _load_posteriors_from_jsonl(Path('posteriors.jsonl'))
    """
    result: Dict[int, Dict[int, Dict[str, float]]] = {}
    p = Path(jsonl_path)
    if not p.exists():
        return result

    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue

            # Skip metadata lines if your exporter uses them
            if isinstance(obj, dict) and obj.get("__metadata__"):
                continue

            try:
                node_id = int(obj.get("node_id") or 0)
            except Exception:
                continue

            positions = obj.get("positions") or {}
            node_map = result.setdefault(node_id, {})

            if isinstance(positions, dict):
                for pos_s, aa_map in positions.items():
                    try:
                        pos = int(pos_s)
                    except Exception:
                        continue
                    if isinstance(aa_map, dict):
                        node_map[pos] = aa_map

    if result:
        logger.info(f"Loaded posteriors JSONL from {p} ({len(result)} nodes)")
    return result


def init_db(db_path: Path) -> None:
    """Create the SQLite database schema if needed.

    :param db_path: Path to the SQLite database file to initialize.
    :type db_path: Path
    :returns: None
    :rtype: None
    """
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("PRAGMA journal_mode=WAL;")
        cur.execute("PRAGMA synchronous=NORMAL;")

        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS gene_alignment (
                gene TEXT PRIMARY KEY,
                seq_by_id_json TEXT,
                seq_by_species_json TEXT,
                taxid_to_species_json TEXT,
                species_to_taxid_json TEXT,
                alignment_extras_json TEXT,
                alignment_path TEXT,
                num_sequences INTEGER,
                alignment_len INTEGER
            )
            """
        )

        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene TEXT,
                msa_pos INTEGER,
                position INTEGER,
                pair_count INTEGER,
                result_json TEXT
            )
            """
        )

        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_results_gene_msa ON results(gene, msa_pos);"
        )
        conn.commit()
    finally:
        conn.close()


def get_connection(db_path: Path) -> sqlite3.Connection:
    """Open and return a SQLite connection with WAL pragmas enabled.

    :param db_path: Path to the SQLite database file.
    :type db_path: Path
    :returns: SQLite connection object with appropriate pragmas set.
    :rtype: sqlite3.Connection
    """
    conn = sqlite3.connect(str(db_path), timeout=30)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    return conn


def insert_gene_alignment(conn: sqlite3.Connection, gene: str, alignment_obj: Dict[str, Any]) -> None:
    """Insert or replace per-gene alignment metadata in the `gene_alignment` table.

    :param conn: Open SQLite connection where the record will be stored.
    :type conn: sqlite3.Connection
    :param gene: Gene name used as primary key.
    :type gene: str
    :param alignment_obj: JSON-serializable alignment metadata dict which may contain keys like
        `seq_by_id`, `seq_by_species`, `taxid_to_species`, `species_to_taxid`,
        `alignment_extras`, `alignment_path`.
    :type alignment_obj: Dict[str, Any]
    :returns: None
    :rtype: None
    """
    cur = conn.cursor()

    seq_by_id = alignment_obj.get("seq_by_id")
    seq_by_species = alignment_obj.get("seq_by_species")
    taxid_to_species = alignment_obj.get("taxid_to_species")
    species_to_taxid = alignment_obj.get("species_to_taxid")
    alignment_extras = alignment_obj.get("alignment_extras")
    alignment_path = alignment_obj.get("alignment_path")

    seq_by_id_json = json.dumps(_sanitize_for_json(seq_by_id)) if seq_by_id is not None else None
    seq_by_species_json = json.dumps(_sanitize_for_json(seq_by_species)) if seq_by_species is not None else None
    taxid_to_species_json = json.dumps(_sanitize_for_json(taxid_to_species)) if taxid_to_species is not None else None
    species_to_taxid_json = json.dumps(_sanitize_for_json(species_to_taxid)) if species_to_taxid is not None else None
    alignment_extras_json = json.dumps(_sanitize_for_json(alignment_extras)) if alignment_extras is not None else None

    num_sequences: Optional[int] = None
    alignment_len: Optional[int] = None
    try:
        if isinstance(seq_by_id, dict) and seq_by_id:
            num_sequences = len(seq_by_id)
            first_seq = next(iter(seq_by_id.values()))
            if isinstance(first_seq, str):
                alignment_len = len(first_seq)
    except Exception:
        num_sequences = None
        alignment_len = None

    cur.execute(
        """
        INSERT OR REPLACE INTO gene_alignment
            (gene, seq_by_id_json, seq_by_species_json, taxid_to_species_json,
             species_to_taxid_json, alignment_extras_json, alignment_path,
             num_sequences, alignment_len)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            gene,
            seq_by_id_json,
            seq_by_species_json,
            taxid_to_species_json,
            species_to_taxid_json,
            alignment_extras_json,
            alignment_path,
            num_sequences,
            alignment_len,
        ),
    )


def insert_result(
    conn: sqlite3.Connection,
    gene: str,
    msa_pos: int,
    position: int,
    result_obj: Any,
) -> None:
    """Insert a single result row into the `results` table.

    The preferred input is a JSON-serializable dict produced by upstream code. As a
    fallback, objects with a `__dict__` attribute will be shallowly serialized.

    :param conn: Open SQLite connection to use for insertion.
    :type conn: sqlite3.Connection
    :param gene: Gene name.
    :type gene: str
    :param msa_pos: Zero-based integer position in the MSA.
    :type msa_pos: int
    :param position: One-based position (or -1 if unknown).
    :type position: int
    :param result_obj: JSON-serializable data structure or object with __dict__.
    :type result_obj: Any
    :returns: None
    :rtype: None
    :raises TypeError: If result_obj is not serializable and no __dict__ can be obtained.
    """
    cur = conn.cursor()

    def _extract_pair_count(r: Any) -> int:
        try:
            if isinstance(r, dict):
                if r.get("pair_details"):
                    return int(len(r.get("pair_details") or []))
                ns = r.get("node_state_details") or {}
                if isinstance(ns, dict) and ns.get("focal_states") is not None:
                    return int(len(ns.get("focal_states") or []))
                nm = r.get("node_mapping") or {}
                if isinstance(nm, dict):
                    focal_nodes = nm.get("focal_nodes")
                    if isinstance(focal_nodes, (list, tuple)):
                        return int(len(focal_nodes))
                    count = sum(1 for k in nm.keys() if str(k).startswith("focal_"))
                    return int(count) if count else 1
                return 1

            if getattr(r, "pair_details", None):
                return int(len(getattr(r, "pair_details") or [])) or 1
            if getattr(r, "node_state_details", None):
                ns = getattr(r, "node_state_details") or {}
                if isinstance(ns, dict):
                    return int(len(ns.get("focal_states") or [])) or 1
            if getattr(r, "node_mapping", None):
                nm = getattr(r, "node_mapping") or {}
                if isinstance(nm, dict):
                    focal_nodes = nm.get("focal_nodes")
                    if isinstance(focal_nodes, (list, tuple)):
                        return int(len(focal_nodes))
                    count = sum(1 for k in nm.keys() if str(k).startswith("focal_"))
                    return int(count) if count else 1
            return 1
        except Exception:
            return 1

    pair_count = _extract_pair_count(result_obj)

    if not isinstance(result_obj, dict):
        try:
            result_obj = dict(getattr(result_obj, "__dict__", {}) or {})
        except Exception as e:
            raise TypeError(
                "insert_result expects a dict or an object with a __dict__. "
                "Upstream code should pre-convert results before DB insertion."
            ) from e

    result_json = json.dumps(_sanitize_for_json(result_obj))

    cur.execute(
        "INSERT INTO results (gene, msa_pos, position, pair_count, result_json) VALUES (?, ?, ?, ?, ?)",
        (gene, int(msa_pos), int(position) if position is not None else -1, int(pair_count), result_json),
    )


def fetch_alignment_for_gene(
    conn: sqlite3.Connection, gene: str, load_posteriors: bool = False
) -> Optional[Dict[str, Any]]:
    """Fetch alignment metadata for a gene.

    :param conn: SQLite connection to query.
    :type conn: sqlite3.Connection
    :param gene: Gene name to fetch metadata for.
    :type gene: str
    :param load_posteriors: If True and `alignment_extras` includes a `posterior_dump_jsonl` path,
        attempt to load those posteriors into `alignment_extras['posterior_data']`.
    :type load_posteriors: bool
    :returns: A dict of alignment metadata or None if the row is not found.
    :rtype: Optional[Dict[str, Any]]
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT
            seq_by_id_json, seq_by_species_json, taxid_to_species_json,
            species_to_taxid_json, alignment_extras_json,
            alignment_path, num_sequences, alignment_len
        FROM gene_alignment
        WHERE gene=?
        """,
        (gene,),
    )
    row = cur.fetchone()
    if not row:
        return None

    (
        seq_by_id_json, seq_by_species_json, taxid_to_species_json,
        species_to_taxid_json, alignment_extras_json,
        alignment_path, num_sequences, alignment_len
    ) = row

    alignment_extras = json.loads(alignment_extras_json) if alignment_extras_json else None

    # Optionally load posteriors from JSONL if caller asks for it
    if load_posteriors and alignment_extras and not alignment_extras.get("posterior_data"):
        dump_path = alignment_extras.get("posterior_dump_jsonl")
        if dump_path:
            try:
                posteriors = _load_posteriors_from_jsonl(Path(dump_path))
                if posteriors:
                    alignment_extras["posterior_data"] = posteriors
            except Exception:
                pass

    return {
        "seq_by_id": json.loads(seq_by_id_json) if seq_by_id_json else None,
        "seq_by_species": json.loads(seq_by_species_json) if seq_by_species_json else None,
        "taxid_to_species": json.loads(taxid_to_species_json) if taxid_to_species_json else None,
        "species_to_taxid": json.loads(species_to_taxid_json) if species_to_taxid_json else None,
        "alignment_extras": alignment_extras,
        "alignment_path": alignment_path,
        "num_sequences": num_sequences,
        "alignment_len": alignment_len,
    }


def iter_group_keys(conn: sqlite3.Connection) -> Iterator[Tuple[str, int]]:
    """Yield (gene, msa_pos) tuples ordered for streaming grouping.

    :param conn: SQLite connection to query.
    :type conn: sqlite3.Connection
    :returns: Yields (gene, msa_pos) tuples.
    :rtype: Iterator[Tuple[str, int]]
    """
    cur = conn.cursor()
    for gene, msa_pos in cur.execute(
        "SELECT gene, msa_pos FROM results GROUP BY gene, msa_pos ORDER BY gene, msa_pos;"
    ):
        yield gene, msa_pos


def iter_results_for_group(
    conn: sqlite3.Connection, gene: str, msa_pos: int
) -> Iterator[Optional[Dict[str, Any]]]:
    """Yield result dict objects for a given (gene, msa_pos).

    Supports an optional legacy envelope that contains a pickled object like::

        {"__pickled": true, "blob": "<hex>"}

    :param conn: SQLite connection to query.
    :type conn: sqlite3.Connection
    :param gene: Gene name.
    :type gene: str
    :param msa_pos: MSA zero-based position.
    :type msa_pos: int
    :returns: Yields result dicts or None when parsing fails.
    :rtype: Iterator[Optional[Dict[str, Any]]]
    """
    cur = conn.cursor()
    for (result_json,) in cur.execute(
        "SELECT result_json FROM results WHERE gene=? AND msa_pos=? ORDER BY id",
        (gene, msa_pos),
    ):
        if not result_json:
            yield None
            continue

        try:
            obj = json.loads(result_json)
        except Exception as e:
            logger.warning(
                f"Failed to parse JSON result for gene {gene} at msa_pos {msa_pos}: {e}"
            )
            yield None
            continue

        if isinstance(obj, dict) and obj.get("__pickled") and obj.get("blob"):
            try:
                return_obj = pickle.loads(bytes.fromhex(obj["blob"]))
                if isinstance(return_obj, dict):
                    yield return_obj
                else:
                    yield dict(getattr(return_obj, "__dict__", {}) or {})
                continue
            except Exception:
                yield obj
                continue

        yield obj
