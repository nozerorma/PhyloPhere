"""Validated, disk-backed discovery grouping; memory is bounded by one gene."""
import csv
import hashlib
import json
import re
import sqlite3
import tempfile
from collections.abc import Mapping
from pathlib import Path

from src.utils.amino import normalize_amino_list
from src.core.labelings import hyp_id


def position_digest(positions):
    return hashlib.sha256(''.join(f'{p}\n' for p in sorted(set(positions))).encode()).hexdigest()


def validate_row(row, design=None):
    gene = row.get('gene') or ''
    if not gene or gene != gene.strip() or any(c in gene for c in '/\\\t\r\n'):
        raise ValueError(f'Invalid gene identifier: {gene!r}')
    if not re.fullmatch(r'[0-9]+', row.get('position') or ''):
        raise ValueError(f'Invalid discovery position: {row.get("position")}')
    if row.get('pattern') not in {'1','2','3','4'}:
        raise ValueError(f'Invalid discovery pattern: {row.get("pattern")}')
    parts = (row.get('caas') or '').split('/')
    if len(parts) != 2 or not all(normalize_amino_list(part) for part in parts):
        raise ValueError(f'Invalid CAAS residues: {row.get("caas")}')
    overlap = len(set(normalize_amino_list(parts[0])) & set(normalize_amino_list(parts[1])))
    supplied = row.get('shared_residue_count')
    if supplied not in (None, '') and (not supplied.isdigit() or int(supplied) != overlap):
        raise ValueError('shared_residue_count disagrees with the residues')
    supplied = row.get('is_conserved_meta')
    if supplied not in (None, ''):
        if supplied.lower() not in {'true','false','1','0'} or (supplied.lower() in {'true','1'}) != bool(overlap):
            raise ValueError('is_conserved_meta disagrees with residue-type overlap')
    conserved = row.get('conserved_pair') or ''
    if conserved:
        prefix, sep, ids = conserved.partition(':')
        values = (ids if sep else conserved).split(',') if (ids if sep else conserved) else []
        if any(not re.fullmatch('[1-9][0-9]*', value) for value in values) or len(set(values)) != len(values):
            raise ValueError(f'Invalid conserved pair IDs: {conserved}')
        if sep and (not prefix.isdigit() or int(prefix) != len(values)):
            raise ValueError('Conserved-pair count does not match its IDs')
    raw = row.get('trait') or ''
    if design is not None:
        match = re.search(r'H([0-9]+)', raw)
        if (match and int(match[1]) not in design) or (not match and len(design) != 1):
            raise ValueError(f'Discovery hypothesis is absent from the design: {raw!r}')
    return gene, int(row['position']), hyp_id(raw)


class DiscoveryRows(Mapping):
    def __init__(self, source, scratch_dir=None, design=None):
        self.tmp = tempfile.TemporaryDirectory(prefix='.discovery-', dir=scratch_dir or Path.cwd())
        self.db = sqlite3.connect(str(Path(self.tmp.name) / 'rows.sqlite'), check_same_thread=False)
        self.db.execute('CREATE TABLE rows (gene TEXT, position INTEGER, hyp TEXT, row TEXT, UNIQUE(gene,position,hyp))')
        source = Path(source)
        self.paths = [source] if source.is_file() else sorted(source.glob('*.b0.discovery.tsv'))
        if not source.exists():
            raise FileNotFoundError(source)
        for path in self.paths:
            with path.open(newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                fields = reader.fieldnames or []
                if len(set(fields)) != len(fields) or not {'gene','position','caas','pattern'} <= set(fields):
                    raise ValueError(f'Invalid discovery header: {path}')
                for row in reader:
                    if None in row or any(v is None for v in row.values()):
                        raise ValueError(f'Wrong discovery column count: {path}')
                    gene, position, hyp = validate_row(row, design)
                    try:
                        self.db.execute('INSERT INTO rows VALUES (?,?,?,?)', (gene, position, hyp, json.dumps(row)))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(f'Duplicate discovery gene/position/hypothesis: {gene}, {position}, {hyp}') from exc
        self.db.execute('CREATE INDEX genes ON rows(gene)')
        self.db.commit()

    def __iter__(self):
        return (row[0] for row in self.db.execute('SELECT DISTINCT gene FROM rows ORDER BY gene'))

    def __len__(self):
        return self.db.execute('SELECT count(DISTINCT gene) FROM rows').fetchone()[0]

    def __getitem__(self, gene):
        rows = [json.loads(row[0]) for row in self.db.execute('SELECT row FROM rows WHERE gene=? ORDER BY rowid', (gene,))]
        if not rows:
            raise KeyError(gene)
        return rows

    def close(self):
        self.db.close()
        self.tmp.cleanup()
