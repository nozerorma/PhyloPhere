"""Literal input formats and content identities shared by POSENRICH components."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import pandas as pd


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def read_table(path, *, quoted=True, sep='\t'):
    """Validate physical rows before pandas; preserve literal IDs and empty strings."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as stream:
        reader = csv.reader(stream, delimiter=sep, quoting=csv.QUOTE_MINIMAL if quoted else csv.QUOTE_NONE, strict=True)
        header = next(reader, None)
        if not header or len(header) != len(set(header)) or any(not c for c in header):
            raise ValueError(f'Missing/duplicate table columns: {path}')
        rows = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f'Malformed table row: {path}:{reader.line_num}')
            rows.append(row)
    return pd.DataFrame(rows, columns=header)


def identifier(value):
    if not isinstance(value, str) or not value or value.strip() != value or any(c in value for c in '\t\r\n'):
        raise ValueError(f'Invalid literal identifier: {value!r}')
    return value


def position(value):
    text = str(value)
    if not text.isascii() or not text.isdigit() or not 0 <= int(text) < 2**31:
        raise ValueError(f'Invalid zero-based position: {value!r}')
    return int(text)


def observed_table(path):
    # Postproc uses standard CSV quoting; SCORING uses literal TSV fields.
    with open(path) as stream:
        header = stream.readline().rstrip('\r\n').split('\t')
    score_col = 'asr_path_score' if 'asr_path_score' in header else 'CAAS_score'
    df = read_table(path, quoted=score_col == 'asr_path_score')
    if not {'Gene', 'Position', 'side', score_col} <= set(df.columns):
        raise ValueError(f'Missing positional columns: {path}')
    keys, none, positive = set(), set(), set()
    for row in df.to_dict('records'):
        g, p, side = identifier(row['Gene']), position(row['Position']), row['side']
        key = g, p, side
        if key in keys:
            raise ValueError(f'Duplicate position side: {key}')
        keys.add(key)
        v = float(row[score_col])
        if not math.isfinite(v) or not 0 <= v <= 1 or side not in ('top', 'bottom', 'none'):
            raise ValueError(f'Invalid side/score: {key}')
        if (side == 'none') != (v == 0) or (score_col == 'CAAS_score' and side == 'none'):
            raise ValueError(f'Incoherent side/score: {key}')
        (none if side == 'none' else positive).add((g, p))
    if none & positive:
        raise ValueError('Conflicting none and positive sides')
    df['Position'] = df['Position'].map(position).astype('int64')
    return df


def file_identity(path):
    stat=Path(path).stat()
    return {'size':stat.st_size,'mtime_ns':stat.st_mtime_ns}


def inventory(path):
    path=Path(path)
    if path.is_file():return file_identity(path)
    return {str(p.relative_to(path)):file_identity(p) for p in sorted(path.rglob('*')) if p.is_file()}


def block_inventory(blocks):
    return sorted([inventory(p) for p in blocks],key=lambda d:json.dumps(d,sort_keys=True))


def check_identity(path,expected):
    path=Path(path)
    if expected is None or not path.is_file():return False
    # New audits use normal metadata. Historical SHA values are retained as
    # provenance; an optional integrity audit can check those separately.
    if isinstance(expected,dict):return path.stat().st_size==expected.get('size')
    if __import__('os').environ.get('PHYLOPHERE_VERIFY_CHECKSUMS')=='1':return sha256_file(path)==expected
    return isinstance(expected,str) and len(expected)==64


def verify_universe(path, manifest_path=None, blocks=None):
    manifest_path = Path(manifest_path) if manifest_path else Path(path).parent/'universe_manifest.json'
    m = json.loads(manifest_path.read_text())
    if m.get('contract_version') != 2 or not check_identity(path,m.get('universe_identity',m.get('universe_sha256'))):
        raise ValueError('Unverified or changed positional universe; regenerate POSENRICH')
    if blocks is not None and (not blocks or (m.get('identity_policy')=='file_metadata' and block_inventory(blocks) != m.get('null', {}).get('blocks'))):
        raise ValueError('Null blocks absent or differ from validated universe inputs')
    return m


def finish_manifest(output_dir, manifest_name, meta, output_names=None):
    root = Path(output_dir)
    meta['source_code_sha256'] = {p.name: sha256_file(p) for p in sorted(Path(__file__).parent.glob('*.py'))}
    meta['identity_policy']='file_metadata'
    meta['outputs'] = {p.name: file_identity(p) for p in sorted(root.iterdir())
                       if p.is_file() and p.name != manifest_name and (output_names is None or p.name in output_names)}
    (root/manifest_name).write_text(json.dumps(meta, indent=2, sort_keys=True)+'\n')


def verify_registry(root):
    root = Path(root)
    m = json.loads((root/'manifest.json').read_text())
    if m.get('contract_version') != 3 or not m.get('outputs'):
        raise ValueError('Unverified registry; regenerate POSENRICH')
    for name, expected in m['outputs'].items():
        if Path(name).name != name or not check_identity(root/name,expected):
            raise ValueError(f'Changed registry input: {name}')
    return m
