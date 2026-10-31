"""Strict paired-design and labeling readers shared by CT preparation and replay."""
import csv
import math
from pathlib import Path
import re

TAG = re.compile(r'b_(0|[1-9][0-9]*)(?:~(H[1-9][0-9]*))?\Z')
HYP_FILE = re.compile(r'traitfile_(H[1-9][0-9]*)\.tab\Z')


def integer(value, name, minimum=0):
    if not re.fullmatch(r'[0-9]+', str(value)) or int(value) < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}: {value!r}')
    return int(value)


def trait_pairs(path):
    pairs, seen = {}, set()
    with Path(path).open(newline='', encoding='utf-8-sig') as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if not row:
                continue
            if len(row) != 3 or row[1] not in {'0', '1'} or not row[0].strip():
                raise ValueError(f'Invalid species/group/pair row in {path}: {row}')
            sp = row[0].strip()
            if ',' in sp or any(ch.isspace() for ch in sp):
                raise ValueError(f'Species names must have no commas or whitespace: {sp!r}')
            pair = integer(row[2], 'pair', 1)
            if sp in seen or row[1] in pairs.setdefault(pair, {}):
                raise ValueError(f'Duplicate species or pair side in {path}: {row}')
            seen.add(sp)
            pairs[pair][row[1]] = sp
    if not pairs or any(set(sides) != {'0', '1'} for sides in pairs.values()):
        raise ValueError(f'Empty or incomplete paired design: {path}')
    if set(pairs) != set(range(1, len(pairs) + 1)):
        raise ValueError(f'Pair IDs must be contiguous from 1 in {path}')
    return dict(sorted(pairs.items()))


def design_files(path, multi=True):
    p = Path(path)
    if p.is_file():
        return {'H1': p}
    if not p.is_dir():
        raise ValueError(f'Design not found: {p}')
    files = {}
    for f in sorted(p.glob('traitfile_H*.tab')):
        match = HYP_FILE.fullmatch(f.name)
        if match is None:
            raise ValueError(f'Invalid hypothesis filename: {f.name}')
        if multi or match[1] == 'H1':
            files[match[1]] = f
    if not files and (p / 'traitfile.tab').is_file():
        files = {'H1': p / 'traitfile.tab'}
    if 'H1' not in files:
        raise ValueError(f'Canonical H1 is missing from {p}')
    if multi and {int(h[1:]) for h in files} != set(range(1, len(files)+1)):
        raise ValueError('Hypothesis IDs must be contiguous from H1')
    return files


def read_design(path, multi=True):
    result = {h: trait_pairs(f) for h, f in design_files(path, multi).items()}
    if len({len(pairs) for pairs in result.values()}) != 1:
        raise ValueError('Every hypothesis must have the canonical pair count')
    return result


def read_labelings(files):
    result = {}
    for path in files:
        with Path(path).open(newline='') as fh:
            for row in csv.reader(fh, delimiter='\t'):
                if not row:
                    continue
                if len(row) != 3 or not TAG.fullmatch(row[0]):
                    raise ValueError(f'Invalid labeling row in {path}: {row}')
                tag, fg, bg = row[0], tuple(row[1].split(',')), tuple(row[2].split(','))
                if (not fg or len(fg) != len(bg) or any(not s or s != s.strip() for s in fg + bg)
                        or len(set(fg + bg)) != len(fg + bg)):
                    raise ValueError(f'Invalid paired species lists for {tag} in {path}')
                if tag in result:
                    raise ValueError(f'Duplicate labeling tag: {tag}')
                result[tag] = (fg, bg)
    return result


def pair_metadata(path, labelings, *, observed=False):
    """Validate exact coverage, species orientation and finite nonnegative weights."""
    expected = {}
    for tag, (fg, bg) in labelings.items():
        m = TAG.fullmatch(tag)
        for i, pair in enumerate(zip(fg, bg), 1):
            expected[(f'b_{m[1]}', m[2] or 'H1', i)] = pair
    out = {}
    with Path(path).open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = {'hypothesis_id', 'pair', 'species1', 'species2', 'pss_score'}
        if not observed:
            fields.add('cycle')
        if not fields <= set(reader.fieldnames or []):
            raise ValueError(f'Missing pair metadata columns in {path}')
        for row in reader:
            key = ('b_0' if observed else row['cycle'], row['hypothesis_id'], integer(row['pair'], 'metadata pair', 1))
            if key not in expected:
                continue  # metadata for inactive hypotheses or unselected cycles
            if key in out or (row['species1'], row['species2']) != expected[key]:
                raise ValueError(f'Duplicate or inconsistent pair metadata: {key}')
            weight = float(row['pss_score'])
            if not math.isfinite(weight) or weight < 0:
                raise ValueError(f'Invalid PSS weight: {key}')
            out[key] = weight
    if set(out) != set(expected):
        raise ValueError(f'Pair metadata does not cover active labelings: {path}')
    return out
