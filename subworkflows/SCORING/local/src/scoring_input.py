"""Strict canonical observed rows shared by scoring and null counters."""
import csv
import math


def read_rows(path):
    out = {}
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)) or not {'Gene','Position','side','asr_path_score'} <= set(reader.fieldnames):
            raise ValueError('Missing/duplicate observed scoring columns')
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError('Malformed observed row')
            g, p, side = row['Gene'], row['Position'], row['side']
            if not g or g.strip()!=g or any(c in g for c in '\t\r\n') or not p.isdigit() or not 0<=int(p)<2**31:
                raise ValueError('Invalid gene or zero-based position')
            value = float(row['asr_path_score'])
            if not math.isfinite(value) or not 0<=value<=1 or side not in {'top','bottom','none'} or (side=='none') != (value==0):
                raise ValueError('Invalid side/asr_path_score')
            key = (g,int(p),side)
            if key in out:
                raise ValueError(f'duplicate row for {key}')
            out[key] = value
    none = {(g,p) for g,p,s in out if s=='none'}
    if any((g,p) in none for g,p,s in out if s!='none'):
        raise ValueError('Conflicting none and positive sides')
    return out
