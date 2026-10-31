#!/usr/bin/env python3
"""Select exactly N validated null cycles and construct the observed b_0 block."""
import argparse
import csv
import json
import hashlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.design_io import TAG, integer, read_design, read_labelings, pair_metadata

PAIR_HEADER = ['cycle', 'hypothesis_id', 'pair', 'species1', 'species2', 'pss_score']


def labelings_of_design(design, fop):
    return {('b_0~' + h if fop else 'b_0'):
            (tuple(s['1'] for s in pairs.values()), tuple(s['0'] for s in pairs.values()))
            for h, pairs in design.items()}


def write_block(directory, labelings, weights=None):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'resample_perms.tab').open('w') as out:
        for tag, (fg, bg) in labelings.items():
            out.write(f"{tag}\t{','.join(fg)}\t{','.join(bg)}\n")
    if weights is not None:
        with (directory / 'fop_pairs.tsv').open('w', newline='') as out:
            writer = csv.writer(out, delimiter='\t', lineterminator='\n')
            writer.writerow(PAIR_HEADER)
            for tag, (fg, bg) in labelings.items():
                m = TAG.fullmatch(tag)
                for i, (f, b) in enumerate(zip(fg, bg), 1):
                    key = (f'b_{m[1]}', m[2] or 'H1', i)
                    writer.writerow([*key, f, b, weights[key]])


def prepare(source, config, n, block_size, multi, allow_relaxed=False):
    n = integer(n, 'caas_full_perms')
    block_size = integer(block_size, 'caas_cycle_block_size', 1)
    design = read_design(config, multi)
    active = set(design)
    k = len(design['H1'])
    source = Path(source)
    if source.is_dir():
        fop_files = list(source.rglob('fop_labelings.tab'))
        if len(fop_files) > 1:
            raise ValueError('Ambiguous harvest: multiple fop_labelings.tab files')
        files = fop_files or sorted(source.rglob('resample_*.tab'))
        metadata = fop_files[0].parent / 'fop_pairs.tsv' if fop_files else None
    elif source.is_file():
        files = [source]
        metadata = source.parent / 'fop_pairs.tsv'
    else:
        raise ValueError(f'Resample input not found: {source}')
    raw = read_labelings(files)
    source_fop = any('~' in tag for tag in raw)
    use_fop = source_fop or len(active) > 1
    grouped = {}
    for tag, labeling in raw.items():
        m = TAG.fullmatch(tag)
        cycle = int(m[1]); hyp = m[2] or 'H1'
        if cycle == 0 or hyp not in active:
            continue
        if len(labeling[0]) != k:
            raise ValueError(f'Pair count differs from observed design: {tag}')
        if hyp in grouped.setdefault(cycle, {}):
            raise ValueError(f'Duplicate base-cycle hypothesis: {tag}')
        grouped[cycle][hyp] = labeling
    if len(grouped) < n:
        raise ValueError(f'Requested {n} distinct null cycles; only {len(grouped)} available')
    selected, reasons = {}, []
    for cycle in sorted(grouped)[:n]:
        hypotheses = grouped[cycle]
        if set(hypotheses) != active:
            if not allow_relaxed:
                raise ValueError(f'Cycle b_{cycle} does not match active hypotheses {sorted(active)}')
            reasons.append(f'b_{cycle}: incomplete hypothesis harvest')
        if 'H1' not in hypotheses:
            raise ValueError(f'Canonical H1 missing from cycle b_{cycle}')
        for hyp in sorted(hypotheses, key=lambda h: int(h[1:])):
            tag = f'b_{cycle}~{hyp}' if use_fop else f'b_{cycle}'
            selected[tag] = hypotheses[hyp]
    provenance_dir = source if source.is_dir() else source.parent
    provenances = list(provenance_dir.rglob('permulation_provenance.json'))
    if len(provenances) > 1:
        raise ValueError('Ambiguous harvest provenance')
    harvest_provenance = json.loads(provenances[0].read_text()) if provenances else None
    if n and harvest_provenance:
        if harvest_provenance.get('relaxed_reasons'):
            if not allow_relaxed:
                raise ValueError('Reused harvest contains relaxed null draws; enable ct_allow_relaxed_null explicitly')
            reasons.extend(harvest_provenance['relaxed_reasons'])
        canonical = ''.join(f"{sides[group]}\t{group}\t{pair}\n"
                            for pair, sides in design['H1'].items() for group in ('0', '1'))
        expected_signature = hashlib.sha256(canonical.encode()).hexdigest()
        stored = harvest_provenance.get('observed_canonical_signature')
        if stored is not None and stored != expected_signature:
            raise ValueError('Reused harvest belongs to a different observed canonical design')
    weights = None
    if use_fop:
        if selected and metadata is not None and metadata.is_file():
            weights = pair_metadata(metadata, selected)
        elif selected and (source_fop or len(active) > 1):
            if not allow_relaxed:
                raise ValueError('FOP reuse requires matching fop_pairs.tsv weights')
            reasons.append('missing null pair weights; equal weights used')
        if weights is None:
            weights = {(f'b_{TAG.fullmatch(t)[1]}', TAG.fullmatch(t)[2] or 'H1', i): 1.0
                       for t, (fg, _) in selected.items() for i in range(1, len(fg)+1)}
    observed = labelings_of_design(design, use_fop)
    observed_weights = None
    if use_fop:
        p = Path(config)
        metadata = (p if p.is_dir() else p.parent) / 'contrast_hypotheses_pairs.tsv'
        if metadata.is_file():
            observed_weights = pair_metadata(metadata, observed, observed=True)
        else:
            if len(active) > 1 and not allow_relaxed:
                raise ValueError('Multiple observed hypotheses require their PSS pair metadata')
            if len(active) > 1:
                reasons.append('missing observed pair weights; equal weights used')
            observed_weights = {('b_0', h, pair): 1.0 for h, pairs in design.items() for pair in pairs}
    write_block(Path('blocks/b0'), observed, observed_weights)
    cycles = sorted(grouped)[:n]
    for start in range(0, n, block_size):
        keep = set(cycles[start:start+block_size])
        labels = {tag: lab for tag, lab in selected.items() if int(TAG.fullmatch(tag)[1]) in keep}
        write_block(Path(f'blocks/blk_{start//block_size+1:04d}'), labels, weights)
    with Path('resample_perms.tab').open('w') as out:
        for tag, (fg, bg) in {**observed, **selected}.items():
            out.write(f"{tag}\t{','.join(fg)}\t{','.join(bg)}\n")
    Path('ct_preparation.json').write_text(json.dumps({
        'schema_version': 1, 'requested_cycles': n, 'selected_cycles': len(cycles),
        'active_hypotheses': sorted(active), 'allow_relaxed_null': allow_relaxed,
        'relaxed_reasons': reasons, 'harvest_provenance': harvest_provenance,
        'harvest_provenance_available': harvest_provenance is not None, 'conservation_semantics': 'within_pair_ids_v2',
    }, indent=2) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--n', required=True)
    ap.add_argument('--block-size', required=True)
    ap.add_argument('--multi', choices=['true', 'false'], required=True)
    ap.add_argument('--allow-relaxed', choices=['true', 'false'], default='false')
    a = ap.parse_args()
    prepare(a.source, a.config, a.n, a.block_size, a.multi == 'true', a.allow_relaxed == 'true')
