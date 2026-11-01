#!/usr/bin/env python3
"""Prepare a declared hybrid analysis: def2 scores/null, current CAAS positional universe.

Never runs ASR/PAML. Missing universe shards are computed by CT's universe-only export.
This does not relabel legacy scientific provenance as v4 or regenerate SCORING p-values.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def canonical(name, genes):
    names = [name, re.sub(r'\.[A-Z][a-z]+_[a-z]+$', '', name)]
    names += [re.sub(r'\.\d+$', '', n) for n in names[:]]
    return next((n for n in names if n in genes), None)


def read_universes(paths, genes):
    found = {}
    for path in paths:
        for line in path.read_text().splitlines():
            raw, text = line.split('\t')
            gene = canonical(raw, genes)
            if gene is None:
                continue
            values = tuple(sorted(map(int, text.split(',')))) if text not in ('', 'NULL') else ()
            if len(set(values)) != len(values) or any(v < 0 for v in values):
                raise ValueError(f'Invalid universe: {path}')
            if gene in found and found[gene] != values:
                raise ValueError(f'Conflicting universe shards for {gene}')
            found[gene] = values
    return found


def complete_universe(repo, work, params, genes, found, output, workers):
    missing = sorted(genes - found.keys())
    if not missing:
        return found
    preps = list(work.glob('*/*/ct_preparation.json'))
    if len(preps) != 1:
        raise ValueError('Expected one new-run preparation; cannot select an observed design safely')
    prep = preps[0].parent
    design = prep / 'traitfiles_ok_dir'
    b0 = prep / 'blocks/b0/resample_perms.tab'
    if not design.is_dir() or not b0.is_file():
        raise FileNotFoundError('New-run design / b0 not found')
    alignment_dir = Path(params['alignment'])
    alignments = {}
    for f in sorted(alignment_dir.iterdir()):
        if f.suffix.lower() not in ('.fa', '.fasta', '.faa', '.fas'):
            continue
        gene = canonical(f.stem, genes)
        if gene:
            if gene in alignments:
                raise ValueError(f'Ambiguous alignment for {gene}')
            alignments[gene] = f
    absent = set(missing) - alignments.keys()
    if absent:
        raise ValueError(f'No alignment for {len(absent)} background genes: {sorted(absent)[:10]}')
    k = len({r.split('\t')[2] for r in (design/'traitfile_H1.tab').read_text().splitlines() if r})
    flags = ['--patterns', str(params['patterns']), '--max_conserved',
             str(int(k*(1-float(params['min_divergent_fraction'])))),
             '--min_non_overlapping', str(params['min_non_overlapping'])]
    if params['miss_pair'] in (True, 'true'):
        flags.append('--miss_pair')
    for flag, key in [('max_fg_gaps','max_fg_gaps_fraction'),('max_bg_gaps','max_bg_gaps_fraction'),
                      ('max_gaps','max_gaps_fraction'),('max_fg_miss','max_fg_miss_fraction'),
                      ('max_bg_miss','max_bg_miss_fraction'),('max_miss','max_miss_fraction')]:
        flags += ['--'+flag, str(int(k*float(params[key])))]
    chunks = [missing[i:i+100] for i in range(0, len(missing), 100)]
    def run(item):
        i, batch = item
        directory = output/'universe_completion'/f'batch_{i:04d}'
        directory.mkdir(parents=True, exist_ok=True)
        manifest = directory/'genes.tsv'
        manifest.write_text(''.join(f'{g}\t{alignments[g]}\n' for g in batch))
        cmd = [sys.executable, str(repo/'subworkflows/CT/local/ct'), 'perm-replay',
               '--manifest', str(manifest), '--alignment_dir', str(alignment_dir),
               '-t', str(design), '-s', str(b0), '--fmt', params['ali_format'],
               '--export_b0_universe', str(directory/'{id}.b0.universe'), *flags]
        with (directory/'ct.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        return list(directory.glob('*.b0.universe'))
    print(f'Universe-only CT: completing {len(missing)} genes; ASR/PAML disabled by construction', flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for paths in pool.map(run, enumerate(chunks)):
            found.update(read_universes(paths, genes))
    if genes - found.keys():
        raise ValueError('Incomplete positional-universe coverage')
    return found


def prepare(old, output, found, genes):
    registry = output/'registry'; registry.mkdir(parents=True, exist_ok=True)
    pairs = [(g, pos) for g in sorted(found) for pos in found[g]]
    universe = pd.DataFrame(pairs, columns=['gene','position'])
    if universe.empty:
        raise ValueError('Empty updated universe')
    ukeys = pd.MultiIndex.from_frame(universe)
    null = pq.read_table(old/'caas_permulation/null_positions.parquet').to_pandas()
    if null.duplicated(['gene','position']).any():
        raise ValueError('Duplicate old null positions')
    null_index = pd.MultiIndex.from_frame(null[['gene','position']])
    universe['nd'] = pd.Series(null.nd.to_numpy(), index=null_index).reindex(ukeys, fill_value=0).to_numpy(dtype=np.int64)
    scores = pd.read_csv(old/'scoring/position_scores.tsv', sep='\t', dtype={'Gene':str})
    if scores.duplicated(['Gene','Position','side']).any():
        raise ValueError('Duplicate old observed position side')
    scores['CAAS_score'] = pd.to_numeric(scores.CAAS_score, errors='raise')
    if not np.isfinite(scores.CAAS_score).all() or not scores.CAAS_score.between(0,1).all():
        raise ValueError('Invalid old observed score')
    inside = pd.MultiIndex.from_frame(scores[['Gene','Position']]).isin(ukeys)
    excluded = scores.loc[~inside].copy()
    excluded.to_csv(output/'observed_outside_updated_universe.tsv', sep='\t', index=False)
    # Retain scientific scores, never advertise old p/q as recalibrated for the new U.
    scores = scores.loc[inside, ['Gene','Position','side','CAAS_score']].copy()
    scores.to_csv(output/'position_scores.tsv', sep='\t', index=False)
    for side in ('top','bottom'):
        sub=scores[(scores.side==side)&(scores.CAAS_score>0)]
        universe['observed_'+side]=ukeys.isin(pd.MultiIndex.from_frame(sub[['Gene','Position']])).astype(np.int8)
    universe['observed_positive']=((universe.observed_top+universe.observed_bottom)>0).astype(np.int8)
    ud = output/'universe';ud.mkdir(exist_ok=True)
    universe.to_csv(ud/'universe.tsv.gz',sep='\t',index=False)
    digest=hashlib.sha256('\n'.join(sorted(genes)).encode()).hexdigest()
    manifest={'contract_version':1,'analysis_type':'legacy_ASR_scores_with_updated_CAAS_universe',
              'gene_set_sha256':digest,'metrics':{'n_genes':len(genes),'n_universe':len(universe),
              'n_observed_positive':int(universe.observed_positive.sum())}}
    (ud/'universe_manifest.json').write_text(json.dumps(manifest,indent=2))
    original=old/'posenrich/registry'
    sites=pd.read_csv(original/'sites.tsv',sep='\t',dtype=str,keep_default_na=False)
    sites=sites[(sites.collection_id=='custom')&sites.term_id.str.startswith('PP2A:')]
    if sites.empty:raise ValueError('No mapped PP2A sites in def2')
    sites.to_csv(registry/'sites.tsv',sep='\t',index=False)
    collections=pd.read_csv(original/'collections.tsv',sep='\t',dtype=str)
    collections=collections[(collections.collection_id=='custom')&collections.term_id.str.startswith('PP2A:')]
    collections.to_csv(registry/'collections.tsv',sep='\t',index=False)
    first=True
    for positions in pd.read_csv(original/'positions.tsv',sep='\t',dtype={'Gene':str},chunksize=100000):
        positions=positions[pd.MultiIndex.from_frame(positions[['Gene','Position']]).isin(ukeys)]
        positions.to_csv(registry/'positions.tsv',sep='\t',index=False,header=first,mode='w' if first else 'a')
        first=False
    scores.to_csv(registry/'position_evidence.tsv',sep='\t',index=False)
    pd.DataFrame(columns=['collection_id','term_id','term_label','Gene','start','end',
                          'source_record_id','mapping_status','detail']).to_csv(registry/'regions.tsv',sep='\t',index=False)
    rm=json.loads((original/'manifest.json').read_text())
    rm['recovery']={'analysis_type':manifest['analysis_type'],'only_collection':'custom:PP2A',
                    'old_position_pvalues_reused':False}
    (registry/'manifest.json').write_text(json.dumps(rm,indent=2))
    audit={'analysis_type':manifest['analysis_type'],'old_root':str(old),'genes':len(genes),
           'universe_positions':len(universe),'observed_rows_kept':len(scores),
           'observed_rows_outside_updated_universe':len(excluded),'pp2a_terms':len(collections),
           'new_universe_positions_absent_from_old_null':int((universe.nd==0).sum()),
           'limits':['Preserves def2 ASR, positive detections, scores and gene/null filters.',
                     'Updates only positional opportunity within def2 cleaned background genes.',
                     'Does not claim v4 scoring or recomputed position p-values/FDR.',
                     'New null support is old scored-position support restricted to U; zero-score rows are implicit.']}
    (output/'recovery_audit.json').write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','old-root','new-work','new-params','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=False)
    genes=set((a.old_root/'postproc/cleaned_backgrounds/cleaned_background_main.txt').read_text().splitlines())-{'Gene','gene',''}
    paths=[]
    for d in a.new_work.glob('*/*/universe'):
        exitfile=d.parent/'.exitcode'
        if exitfile.is_file() and exitfile.read_text().strip()=='0':paths.extend(d.glob('*.b0.universe'))
    found=read_universes(paths,genes)
    reused=len(found)
    found=complete_universe(a.repo.resolve(),a.new_work.resolve(),json.loads(a.new_params.read_text()),genes,found,a.output,a.workers)
    prepare(a.old_root.resolve(),a.output,found,genes)
    audit=json.loads((a.output/'recovery_audit.json').read_text());audit.update(universe_genes_reused=reused,universe_genes_computed=len(genes)-reused)
    (a.output/'recovery_audit.json').write_text(json.dumps(audit,indent=2))


if __name__=='__main__':main()
