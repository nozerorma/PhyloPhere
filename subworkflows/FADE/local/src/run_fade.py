#!/usr/bin/env python3
"""Shared single/batch FADE runner: exact annotation, audited exclusions, fail closed."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]/'SELECTION/local/src'))
from annotate_tree_fg import annotate
from alignment_io import read_alignment
from Bio import SeqIO
from prepare_alignments import validate_gene
from fade_json_to_tables import reduce_json, write_tables


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as fh:
        for block in iter(lambda:fh.read(1048576),b''):h.update(block)
    return h.hexdigest()


def run_group(jobs, settings, cpu, container):
    genes,sites,audits=[],[],[]
    gene=jobs[0][0]
    cache=Path(f'{gene}.baseline.cache').absolute()
    cache.unlink(missing_ok=True)
    shared = (len(jobs)>1 and settings['background_scope']=='all'
              and len({(digest(row[2]),digest(row[3])) for row in jobs})==1)
    for _,direction,fasta,tree,fg,bg in jobs:
        stem=f'{gene}.{direction}'
        log=Path(stem+'.hyphy.log')
        output=Path(stem+'.FADE.json')
        audit=dict(gene=gene,direction=direction,status='technical_failure',
                   code_sha256=Path('inputs.sha256').read_text().strip().split(':')[0] if Path('inputs.sha256').exists() else None,
                   input_sha256={k:digest(p) for k,p in [('alignment',fasta),('tree',tree),('foreground',fg)]+([('background',bg)] if bg else [])})
        try:
            annotated=Path(stem+'_fg.nwk')
            filtered=Path(stem+'.fa')
            result=annotate(tree,fasta,fg,annotated,bg_path=bg,
                            fasta_out=filtered if bg else None, internal_mode=settings['internal_nodes'])
            audit.update(result)
            if result['status']=='eligible':
                # SELECTION_PREP already retained tree taxa in all mode. Reuse
                # that FASTA; only explicit extra-taxon filtering needs a copy.
                if bg is None:
                    if result['excluded_taxa']:
                        retained=set(result['retained_taxa'])
                        SeqIO.write([r for r in read_alignment(fasta) if r.id in retained],filtered,'fasta')
                    else:
                        filtered=Path(fasta)
                cmd=(['/usr/local/bin/_entrypoint.sh','hyphy'] if container else ['hyphy'])
                cmd += [str(HERE/'FADE_phylophere.bf'),'--alignment',str(filtered),
                        '--tree',str(annotated),'--branches','Foreground','--model',settings['model'],
                        '--method',settings['method'],'--grid',str(settings['grid']),
                        '--concentration_parameter',str(settings['concentration']),
                        '--cpu',str(cpu),'--cache',str(cache) if shared else 'NONE','--output',str(output.absolute())]
                if settings['method']!='Variational-Bayes':
                    for flag,key in [('chains','chains'),('chain-length','chain_length'),('burn-in','burn_in'),('samples','samples')]:
                        cmd += ['--'+flag,str(settings[key])]
                env=os.environ.copy()
                for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS']:env[name]=str(cpu)
                with log.open('w') as fh:
                    proc=subprocess.run(cmd,stdout=fh,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env)
                if proc.returncode:
                    raise RuntimeError(f'HyPhy exited {proc.returncode}; see {log}')
                js=json.loads(output.read_text())
                row,hit_sites=reduce_json(js,gene,direction,settings['bf_threshold'],result['alignment_columns'])
                genes.append(row);sites.extend(hit_sites)
                audit.update(status='complete',hyphy_version=js.get('runtime'),
                             baseline_reused='Loaded baseline model fit from cache' in log.read_text(),
                             n_sites_tested=row[2],n_significant_sites=len(hit_sites))
                output.unlink()
        except Exception as exc:
            audit.update(status='technical_failure',error=str(exc))
        Path(stem+'.fade.audit.json').write_text(json.dumps(audit,indent=2)+'\n')
        audits.append(audit)
    cache.unlink(missing_ok=True)
    return genes,sites,audits


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',required=True)
    parser.add_argument('--settings',required=True)
    parser.add_argument('--batch-id',required=True)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--cpu-per-worker',type=int,default=1)
    parser.add_argument('--container',action='store_true')
    args=parser.parse_args()
    if args.workers<1 or args.cpu_per_worker<1:raise ValueError('Invalid worker/CPU count')
    settings=json.loads(Path(args.settings).read_text())
    groups,seen={},set()
    with open(args.manifest,newline='') as fh:
        for row in csv.reader(fh,delimiter='\t',quoting=csv.QUOTE_NONE):
            if len(row)!=6:raise ValueError(f'Expected six manifest columns: {row}')
            gene,direction,*paths=row
            validate_gene(gene)
            if direction not in {'top','bottom'} or (gene,direction) in seen:raise ValueError('Duplicate gene/direction or invalid direction')
            seen.add((gene,direction))
            fasta,tree,fg,bg=paths
            for p in [fasta,tree,fg]+([bg] if bg!='NO_FILE' else []):
                if not Path(p).is_file():raise ValueError(f'Input not found: {p}')
            groups.setdefault(gene,[]).append((gene,direction,fasta,tree,fg,None if bg=='NO_FILE' else bg))
    if not groups:raise ValueError('Empty FADE manifest')
    jobs=[sorted(g,key=lambda r:['top','bottom'].index(r[1])) for g in groups.values()]
    with ThreadPoolExecutor(max_workers=min(args.workers,len(jobs))) as pool:
        results=list(pool.map(lambda g:run_group(g,settings,args.cpu_per_worker,args.container),jobs))
    audits=[a for _,_,aa in results for a in aa]
    failures=[a for a in audits if a['status']=='technical_failure']
    if failures:
        raise RuntimeError('FADE technical failures: '+'; '.join(f"{a['gene']}/{a['direction']}: {a['error']}" for a in failures))
    if len(audits)!=len(seen):raise RuntimeError('Incomplete FADE audit coverage')
    write_tables(args.batch_id+'.fade_genes.tsv',args.batch_id+'.fade_sites.tsv',
                 [r for gg,_,_ in results for r in gg],[r for _,ss,_ in results for r in ss])


if __name__=='__main__':
    main()
