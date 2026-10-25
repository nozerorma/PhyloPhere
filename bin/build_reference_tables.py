#!/usr/bin/env python3
"""Build explicit positional references without changing the source annotations.

project-ppga projects the verified PPGA productive UCR/Pfam coordinates and
validates FUBAR sites against the original MAP; legacy tables are never guessed.
marine-ucr recomputes sliding UCR on the supplied CAAS alignments with the existing
Valdar/sliding implementation and defaults. No taxonomic clade statistics are used.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def write_manifest(output, sources, counts, **extra):
    data=dict(status='validated',coordinate_space='CAAS_0based',counts=dict(counts),
              sources={str(p):sha(p) for p in sources},table_sha256=sha(output),**extra)
    Path(str(output)+'.manifest.json').write_text(json.dumps(data,indent=2)+'\n')


@lru_cache(maxsize=8)
def map_gene(map_dir,gene):
    paths=list(Path(map_dir).glob(gene+'.*.map.tsv'))
    if not paths: paths=list(Path(map_dir).glob(gene+'.map.tsv'))
    if len(paths)!=1:raise ValueError(f'{gene}: expected exactly one MAP, found {len(paths)}')
    raw={};by_site={}
    for r in csv.DictReader(paths[0].open(),delimiter='\t'):
        if r['status']!='selected':continue
        col=int(r['prot_ali_col'])-1
        if col in by_site:raise ValueError(f'{gene}: duplicate MAP column')
        by_site[col]=r;raw[int(r['ori_codon_col'])]=col
    return raw,by_site


def project_ppga(base, output):
    from ppga_position_references import project_ppga as project_verified_ppga
    return project_verified_ppga(base, output)


def ucr_gene(path):
    import compute_variability as cv
    from detect_ucr import detect_sliding
    seqs=cv.read_fasta(path)
    widths={len(s) for _,s in seqs}
    if len(widths)!=1:raise ValueError(f'Unequal alignment lengths: {path}')
    weights=cv.henikoff_weights([s for _,s in seqs])
    per_col,_=cv.compute_per_column(seqs,weights)
    # The regular producer serializes these values with six decimal places.
    rows=[dict(position=i+1,C_trident=round(r[3],6),variability=round(r[4],6),g=round(r[2],6))
          for i,r in enumerate(per_col)]
    gene=Path(path).name.split('.')[0];records=[]
    for u in detect_sliding(gene,rows,10,.85,.3,5):
        for pos in range(u['flank_start'],u['flank_end']+1):
            region='flank_up' if pos<u['start_pos'] else 'flank_down' if pos>u['end_pos'] else 'core'
            r=rows[pos-1]
            records.append(dict(gene=gene,ucr_id=u['ucr_id'],method='sliding',position=pos-1,
                                region_type=region,C_trident=r['C_trident'],variability=r['variability'],
                                g=r['g'],coordinate_space='CAAS_0based'))
    return records,str(path),sha(path)


def marine_ucr(alignment_dir,output,chunk,chunks,cpus):
    output.parent.mkdir(parents=True,exist_ok=True)
    files=sorted(p for p in alignment_dir.iterdir() if p.suffix.lower() in ('.fa','.fasta','.faa','.fas'))[chunk::chunks]
    fields=['gene','ucr_id','method','position','region_type','C_trident','variability','g','coordinate_space']
    fingerprints={};counts=Counter()
    with gzip.open(output,'wt',newline='') as f,ProcessPoolExecutor(max_workers=cpus) as pool:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
        for records,path,digest in pool.map(ucr_gene,files,chunksize=8):
            writer.writerows(records);fingerprints[Path(path).name]=digest
            counts['genes_scanned']+=1;counts['positions']=counts['positions']+len(records)
    write_manifest(output,[],counts,alignment_dir=str(alignment_dir.resolve()),alignments=fingerprints,
                   method='sliding',window_size=10,window_threshold=.85,window_max_pos_gap=.3,flank_size=5,
                   valdar_exponents=[1.,1.,1.])
    print(dict(counts),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['project-ppga','marine-ucr','alignment-ucr'])
    p.add_argument('--base',type=Path);p.add_argument('--alignment-dir',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--chunk',type=int,default=0)
    p.add_argument('--chunks',type=int,default=1);p.add_argument('--cpus',type=int,default=1);a=p.parse_args()
    if a.mode=='project-ppga':project_ppga(a.base,a.output)
    else:marine_ucr(a.alignment_dir,a.output,a.chunk,a.chunks,a.cpus)
