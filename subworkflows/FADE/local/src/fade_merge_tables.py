#!/usr/bin/env python3
"""Validate full gene/direction coverage and merge FADE results for one direction."""
import argparse
import csv
import json
import hashlib
import math
from pathlib import Path
from fade_json_to_tables import AA_LETTERS, GENES_COLUMNS, SITES_COLUMNS, number

SUMMARY=['gene','n_sites_tested','n_sig_site_aa','max_bf','top_target_aa','top_site']
SITES=['gene','position','max_bf','target_aa']


def integer(value, minimum=0):
    if not str(value).isdigit() or int(value)<minimum:raise ValueError(f'Invalid integer: {value!r}')
    return int(value)


def read_rows(paths, columns):
    rows=[]
    for path in paths:
        with open(path,newline='') as fh:
            reader=csv.DictReader(fh,delimiter='\t',quoting=csv.QUOTE_NONE)
            if reader.fieldnames!=columns:raise ValueError(f'Invalid table header: {path}')
            for row in reader:
                if None in row or any(v is None for v in row.values()):raise ValueError(f'Malformed table row: {path}')
                rows.append(row)
    return rows


def validate(genes, sites, audits, expected, threshold):
    keys=[tuple(key) for key in expected]
    if (len(set(keys))!=len(keys) or any(len(k)!=2 or not isinstance(k[0],str) or not k[0]
            or k[1] not in {'top','bottom'} for k in keys)):
        raise ValueError('Duplicate or invalid expected gene/direction keys')
    expected=set(keys)
    status={}
    for audit in audits:
        key=(audit['gene'],audit['direction'])
        if key in status or key not in expected or audit['status'] not in {'complete','excluded'}:
            raise ValueError(f'Duplicate, unexpected or failed audit: {key}')
        integer(audit.get('alignment_columns',''),1)
        retained=audit.get('retained_taxa')
        fg=audit.get('foreground_taxa')
        if (not isinstance(retained,list) or not isinstance(fg,list) or len(set(retained))!=len(retained)
                or not set(fg)<=set(retained)):
            raise ValueError('Invalid audited taxon sets')
        if audit['status']=='complete' and (len(retained)<2 or not fg):
            raise ValueError('Completed run is not biologically eligible')
        if audit['status']=='excluded':
            reason=audit.get('reason')
            if not ((reason=='fewer_than_two_taxa' and len(retained)<2)
                    or (reason=='no_foreground' and len(retained)>=2 and not fg)):
                raise ValueError('Invalid biological exclusion')
        status[key]=audit
    if set(status)!=expected:raise ValueError('Incomplete FADE audit coverage')
    by_gene={}
    for row in genes:
        key=(row['gene'],row['direction'])
        if key in by_gene or key not in status or status[key]['status']!='complete':
            raise ValueError(f'Duplicate or unexpected gene result: {key}')
        n=integer(row['n_sites_tested'],1)
        if n!=status[key]['alignment_columns']:raise ValueError(f'Site coverage mismatch: {key}')
        counts=[integer(row['sigcount_'+aa]) for aa in AA_LETTERS]
        if any(c>n for c in counts) or sum(counts)!=integer(row['n_sig_site_aa']):raise ValueError('Inconsistent site/AA counts')
        bf=number(row['max_bf'])
        hit=bf>=threshold
        if hit!=(sum(counts)>0):raise ValueError('Gene BF and hit counts disagree')
        if hit:
            if row['top_target_aa'] not in AA_LETTERS or not 1<=integer(row['top_site'],1)<=n:raise ValueError('Invalid top site/AA')
        elif row['top_target_aa']!='NA' or row['top_site']!='NA':raise ValueError('Unexpected top hit below threshold')
        by_gene[key]=row
    if set(by_gene)!={key for key,a in status.items() if a['status']=='complete'}:raise ValueError('Incomplete completed-gene coverage')
    seen,by_sites=set(),{}
    for row in sites:
        key=(row['gene'],row['direction']);pos=integer(row['position'])
        if key not in by_gene or (key,pos) in seen:raise ValueError('Duplicate or orphan significant site')
        seen.add((key,pos))
        if pos>=integer(by_gene[key]['n_sites_tested']) or integer(row['site'],1)!=pos+1:raise ValueError('Invalid site coordinate')
        if number(row['max_bf'])<threshold or number(row['max_bf'])>number(by_gene[key]['max_bf']):raise ValueError('Site BF disagrees with gene')
        if row['target_aa'] not in AA_LETTERS:raise ValueError('Invalid target residue')
        by_sites.setdefault(key,[]).append(row)
    for key,g in by_gene.items():
        rows=by_sites.get(key,[])
        n=len(rows);pairs=integer(g['n_sig_site_aa'])
        if n!=status[key]['n_significant_sites'] or not n<=pairs<=20*n:raise ValueError('Inconsistent significant-site coverage')
        if n:
            top=next((s for s in rows if s['site']==g['top_site']),None)
            if not top or number(top['max_bf'])!=number(g['max_bf']) or top['target_aa']!=g['top_target_aa']:
                raise ValueError('Top site and gene maximum disagree')
    return status


def write(path, columns, rows, delimiter='\t'):
    with open(path,'w',newline='') as fh:
        writer=csv.DictWriter(fh,columns,delimiter=delimiter,extrasaction='ignore',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--direction',required=True,choices=['top','bottom'])
    for name in ['genes','sites','audits']:p.add_argument('--'+name,nargs='*',default=[])
    p.add_argument('--expected',required=True);p.add_argument('--settings',required=True)
    a=p.parse_args();settings=json.loads(Path(a.settings).read_text())
    genes=read_rows(a.genes,GENES_COLUMNS);sites=read_rows(a.sites,SITES_COLUMNS)
    audits=[json.loads(Path(f).read_text()) for f in a.audits]
    status=validate(genes,sites,audits,json.loads(Path(a.expected).read_text()),settings['bf_threshold'])
    genes=sorted([g for g in genes if g['direction']==a.direction],key=lambda g:(-int(g['n_sig_site_aa']),-number(g['max_bf']),g['gene']))
    sites=sorted([s for s in sites if s['direction']==a.direction],key=lambda s:(-number(s['max_bf']),s['gene'],int(s['position'])))
    write(f'fade_summary_{a.direction}.tsv',SUMMARY,genes)
    write(f'fade_gene_counts_{a.direction}.tsv',SUMMARY+[f'sigcount_{aa}' for aa in AA_LETTERS],genes)
    write(f'fade_site_bf_{a.direction}.tsv',SITES,sites)
    write(f'fade_sites_{a.direction}.csv',SITES,sites,',')
    audit_rows=[dict(gene=k[0],direction=k[1],status=v['status'],reason=v.get('reason','')) for k,v in sorted(status.items()) if k[1]==a.direction]
    write('fade_audit.tsv',['gene','direction','status','reason'],audit_rows)
    metadata=dict(schema_version=1,complete=True,direction=a.direction,settings=settings,
                  expected=len(audit_rows),tested=len(genes),excluded=sum(r['status']=='excluded' for r in audit_rows))
    metadata['output_sha256']={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in
        [f'fade_gene_counts_{a.direction}.tsv',f'fade_sites_{a.direction}.csv','fade_audit.tsv']}
    foreground={v['input_sha256']['foreground'] for k,v in status.items() if k[1]==a.direction}
    if len(foreground)>1:raise ValueError('Conflicting foreground inputs')
    metadata['foreground_sha256']=next(iter(foreground),None)
    metadata['code_sha256']=Path('code.sha256').read_text().strip() if Path('code.sha256').exists() else None
    Path('fade_run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':main()
