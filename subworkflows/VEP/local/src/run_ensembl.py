#!/usr/bin/env python3
"""Run version-pinned VEP on genomic SNVs; never install or download resources."""
import argparse
import re
import shutil
import subprocess
from pathlib import Path
from build_vep_variants import build
from join_vep_output import join
from vep_common import BASE, finish, integer, output_row, read_targets, write_json

VEP_COLUMNS=['Uploaded_variation','Location','Allele','Ensembl_Gene','Feature','Feature_type',
             'Consequence','Protein_position','Amino_acids','Extra']


def resource(cache,release):
    """Use the standard Ensembl cache layout, without a private manifest."""
    root=Path(cache).resolve()
    if not (root/'homo_sapiens'/f'{release}_GRCh38').is_dir():
        raise ValueError(f'Ensembl cache {release}/GRCh38 is absent; prepare it separately')
    fastas=[p for p in root.rglob('*') if p.is_file() and p.name.endswith(('.fa','.fasta','.fa.gz','.fasta.gz')) and Path(str(p)+'.fai').is_file()]
    if len(fastas)!=1:raise ValueError('Offline VEP needs one preindexed genomic FASTA in its cache directory')
    return fastas[0]


def run(a):
    from resolve_ensembl import resolve
    data=read_targets(a.targets)
    inputs={'targets':str(a.targets),'source_fingerprint':a.fingerprint}
    meta={'mode':a.mode,'species':'homo_sapiens','assembly':'GRCh38'}
    executable=shutil.which('vep')
    ready=any(t['status']=='ready' for t in data['targets'])
    if ready and not executable:raise RuntimeError('Ensembl VEP executable missing')
    if ready and a.mode=='database':data=resolve(data)
    for t in data['targets']:
        if t['status']=='ready' and not t['variants']:
            t.update(status='unmapped_codon',detail='Exact codon bases unavailable for offline VCF; database mode can resolve the declared protein')
    data,variants=build(a.targets,'candidate_variants.vcf',data)
    if not variants:
        finish(a.output,BASE+VEP_COLUMNS,[output_row(t,'ensembl') for t in data['targets']],
               'ensembl',inputs,meta|{'execution':'not_needed_no_resolved_candidate_snvs'})
        return
    p=subprocess.run([executable,'--help'],capture_output=True,text=True,timeout=60)
    Path('vep_version.txt').write_text(p.stdout+p.stderr)
    match=re.search(r'ensembl-vep\s*:\s*(\d+)',p.stdout+p.stderr,re.I)
    if not match:raise ValueError('Cannot determine the installed Ensembl VEP version')
    release=int(match[1]);meta['release']=release
    cmd=[executable,'--input_file','candidate_variants.vcf','--format','vcf',
         '--output_file','vep_tab_output.txt','--tab','--species','homo_sapiens','--assembly','GRCh38',
         '--canonical','--transcript_version','--protein','--check_ref','--force_overwrite','--no_stats',
         '--warning_file','vep_warnings.txt','--no_check_variants_order']
    if release>=105:cmd+=['--mane']
    if a.mode=='offline':
        fasta=resource(a.cache,release)
        cmd+=['--offline','--cache','--cache_version',str(release),'--dir_cache',str(Path(a.cache).resolve()),'--fasta',str(fasta)]
    else:
        if not re.fullmatch(r'[A-Za-z0-9.-]+',a.host):raise ValueError('Invalid Ensembl database host')
        port=integer(a.port)
        if port>65535:raise ValueError('Invalid Ensembl database port')
        cmd+=['--database','--host',a.host,'--port',str(port),'--user','anonymous','--db_version',str(release)]
        meta.update(host=a.host,port=port,cache_reusable=False)
    meta['command']=cmd
    with open('vep_stdout.log','w') as out,open('vep_stderr.log','w') as err:
        p=subprocess.run(cmd,stdout=out,stderr=err)
    if p.returncode:raise RuntimeError(f'Ensembl VEP failed with exit {p.returncode}; see vep_stderr.log')
    if re.search(r'\b(?:ERROR|EXCEPTION)\b',Path('vep_stderr.log').read_text(),re.I):
        raise RuntimeError('Ensembl VEP reported an error/warning on stderr; see vep_stderr.log')
    warnings=Path('vep_warnings.txt')
    meta['warnings_present']=warnings.is_file() and bool(warnings.read_text().strip())
    join('vep_tab_output.txt',data,variants,a.output,inputs,meta)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--targets',required=True);p.add_argument('--output',required=True)
    p.add_argument('--fingerprint',required=True);p.add_argument('--cache',default='')
    p.add_argument('--mode',choices=['database','offline'],required=True)
    p.add_argument('--host',default='ensembldb.ensembl.org');p.add_argument('--port',default='3306')
    a=p.parse_args()
    try:run(a)
    except Exception as exc:
        Path(a.output).unlink(missing_ok=True)
        write_json('ensembl_failure.json',dict(status='technical_failure',error=str(exc)))
        raise

if __name__=='__main__':main()
