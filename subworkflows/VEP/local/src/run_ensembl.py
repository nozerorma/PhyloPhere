#!/usr/bin/env python3
"""Run version-pinned VEP on genomic SNVs; never install or download resources."""
import argparse
import re
import shutil
import subprocess
from pathlib import Path
from build_vep_variants import build
from join_vep_output import join
from vep_common import BASE, digest, finish, integer, output_row, read_json, write_json

VEP_COLUMNS=['Uploaded_variation','Location','Allele','Ensembl_Gene','Feature','Feature_type',
             'Consequence','Protein_position','Amino_acids','Extra']


def resource(cache,release):
    root=Path(cache).resolve(); manifest=root/'vep_resource.json'
    m=read_json(manifest)
    if m.get('schema')!='phylophere.vep_cache.v1' or m.get('species')!='homo_sapiens' or m.get('assembly')!='GRCh38' or m.get('release')!=release:
        raise ValueError('Incompatible offline cache manifest')
    if not m.get('source_uri') or not isinstance(m.get('files'),dict) or not m['files']:
        raise ValueError('Offline cache requires provenance and file checksums')
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=manifest}
    if actual!=set(m['files']):
        raise ValueError('Offline cache has missing or unmanifested files; prepare it before annotation')
    for rel,sha in m['files'].items():
        p=(root/rel).resolve()
        if not p.is_relative_to(root) or not p.is_file() or digest(p)!=sha:
            raise ValueError(f'Offline cache checksum/path mismatch: {rel}')
    fasta=m.get('fasta')
    if fasta not in m['files'] or not Path(str(root/fasta)+'.fai').is_file():
        raise ValueError('Offline cache needs a checksummed, preindexed genomic FASTA')
    if not (root/'homo_sapiens'/f'{release}_GRCh38').is_dir():
        raise ValueError('Expected version/assembly cache directory is missing')
    return root/fasta,digest(manifest)


def run(a):
    data,variants=build(a.targets,'candidate_variants.vcf')
    inputs={'targets_sha256':digest(a.targets),'content_fingerprint':a.fingerprint}
    meta={'mode':a.mode,'species':'homo_sapiens','assembly':'GRCh38'}
    if not variants:
        finish(a.output,BASE+VEP_COLUMNS,[output_row(t,'ensembl') for t in data['targets']],
               'ensembl',inputs,meta|{'execution':'not_needed_no_candidate_snvs'})
        return
    releases={t['reference']['annotation_release'] for t in data['targets'] if t['status']=='ready'}
    if len(releases)!=1:raise ValueError('All VEP references must share one Ensembl release')
    release=releases.pop()
    executable=shutil.which('vep')
    if not executable:raise RuntimeError('Ensembl VEP executable missing')
    p=subprocess.run([executable,'--help'],capture_output=True,text=True,timeout=60)
    Path('vep_version.txt').write_text(p.stdout+p.stderr)
    match=re.search(r'ensembl-vep\s*:\s*(\d+)',p.stdout+p.stderr,re.I)
    if not match or int(match[1])!=release:
        raise ValueError(f'VEP executable and reference release must agree ({release})')
    meta.update(release=release,executable_sha256=digest(executable),version_output_sha256=digest('vep_version.txt'))
    cmd=[executable,'--input_file','candidate_variants.vcf','--format','vcf',
         '--output_file','vep_tab_output.txt','--tab','--species','homo_sapiens','--assembly','GRCh38',
         '--canonical','--transcript_version','--protein','--check_ref','--force_overwrite','--no_stats',
         '--warning_file','vep_warnings.txt','--no_check_variants_order']
    if release>=105:cmd+=['--mane']
    if a.mode=='offline':
        fasta,sha=resource(a.cache,release)
        inputs['cache_manifest_sha256']=sha
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
    if re.search(r'\b(?:WARNING|ERROR|EXCEPTION)\b',Path('vep_stderr.log').read_text(),re.I):
        raise RuntimeError('Ensembl VEP reported an error/warning on stderr; see vep_stderr.log')
    warnings=Path('vep_warnings.txt')
    if warnings.is_file() and warnings.read_text().strip():
        raise RuntimeError('Ensembl VEP reported warnings/skipped or incompatible variants; see vep_warnings.txt')
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
