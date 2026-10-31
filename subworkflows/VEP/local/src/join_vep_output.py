#!/usr/bin/env python3
"""Validate the VEP VCF-ID join, retaining the verified transcript only.

Different consequences for other transcripts are outside the transfer contract.
VEP rows absent for a ready candidate remain unresolved, never negative results.
"""
import re
from pathlib import Path
from vep_common import BASE, digest, finish, output_row, same_id


def join(output_file,data,variants,output,inputs,metadata):
    header=None; records=[]; seen=set(); encountered=set()
    with open(output_file) as f:
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#'):
                if header is not None:raise ValueError('Repeated VEP output header')
                header=line[1:].rstrip('\r\n').split('\t')
                if len(header)!=len(set(header)) or not {'Uploaded_variation','Location','Allele','Gene','Feature','Feature_type','Consequence','Protein_position','Amino_acids','Extra'}<=set(header):
                    raise ValueError('Invalid VEP output schema')
                continue
            if header is None:raise ValueError('VEP output has no header')
            fields=line.rstrip('\r\n').split('\t')
            if len(fields)!=len(header):raise ValueError('Ragged VEP output')
            row=dict(zip(header,fields));vid=row['Uploaded_variation']
            if vid not in variants:raise ValueError(f'Unexpected VEP variant {vid}')
            v=variants[vid]
            loc=row['Location'].removeprefix('chr')
            if loc not in {f"{v['chrom'].removeprefix('chr')}:{v['pos']}",f"{v['chrom'].removeprefix('chr')}:{v['pos']}-{v['pos']}"} or row['Allele']!=v['alt']:
                raise ValueError('VEP location/allele differs from candidate SNV')
            encountered.add(vid)
            key=tuple(fields)
            if key not in seen:records.append(row);seen.add(key)
    if header is None:raise ValueError('VEP output has no header')
    lookup={}
    for r in records:lookup.setdefault(r['Uploaded_variation'],[]).append(r)
    rows=[]
    for t in data['targets']:
        if t['status']!='ready':rows.append(output_row(t,'ensembl'));continue
        matches=[];covered=set()
        for v in t['variants']:
            for r in lookup.get(v['variant_id'],[]):
                if r['Feature_type']!='Transcript' or r['Feature']!=t['reference_transcript_id']:continue
                extra_fields=dict(x.split('=',1) for x in r['Extra'].split(';') if '=' in x)
                protein=r.get('ENSP') or extra_fields.get('ENSP','')
                if not protein or not same_id(protein,t['reference_protein_id']):
                    raise ValueError('VEP protein identity differs from verified reference')
                if not same_id(r['Gene'],t['reference_gene_id']):raise ValueError('VEP transcript has incompatible gene')
                if r['Protein_position']!=str(t['reference_aa_position']) or r['Amino_acids']!=v['ref_aa']+'/'+v['alt_aa']:
                    raise ValueError('VEP protein consequence differs from verified reference')
                terms=r['Consequence'].split(',')
                if 'missense_variant' not in terms:raise ValueError('Candidate missense has incompatible VEP consequence')
                extra=dict(r);extra['Ensembl_Gene']=extra.pop('Gene')
                extra.update(variant_id=v['variant_id'],orientation=v['orientation'])
                matches.append(extra);covered.add(v['variant_id'])
        if matches:
            coverage='complete' if len(covered)==len(t['variants']) and not t.get('unqueryable_alt_aas') else 'partial'
            for r in matches:rows.append(output_row(t,'ensembl','annotated',coverage=coverage,annotated_count=len(covered),**r))
        else:rows.append(output_row(t,'ensembl','unresolved_annotation',coverage='unavailable',detail='No compatible verified-transcript consequence'))
    finish(output,BASE+['Ensembl_Gene' if k=='Gene' else k for k in header],rows,'ensembl',
           inputs|{'raw_output_sha256':digest(output_file)},metadata|{'returned_variants':len(encountered)})
