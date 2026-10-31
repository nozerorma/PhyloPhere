#!/usr/bin/env python3
"""Streaming, transcript-aware joins of human candidate SNVs to fixed snapshots."""
import argparse
import hashlib
import math
import re
from vep_common import (BASE, columns, digest, finish, integer, output_row,
                        read_targets, same_id, source_manifest, table)

DEFAULTS = {
    'primateai':dict(chrom='chr',pos='pos',ref='non_flipped_ref',alt='non_flipped_alt',
                    ref_aa='ref_aa',alt_aa='alt_aa',aa_pos='change_position_1based',
                    transcript_id='gene_name',gene_id='ensembl_gene_id',
                    score='score_PAI3D',percentile='percentile_PAI3D',prediction='prediction'),
    'cosmic':dict(chrom='CHROMOSOME',pos='GENOME_START',ref='GENOMIC_WT_ALLELE',
                  alt='GENOMIC_MUT_ALLELE',aa_change='MUTATION_AA',
                  transcript_id='TRANSCRIPT_ACCESSION',gene_id='GENE_ID')}


def map_database(targets_file, database, manifest_file, output, source):
    data = read_targets(targets_file)
    if data.get('schema')!='phylophere.vep_targets.v1':
        raise ValueError('Invalid prepared VEP targets')
    targets = data['targets']
    manifest = source_manifest(database,manifest_file,source)
    fields = DEFAULTS[source] | manifest.get('fields',{})
    header = columns(database)
    if not set(fields.values()) <= set(header):
        raise ValueError(f'Database lacks declared fields: {set(fields.values())-set(header)}')
    if source=='cosmic' and (not isinstance(manifest.get('scope_gene_ids'),list) or not manifest['scope_gene_ids']):
        raise ValueError('COSMIC snapshot requires explicit nonempty scope_gene_ids (Census coverage)')
    lookup, hits = {}, {i:[] for i in range(len(targets))}
    for i,t in enumerate(targets):
        if t['status']!='ready':continue
        if t['reference']['annotation_release']!=manifest['annotation_release']:
            raise ValueError('Annotation database and reference releases differ')
        for v in t['variants']:
            lookup.setdefault((v['chrom'],v['pos'],v['ref'],v['alt']),[]).append((i,v))
    seen = set()
    predictions = {}
    for row in table(database):
        chrom = row[fields['chrom']]
        if not chrom.startswith('chr'):chrom='chr'+chrom
        if chrom=='chrMT':chrom='chrM'
        pos = integer(row[fields['pos']])
        ref,alt = row[fields['ref']],row[fields['alt']]
        # Other variant classes are valid source rows, but outside this SNV contract.
        matches = lookup.get((chrom,pos,ref,alt))
        if not matches:continue
        if source=='cosmic':
            m = re.fullmatch(r'p\.([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])',row[fields['aa_change']])
            if not m:continue
            wt,aa,mut=m.group(1),int(m.group(2)),m.group(3)
        else:
            wt,mut=row[fields['ref_aa']],row[fields['alt_aa']]
            aa=integer(row[fields['aa_pos']])
        for i,v in matches:
            t=targets[i];r=t['reference']
            aliases = r.get('transcript_aliases',[])
            tx=row[fields['transcript_id']]
            if not any(same_id(tx,x) for x in [r['transcript_id']]+aliases):continue
            gid=row[fields['gene_id']]
            if not same_id(gid,r['gene_id']):continue
            if aa!=t['reference_aa_position'] or wt!=v['ref_aa'] or mut!=v['alt_aa']:
                raise ValueError('Source transcript/variant matches but protein consequence disagrees with verified reference')
            if source=='primateai':
                score,percentile=float(row[fields['score']]),float(row[fields['percentile']])
                if not math.isfinite(score) or not 0<=score<=1 or not math.isfinite(percentile) or not 0<=percentile<=1:
                    raise ValueError('Invalid PrimateAI score/percentile (expected [0,1])')
                if row[fields['prediction']].lower() not in {'benign','pathogenic'}:
                    raise ValueError('Unrecognised PrimateAI gene-specific prediction')
            if source=='primateai':
                signature=(float(row[fields['score']]),float(row[fields['percentile']]),row[fields['prediction']].lower())
                pk=(i,v['variant_id'])
                if pk in predictions and predictions[pk]!=signature:
                    raise ValueError('Conflicting PrimateAI predictions for the same verified variant/transcript')
                predictions[pk]=signature
            fingerprint=hashlib.sha256('\t'.join(row[k] for k in header).encode()).hexdigest()
            key=(i,v['variant_id'],fingerprint)
            if key in seen:continue
            seen.add(key)
            extra={'db_'+k:value for k,value in row.items()}
            if source=='primateai':
                extra.update(score_PAI3D=row[fields['score']],percentile_PAI3D=row[fields['percentile']],
                             prediction=row[fields['prediction']].lower())
            extra.update(source_record_sha256=fingerprint,variant_id=v['variant_id'],orientation=v['orientation'])
            hits[i].append(extra)
    rows=[]
    for i,t in enumerate(targets):
        if t['status']!='ready':
            rows.append(output_row(t,source));continue
        scope = source!='cosmic' or any(same_id(g,t['reference_gene_id']) for g in manifest['scope_gene_ids'])
        if not scope:
            if hits[i]:raise ValueError('COSMIC hit outside declared Census scope')
            rows.append(output_row(t,source,'outside_source_scope'));continue
        covered=len({h['variant_id'] for h in hits[i]})
        if hits[i]:
            coverage='complete' if covered==len(t['variants']) and not t.get('unqueryable_alt_aas') else 'partial'
            for h in hits[i]:rows.append(output_row(t,source,'annotated',coverage=coverage,annotated_count=covered,**h))
        else:
            # No source record is not a benign prediction or absence of cancer relevance.
            status='no_hit' if source=='cosmic' else 'missing_annotation'
            rows.append(output_row(t,source,status,coverage='snapshot_searched' if source=='cosmic' else 'unavailable'))
    extras=(['score_PAI3D','percentile_PAI3D','prediction'] if source=='primateai' else [])+['source_record_sha256']+['db_'+k for k in header]
    finish(output,BASE+extras,rows,source,
           dict(targets_sha256=digest(targets_file),database_sha256=manifest['sha256'],manifest_sha256=manifest['_sha256']),
           {k:manifest[k] for k in ['version','source_uri','annotation_release','species','assembly']})


def main(source):
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['targets','database','manifest','output']:p.add_argument('--'+arg,required=True)
    a=p.parse_args()
    map_database(a.targets,a.database,a.manifest,a.output,source)
