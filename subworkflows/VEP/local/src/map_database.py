#!/usr/bin/env python3
"""Read original PrimateAI/COSMIC tables without private snapshot prerequisites.

A match is descriptive human variant evidence at a MAP locus/residue. Source
transcripts/proteins are retained, checked when the input declares their identity,
and never relabelled as an ASR event in a nonhuman species.
"""
import argparse
import math
import re
from collections import defaultdict
from pathlib import Path
from vep_common import (BASE, candidate_orientation, columns, finish, integer,
                        output_row, read_json, read_targets, same_id, table)

FIELDS={
 'primateai':dict(chrom='chr',pos='pos',ref='non_flipped_ref',alt='non_flipped_alt',
    wt='ref_aa',mut='alt_aa',aa_pos='change_position_1based',transcript='gene_name',score='score_PAI3D'),
 'cosmic':dict(chrom='CHROMOSOME',pos='GENOME_START',ref='GENOMIC_WT_ALLELE',alt='GENOMIC_MUT_ALLELE',
    change='MUTATION_AA',transcript='TRANSCRIPT_ACCESSION')}


def map_database(targets_file,database,manifest_file,output,source):
    data=read_targets(targets_file);targets=data['targets']
    # Optional metadata records source version/field aliases; it is not a run prerequisite.
    meta=read_json(manifest_file) if manifest_file and Path(manifest_file).is_file() else {}
    if meta.get('assembly') and meta['assembly']!='GRCh38':raise ValueError('Annotation database assembly differs from MAP hg38 coordinates')
    fields=FIELDS[source] | meta.get('fields',{})
    header=columns(database)
    missing=set(fields.values())-set(header)
    if missing:raise ValueError(f'{source} database lacks required columns: {sorted(missing)}')
    lookup=defaultdict(list);hits=defaultdict(list);seen=set();predictions={}
    for i,t in enumerate(targets):
        if t['status']=='ready':
            for coord in t['genomic_positions']:lookup[coord].append(i)
    scanned=0;skipped_unlocated=0
    # No viable targets: validate the schema, but do not decompress a multi-GB database.
    for row in (table(database) if lookup else []):
        scanned+=1
        chrom=row[fields['chrom']]
        if not chrom.startswith('chr'):chrom='chr'+chrom
        if chrom=='chrMT':chrom='chrM'
        if row[fields['pos']] in {'','NA','-'}:
            skipped_unlocated+=1;continue
        pos=integer(row[fields['pos']]);coord=f'{chrom}:{pos}'
        if coord not in lookup:continue
        ref,alt=row[fields['ref']],row[fields['alt']]
        if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt:continue
        if source=='primateai':
            wt,mut=row[fields['wt']],row[fields['mut']];aa=integer(row[fields['aa_pos']])
        else:
            match=re.fullmatch(r'p\.([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])',row[fields['change']])
            if not match:continue
            wt,aa,mut=match[1],int(match[2]),match[3]
        vid=f'{coord}:{ref}:{alt}'
        for i in lookup[coord]:
            t=targets[i]
            if aa!=t['reference_aa_position']:continue
            if t.get('hg38_ref_aa') and wt!=t['hg38_ref_aa']:continue
            orientation=candidate_orientation(t,wt,mut)
            if not orientation:continue
            tx=row[fields['transcript']]
            if t.get('reference_transcript_id') and not same_id(tx,t['reference_transcript_id']):continue
            if source=='cosmic' and t.get('reference_protein_id') and row.get('HGVSP'):
                pid=row['HGVSP'].split(':',1)[0]
                if pid.startswith('ENSP') and not same_id(pid,t['reference_protein_id']):continue
            # With complete codon bases, also enforce the exact nucleotide change.
            if t['candidate_count_known'] and vid not in {v['variant_id'] for v in t['variants']}:continue
            extra={'db_'+k:v for k,v in row.items()}
            if source=='primateai':
                score=float(row[fields['score']])
                if not math.isfinite(score) or not 0<=score<=1:raise ValueError('Invalid PrimateAI score')
                percentile=row.get('percentile_PAI3D','')
                if percentile and (not math.isfinite(float(percentile)) or not 0<=float(percentile)<=1):raise ValueError('Invalid PrimateAI percentile')
                prediction=row.get('prediction','').lower()
                if prediction and prediction not in {'benign','pathogenic'}:raise ValueError('Unknown PrimateAI prediction')
                signature=(score,percentile,prediction)
                pk=(i,vid,tx)
                if pk in predictions and predictions[pk]!=signature:raise ValueError('Conflicting PrimateAI records for a variant/transcript')
                predictions[pk]=signature
                extra.update(score_PAI3D=row[fields['score']],percentile_PAI3D=percentile,prediction=prediction)
            raw=tuple(row[k] for k in header)
            key=(i,vid,raw)
            if key in seen:continue
            seen.add(key)
            extra.update(variant_id=vid,orientation=orientation,hg38_ref_aa=wt,caas_alt_aas=mut)
            hits[i].append(extra)
    rows=[]
    for i,t in enumerate(targets):
        if t['status']!='ready':rows.append(output_row(t,source));continue
        if hits[i]:
            n=len({h['variant_id'] for h in hits[i]})
            complete=t['candidate_count_known'] and n==len(t['variants']) and not t.get('unqueryable_alt_aas')
            for h in hits[i]:
                rows.append(output_row(t,source,'annotated',coverage='complete' if complete else 'partial',
                    candidate_count=len(t['variants']) if t['candidate_count_known'] else n,annotated_count=n,**h))
        else:
            # Absence in a source with unknown biological/isoform coverage is not a negative call.
            rows.append(output_row(t,source,'missing_annotation',coverage='searched_map_coordinates',
                        detail='No compatible source record at declared MAP coordinates; source coverage is not assumed'))
    extras=(['score_PAI3D','percentile_PAI3D','prediction'] if source=='primateai' else [])+['db_'+k for k in header]
    finish(output,BASE+extras,rows,source,dict(targets=str(targets_file),database=str(database)),
           dict(version=meta.get('version','not_recorded'),source_uri=meta.get('source_uri',str(database)),
                scanned_rows=scanned,unlocated_source_rows=skipped_unlocated,
                limitation='MAP anchors alone cover only the declared first base; inferred complete codons require adjacent reference anchors or explicit coordinates.'))


def main(source):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['targets','database','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--manifest')
    a=p.parse_args();map_database(a.targets,a.database,a.manifest,a.output,source)
