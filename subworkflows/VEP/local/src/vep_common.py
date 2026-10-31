#!/usr/bin/env python3
"""Strict literal TSV IO and verified human-reference projection for VEP.

SCORING's pooled residue sets describe candidate states, never paired ASR events.
All genomic coordinates here are one-based, forward-genome coordinates.
"""
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

AA = set('ACDEFGHIKLMNPQRSTVWY')
MISSING = {'', 'NA', '-'}
BASE = ['Gene', 'Position', 'side', 'status', 'coverage', 'detail', 'source',
        'species', 'assembly', 'reference_gene_id', 'reference_transcript_id',
        'reference_protein_id', 'reference_aa_position', 'hg38_ref_aa',
        'caas_alt_aas', 'caas_change', 'transfer', 'orientation', 'evidence_basis',
        'variant_id', 'candidate_count', 'annotated_count', 'unqueryable_alt_aas']
# NCBI standard genetic code, in T,C,A,G order at each base.
_CODONS = [a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
CODE = dict(zip(_CODONS, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    def pairs(items):
        out = {}
        for k, v in items:
            if k in out:
                raise ValueError(f'Duplicate JSON key {k} in {path}')
            out[k] = v
        return out
    with open(path) as f:
        return json.load(f, object_pairs_hook=pairs)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def open_text(path):
    return gzip.open(path, 'rt', newline='') if str(path).endswith('.gz') else open(path, newline='')


def table(path):
    """Yield literal fields. No quoting, NA conversion, ID stripping or inference."""
    with open_text(path) as f:
        header = f.readline().rstrip('\r\n').split('\t')
        if not all(header) or len(header) != len(set(header)):
            raise ValueError(f'Empty or duplicate header in {path}')
        for n, line in enumerate(f, 2):
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != len(header):
                raise ValueError(f'Ragged TSV {path}:{n}')
            yield dict(zip(header, fields))


def columns(path):
    with open_text(path) as f:
        names = f.readline().rstrip('\r\n').split('\t')
    if not all(names) or len(names) != len(set(names)):
        raise ValueError(f'Empty or duplicate header in {path}')
    return names


def require(path, names):
    missing = set(names)-set(columns(path))
    if missing:
        raise ValueError(f'Missing columns in {path}: {sorted(missing)}')


def write_table(path, header, rows):
    if len(header) != len(set(header)):
        raise ValueError('Duplicate output columns')
    with open(path, 'w') as out:
        out.write('\t'.join(header)+'\n')
        for row in rows:
            values = [str(row.get(k, '')) for k in header]
            if any(any(c in v for c in '\t\r\n') for v in values):
                raise ValueError('Tabs/newlines forbidden in literal TSV fields')
            out.write('\t'.join(values)+'\n')


def integer(value, minimum=1):
    if isinstance(value, bool) or not re.fullmatch(r'0|[1-9][0-9]*', str(value)) or int(value)<minimum:
        raise ValueError(f'Invalid integer coordinate/count: {value!r}')
    return int(value)


def residue_set(text):
    if not text:
        return set()
    values = text.split(',')
    if len(values) != len(set(values)) or any(len(v)!=1 or v not in AA | set('XBZJUO*-?') for v in values):
        raise ValueError(f'Invalid residue set {text!r}')
    return set(values)


def anc_der_from_row(ancestral_aa, derived_aa, caas=None, side=None):
    """No ancestral inference from the opposite phenotype group."""
    return residue_set(ancestral_aa), residue_set(derived_aa)


def load_scores(path):
    require(path, ['Gene','Position','side','CAAS_score','ancestral_aa','derived_aa','caas'])
    seen = set()
    rows = []
    for row in table(path):
        gene = row['Gene']
        if not gene or gene != gene.strip() or '/' in gene or '\\' in gene or gene in {'.','..'}:
            raise ValueError(f'Invalid literal gene ID {gene!r}')
        pos = integer(row['Position'], 0)
        if row['side'] not in {'top','bottom'}:
            raise ValueError('SCORING side must be top or bottom')
        key = (gene, pos, row['side'])
        if key in seen:
            raise ValueError(f'Duplicate SCORING key {key}')
        seen.add(key)
        score = float(row['CAAS_score'])
        if not math.isfinite(score) or not 0<score<=1:
            raise ValueError('SCORING score must be finite and in (0,1]')
        residue_set(row['ancestral_aa']); residue_set(row['derived_aa'])
        rows.append(dict(row, Position=pos))
    return rows


def map_path(gene, directory):
    # Literal matching: dots, quotes, brackets and glob metacharacters are IDs.
    paths = [p for p in Path(directory).iterdir() if p.is_file() and
             (p.name == gene+'.map.tsv' or (p.name.startswith(gene+'.') and p.name.endswith('.map.tsv')))]
    if len(paths)>1:
        raise ValueError(f'Multiple MAP references for {gene}: {[p.name for p in paths]}')
    return paths[0] if paths else None


def load_map(path):
    require(path, ['status','prot_ali_col','hg38_aa_pos','hg38_nt_coord'])
    result, residues = {}, set()
    for r in table(path):
        if r['status'] not in {'selected','removed'}:
            raise ValueError(f'Unknown MAP status in {path}')
        if r['status'] != 'selected':
            if r['prot_ali_col'] not in MISSING:
                raise ValueError('Removed MAP row carries a trimmed column')
            continue
        col = integer(r['prot_ali_col'])-1
        if col in result:
            raise ValueError('Repeated MAP alignment column')
        aa = r['hg38_aa_pos']
        aa = None if aa in MISSING else integer(aa)
        if aa is not None and aa in residues:
            raise ValueError('Repeated MAP reference residue')
        if aa is not None:
            residues.add(aa)
        coord = r['hg38_nt_coord']
        if coord not in MISSING and not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):[1-9][0-9]*',coord):
            raise ValueError('Invalid MAP genomic coordinate')
        result[col] = {'aa_position':aa, 'anchor':None if coord in MISSING else coord}
    if sorted(result) != list(range(len(result))):
        raise ValueError('Selected MAP columns must be contiguous 1..n')
    return result


def identity(value, prefix):
    if not isinstance(value, str) or not re.fullmatch(prefix+r'[0-9]+\.[1-9][0-9]*',value):
        raise ValueError(f'Versioned {prefix} identity required: {value!r}')
    return value


def same_id(actual, expected):
    """Stable IDs may omit a version; supplied versions must agree."""
    return actual == expected or ('.' not in actual and actual == expected.rsplit('.',1)[0])


def load_reference(path, gene, mapping):
    meta_path = Path(str(path)+'.vep.json')
    if not meta_path.is_file():
        return None
    m = read_json(meta_path)
    if m.get('schema')!='phylophere.vep_reference.v1' or m.get('gene')!=gene:
        raise ValueError('Invalid VEP reference manifest identity/schema')
    if m.get('species')!='homo_sapiens' or m.get('assembly')!='GRCh38':
        raise ValueError('VEP reference must be homo_sapiens/GRCh38')
    if m.get('map_sha256') != digest(path):
        raise ValueError('MAP checksum differs from verified reference')
    integer(m.get('annotation_release'))
    for k, prefix in [('gene_id','ENSG'),('transcript_id','ENST'),('protein_id','ENSP')]:
        identity(m.get(k),prefix)
    if m.get('translation_table',1) != 1:
        raise ValueError('Only nuclear standard-code references are supported')
    if m.get('strand') not in {'+','-'} or not m.get('source_uri'):
        raise ValueError('Reference strand and provenance source_uri required')
    sequence = m.get('protein_sequence','')
    if not sequence or set(sequence)-AA or hashlib.sha256(sequence.encode()).hexdigest()!=m.get('protein_sha256'):
        raise ValueError('Invalid reference protein sequence/checksum')
    aliases = m.get('transcript_aliases', [])
    if not isinstance(aliases, list) or any(not isinstance(x,str) or not re.fullmatch(r'(?:ENST[0-9]+|N[MR]_[0-9]+)\.[1-9][0-9]*',x) for x in aliases) or len(aliases)!=len(set(aliases)):
        raise ValueError('Invalid versioned transcript aliases')
    codons = m.get('codons')
    if not isinstance(codons, dict):
        raise ValueError('Reference codons object required')
    previous = None
    seen_nt = set()
    for key,c in sorted(codons.items(),key=lambda item:integer(item[0])):
        aa = integer(key)
        positions = c.get('positions',[])
        triplet = c.get('coding_bases','')
        if aa>len(sequence) or CODE.get(triplet)!=sequence[aa-1] or len(positions)!=3:
            raise ValueError('Codon translation differs from reference protein')
        if any(not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):[1-9][0-9]*',x) for x in positions):
            raise ValueError('Invalid codon coordinates')
        chroms = {x.split(':')[0] for x in positions}
        nums = [int(x.split(':')[1]) for x in positions]
        sign = 1 if m['strand']=='+' else -1
        if len(chroms)!=1 or any((b-a)*sign<=0 for a,b in zip(nums,nums[1:])):
            raise ValueError('Codon coordinate order disagrees with strand')
        if previous and (positions[0].split(':')[0]!=previous[0] or (nums[0]-previous[1])*sign<=0):
            raise ValueError('Reference codons not ordered on a single chromosome')
        previous = (positions[-1].split(':')[0],nums[-1])
        if seen_nt.intersection(positions):
            raise ValueError('Repeated nucleotide across reference codons')
        seen_nt.update(positions)
    for r in mapping.values():
        c = codons.get(str(r['aa_position']))
        if c and r['anchor'] and c['positions'][0]!=r['anchor']:
            raise ValueError('MAP anchor differs from verified codon')
    m['_manifest_sha256'] = digest(meta_path)
    return m


def variants_for(ref, aa_position, alternatives):
    c = ref['codons'].get(str(aa_position))
    if not c:
        return []
    result = []
    for i, base in enumerate(c['coding_bases']):
        for alt in 'ACGT':
            triplet = c['coding_bases'][:i]+alt+c['coding_bases'][i+1:]
            alt_aa = CODE[triplet]
            if alt==base or alt_aa not in alternatives:
                continue
            chrom, pos = c['positions'][i].split(':')
            wt, mut = (base,alt) if ref['strand']=='+' else (base.translate(COMPLEMENT),alt.translate(COMPLEMENT))
            vid = f'{chrom}:{pos}:{wt}:{mut}'
            result.append(dict(variant_id=vid,chrom=chrom,pos=int(pos),ref=wt,alt=mut,
                               ref_aa=ref['protein_sequence'][aa_position-1],alt_aa=alt_aa))
    return result


def prepare(scores, map_dir, proteins=None):
    rows = load_scores(scores)
    proteins_map = None
    if proteins:
        require(proteins,['gene','human_protein_id'])
        proteins_map = {}
        for r in table(proteins):
            if not r['gene'] or r['gene'] in proteins_map:
                raise ValueError('Empty or duplicate gene in protein mapping')
            proteins_map[r['gene']] = r['human_protein_id']
    cache, targets = {}, []
    for row in rows:
        gene, position, side = row['Gene'], row['Position'], row['side']
        t = dict(Gene=gene,Position=position,side=side,status='',coverage='not_evaluated',detail='',
                 species='homo_sapiens',assembly='GRCh38',caas_change=row['caas'],
                 transfer='human_reference_candidate',evidence_basis='pooled_ASR_state_candidate',
                 variants=[],reference=None)
        if gene not in cache:
            path = map_path(gene,map_dir)
            mapping = load_map(path) if path else {}
            ref = load_reference(path,gene,mapping) if path else None
            cache[gene] = (mapping,ref)
        mapping, ref = cache[gene]
        if mapping and position not in mapping:
            raise ValueError(f'SCORING position outside selected MAP: {gene}:{position}')
        anc, der = anc_der_from_row(row['ancestral_aa'],row['derived_aa'])
        if not anc or not der:
            t['status'] = 'unresolved_asr'
        elif (anc|der)-AA:
            t['status'] = 'excluded_gap' if (anc|der)&set('-?') else 'excluded_ambiguous_residue'
        elif not mapping:
            t['status'] = 'unmapped_no_map'
        elif position not in mapping:
            raise ValueError(f'SCORING position outside selected MAP: {gene}:{position}')
        elif mapping[position]['aa_position'] is None:
            t['status'] = 'unmapped_reference_gap'
        elif ref is None:
            t['status'] = 'unverified_reference'
        else:
            aa = mapping[position]['aa_position']
            if aa>len(ref['protein_sequence']):
                raise ValueError('MAP residue outside reference protein')
            if proteins_map is not None and (gene not in proteins_map or proteins_map[gene] in MISSING):
                t['status'] = 'unmapped_protein_id'
            elif proteins_map is not None and not same_id(proteins_map[gene],ref['protein_id']):
                raise ValueError(f'Protein mapping disagrees with reference: {gene}')
            else:
                wt = ref['protein_sequence'][aa-1]
                # Pooled sets cannot recover domain-paired substitutions. Label every candidate.
                direct, reverse = (der-{wt} if wt in anc else set()), (anc-{wt} if wt in der else set())
                alternatives = direct|reverse
                t.update(reference=ref,reference_gene_id=ref['gene_id'],reference_transcript_id=ref['transcript_id'],
                         reference_protein_id=ref['protein_id'],reference_aa_position=aa,hg38_ref_aa=wt,
                         caas_alt_aas=','.join(sorted(alternatives)))
                if not alternatives:
                    t['status'] = 'reference_state_mismatch' if wt not in anc|der else 'excluded_no_change'
                elif str(aa) not in ref['codons']:
                    t['status'] = 'unmapped_codon'
                else:
                    variants = variants_for(ref,aa,alternatives)
                    reachable = {v['alt_aa'] for v in variants}
                    t['unqueryable_alt_aas'] = ','.join(sorted(alternatives-reachable))
                    if alternatives-reachable:
                        t['detail'] = 'Candidate amino acids outside single-nucleotide scope: '+t['unqueryable_alt_aas']
                    for v in variants:
                        v['orientation'] = 'direct_and_reverse_candidate' if v['alt_aa'] in direct&reverse else ('direct_candidate' if v['alt_aa'] in direct else 'reverse_candidate')
                    t['variants'] = variants
                    t['status'] = 'ready' if variants else 'excluded_requires_multiple_nucleotides'
        targets.append(t)
    return targets


def output_row(target, source, status=None, **extra):
    result = {k:target.get(k,'') for k in BASE}
    result.update(source=source,status=status or target['status'],candidate_count=len(target['variants']),annotated_count=0)
    result.update(extra)
    return result


def source_manifest(database, manifest_path, source):
    m = read_json(manifest_path)
    if m.get('schema')!='phylophere.vep_database.v1' or m.get('source')!=source:
        raise ValueError('Invalid annotation database manifest')
    if m.get('species')!='homo_sapiens' or m.get('assembly')!='GRCh38' or m.get('coordinates')!='genomic_1based_forward':
        raise ValueError('Incompatible annotation database species/assembly/coordinates')
    if not m.get('version') or not m.get('source_uri') or m.get('sha256')!=digest(database):
        raise ValueError('Annotation database version/provenance/checksum missing or mismatched')
    integer(m.get('annotation_release'))
    m['_sha256'] = digest(manifest_path)
    return m


def finish(output, header, rows, source, inputs, metadata=None):
    write_table(output,header,rows)
    statuses = {}
    # Count unique objectives, not transcript/variant rows.
    for r in rows:
        statuses.setdefault(r['status'],set()).add((r['Gene'],r['Position'],r['side']))
    write_json(str(output)+'.manifest.json',dict(schema='phylophere.vep_output.v1',source=source,
               complete=True,inputs=inputs,metadata=metadata or {},output_sha256=digest(output),
               implementation_sha256={p.name:digest(p) for p in sorted(Path(__file__).parent.glob('*.py'))},
               targets_by_status={k:len(v) for k,v in statuses.items()}))


def pack_targets(targets, inputs):
    """Store each protein/codon reference once, avoiding positions x protein-size JSON."""
    references = {t['Gene']:t['reference'] for t in targets if t.get('reference')}
    compact = [{k:v for k,v in t.items() if k!='reference'} | {'reference_id':t['Gene'] if t.get('reference') else None} for t in targets]
    return dict(schema='phylophere.vep_targets.v1',targets=compact,references=references,inputs=inputs)


def read_targets(path):
    data=read_json(path)
    if data.get('schema')!='phylophere.vep_targets.v1' or not isinstance(data.get('references'),dict):
        raise ValueError('Invalid prepared VEP targets')
    for t in data['targets']:
        rid=t.pop('reference_id')
        if rid is not None and rid not in data['references']:
            raise ValueError('Unknown reference in prepared VEP targets')
        t['reference']=data['references'].get(rid)
    return data
