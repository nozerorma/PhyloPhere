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
from functools import lru_cache
import warnings

AA = set('ACDEFGHIKLMNPQRSTVWY')
MISSING = {'', 'NA', '-'}
BASE = ['Gene', 'Position', 'side', 'status', 'coverage', 'detail', 'source',
        'species', 'assembly', 'reference_gene_id', 'reference_transcript_id',
        'reference_protein_id', 'reference_aa_position', 'hg38_ref_aa',
        'caas_alt_aas', 'caas_change', 'transfer', 'orientation', 'evidence_basis',
        'variant_id', 'candidate_count', 'annotated_count', 'unqueryable_alt_aas',
        'reference_status', 'mapping_basis', 'candidate_count_known']
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


@lru_cache(maxsize=8)
def map_index(directory):
    exact, prefixed = {}, {}
    for p in sorted(Path(directory).iterdir()):
        if not p.is_file() or not p.name.endswith('.map.tsv'):
            continue
        stem = p.name.removesuffix('.map.tsv')
        exact[stem] = p
        for i,c in enumerate(stem):
            if c == '.':
                prefixed.setdefault(stem[:i], []).append(p)
    return exact, prefixed


def map_path(gene, directory):
    exact, prefixed = map_index(str(Path(directory).resolve()))
    if gene in exact:
        return exact[gene]
    paths = prefixed.get(gene, [])
    if len(paths)>1:
        raise ValueError(f'Multiple MAP references for {gene}: {[p.name for p in paths]}')
    return paths[0] if paths else None


def genomic_coord(text):
    if text in MISSING:
        return None
    match = re.fullmatch(r'([^:\s]+):([1-9][0-9]*)', text)
    if not match:
        raise ValueError(f'Invalid MAP genomic coordinate: {text!r}')
    chrom=match[1]
    if not chrom.startswith('chr'):chrom='chr'+chrom
    if chrom=='chrMT':chrom='chrM'
    return f'{chrom}:{match[2]}'


def load_map(path):
    """Read declared MAP coordinates; do not require reference sidecars.

    A codon is contiguous only when the next reference residue starts exactly
    three genomic bases away in consistent coding order. Otherwise only its
    declared first base is usable. Extended producer columns supply split codons.
    """
    require(path, ['status','prot_ali_col','hg38_aa_pos','hg38_nt_coord'])
    result, all_residues = {}, {}
    for r in table(path):
        if r['status'] not in {'selected','removed'}:
            raise ValueError(f'Unknown MAP status in {path}')
        aa = None if r['hg38_aa_pos'] in MISSING else integer(r['hg38_aa_pos'])
        anchor = genomic_coord(r['hg38_nt_coord'])
        explicit = [anchor, genomic_coord(r.get('hg38_nt_coord_2','')), genomic_coord(r.get('hg38_nt_coord_3',''))]
        if any(explicit[1:]) and not all(explicit):
            raise ValueError('MAP codon must provide all three coordinates or only its anchor')
        strand = r.get('strand','')
        if strand not in {'','NA','+','-'}:raise ValueError('Invalid MAP strand')
        bases = r.get('hg38_codon','').upper()
        if bases not in MISSING and (len(bases)!=3 or set(bases)-set('ACGTN')):
            raise ValueError('Invalid MAP coding triplet')
        entry = dict(aa_position=aa,anchor=anchor,positions=explicit if all(explicit) else ([anchor] if anchor else []),
                     strand=strand if strand in {'+','-'} else None, coding_bases=bases if bases not in MISSING else '',
                     protein_id=r.get('human_protein_id',''))
        if aa is not None:
            if aa in all_residues:raise ValueError('Repeated MAP reference residue')
            all_residues[aa]=entry
        if r['status'] != 'selected':
            if r['prot_ali_col'] not in MISSING:raise ValueError('Removed MAP row carries a trimmed column')
            continue
        col=integer(r['prot_ali_col'])-1
        if col in result:raise ValueError('Repeated MAP alignment column')
        result[col]=entry
    if sorted(result) != list(range(len(result))):raise ValueError('Selected MAP columns must be contiguous 1..n')
    anchors=[e['anchor'] for _,e in sorted(all_residues.items()) if e['anchor']]
    chroms={x.split(':')[0] for x in anchors}
    if len(chroms)>1:raise ValueError('A MAP reference spans different chromosomes')
    steps={1 if int(b.split(':')[1])>int(a.split(':')[1]) else -1 for a,b in zip(anchors,anchors[1:]) if a!=b}
    trend=('+' if steps=={1} else '-' if steps=={-1} else None)
    for aa,e in all_residues.items():
        strand=e['strand'] or trend
        if e['strand'] and trend and strand!=trend:raise ValueError('MAP strand disagrees with coordinate order')
        e['strand']=strand
        if len(e['positions'])==3:
            nums=[int(x.split(':')[1]) for x in e['positions']]
            if len({x.split(':')[0] for x in e['positions']})!=1 or len(set(nums))!=3:
                raise ValueError('Invalid explicit MAP codon')
            order='+' if nums[0]<nums[1]<nums[2] else '-' if nums[0]>nums[1]>nums[2] else None
            if order is None or (strand and strand!=order):raise ValueError('MAP codon order disagrees with strand')
            e['strand']=order
        elif e['anchor'] and strand:
            nxt=all_residues.get(aa+1,{})
            first=int(e['anchor'].split(':')[1]); step=1 if strand=='+' else -1
            if nxt.get('anchor')==e['anchor'].split(':')[0]+':'+str(first+3*step):
                e['positions']=[e['anchor'].split(':')[0]+':'+str(first+j*step) for j in range(3)]
    return result


def identity(value, prefix):
    if not isinstance(value, str) or not re.fullmatch(prefix+r'[0-9]+\.[1-9][0-9]*',value):
        raise ValueError(f'Versioned {prefix} identity required: {value!r}')
    return value


def same_id(actual, expected):
    """Stable IDs may omit a version; supplied versions must agree."""
    if not actual or not expected:return False
    if actual==expected:return True
    a=actual.rsplit('.',1);b=expected.rsplit('.',1)
    return a[0]==b[0] and (len(a)==1 or len(b)==1)


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


def gene_aliases(gene):
    """Aliases for producer species suffixes, retaining literal IDs and versions."""
    out = {gene}
    stem = gene.removesuffix('.prot')
    out.add(stem)
    match = re.fullmatch(r'(.+)\.[A-Z][a-z]+_[a-z]+(?:_[A-Za-z0-9]+)*', stem)
    if match:
        core = match.group(1)
        out.add(core)
        # A numbered suffix is an alias only in the documented species-labelled form.
        out.add(re.sub(r'\.[0-9]+$', '', core))
    return out


_gene_index_cache = None


def _gene_alias_index(identifiers):
    global _gene_index_cache
    if _gene_index_cache is None or _gene_index_cache[0] is not identifiers:
        index = {}
        for name in identifiers:
            for alias in gene_aliases(name):
                index.setdefault(alias, set()).add(name)
        _gene_index_cache = (identifiers, index)
    return _gene_index_cache[1]


def resolve_gene(gene,identifiers):
    if gene in identifiers:return gene
    index = _gene_alias_index(identifiers)
    matches = sorted(set().union(*(index.get(alias, set()) for alias in gene_aliases(gene))))
    def version(name):
        stem = name.removesuffix('.prot')
        stem = re.sub(r'\.[A-Z][a-z]+_[a-z]+(?:_[A-Za-z0-9]+)*$', '', stem)
        match = re.search(r'\.([0-9]+)$', stem)
        return match.group(1) if match else None
    query_version = version(gene)
    matches = [name for name in matches if not (query_version and version(name) and query_version != version(name))]
    if len(matches)>1:raise ValueError(f'Ambiguous producer gene alias: {gene}: {matches}')
    return matches[0] if matches else gene


def protein_mapping(path):
    if not path:return {}
    require(path,['gene','human_protein_id'])
    mapping={}
    for row in table(path):
        if not row['gene'] or row['gene'] in mapping:raise ValueError('Empty or duplicate gene in protein mapping')
        mapping[row['gene']]=row
    return mapping


@lru_cache(maxsize=8)
def alignment_index(directory):
    out={}
    for p in Path(directory).iterdir():
        if p.is_file() and p.suffix.lower() in {'.fa','.fasta','.faa','.fas'}:
            stem=p.stem
            out[stem]=p
    return out


def alignment_reference(directory,gene):
    if not directory:return None
    index=alignment_index(str(Path(directory).resolve()))
    paths=[p for name,p in index.items() if name==gene or name.startswith(gene+'.')]
    if len(paths)>1:
        if gene in index:paths=[index[gene]]
        else:raise ValueError(f'Multiple productive alignments for {gene}')
    if not paths:return None
    refs=[];header='';seq=[]
    def keep():
        if 'homo_sapiens' in header.lower() or 'hg38' in header.lower():refs.append(''.join(seq).upper())
    for line in paths[0].read_text().splitlines():
        if line.startswith('>'):
            keep();header=line[1:];seq=[]
        else:seq.append(line.strip())
    keep()
    if len(refs)>1:raise ValueError(f'Multiple human references in alignment for {gene}')
    return refs[0] if refs else None


def prepare(scores,map_dir,proteins=None,alignments=None):
    rows=load_scores(scores); proteins_map=protein_mapping(proteins)
    cache,targets={},[]
    for row in rows:
        gene,position,side=row['Gene'],row['Position'],row['side']
        t=dict(Gene=gene,Position=position,side=side,status='',coverage='not_evaluated',detail='',
               species='homo_sapiens',assembly='GRCh38',caas_change=row['caas'],
               transfer='human_reference_candidate',evidence_basis='pooled_ASR_state_candidate',
               variants=[],reference=None,candidate_count_known=False)
        if gene not in cache:
            path=map_path(gene,map_dir)
            mapping=load_map(path) if path else {}
            ref=load_reference(path,gene,mapping) if path else None
            sequence=alignment_reference(alignments,gene)
            if sequence is not None and len(sequence)!=len(mapping):raise ValueError(f'MAP/productive alignment length differs: {gene}')
            cache[gene]=(mapping,ref,sequence)
        mapping,ref,sequence=cache[gene]
        if mapping and position not in mapping:raise ValueError(f'SCORING position outside selected MAP: {gene}:{position}')
        anc,der=anc_der_from_row(row['ancestral_aa'],row['derived_aa'])
        t.update(ancestral_aas=sorted(anc),derived_aas=sorted(der))
        if not anc or not der:t['status']='unresolved_asr'
        elif (anc|der)-AA:t['status']='excluded_gap' if (anc|der)&set('-?') else 'excluded_ambiguous_residue'
        elif not mapping:t['status']='unmapped_no_map'
        elif mapping[position]['aa_position'] is None:t['status']='unmapped_reference_gap'
        else:
            e=mapping[position];aa=e['aa_position'];protein=proteins_map.get(resolve_gene(gene,proteins_map),{})
            pid=e.get('protein_id') if e.get('protein_id') not in MISSING|{None} else protein.get('human_protein_id','')
            if pid in MISSING:pid=''
            t.update(reference_aa_position=aa,reference_protein_id=pid,
                     reference_gene_id=protein.get('ensembl_gene_id',''),reference_transcript_id=protein.get('human_transcript_id',''),
                     genomic_positions=e['positions'],strand=e['strand'],map_file=str(map_path(gene,map_dir)),
                     reference_status='map_declared',mapping_basis='MAP protein position and genomic coordinates')
            wt=None
            if ref:
                wt=ref['protein_sequence'][aa-1]
                t.update(reference=ref,reference_gene_id=ref['gene_id'],reference_transcript_id=ref['transcript_id'],
                         reference_protein_id=ref['protein_id'],reference_status='reference_sequence_checked')
                if pid and not same_id(pid,ref['protein_id']):raise ValueError(f'Protein mapping disagrees with reference: {gene}')
                c=ref['codons'].get(str(aa))
                if c:t['genomic_positions']=c['positions'];t['strand']=ref['strand']
            elif sequence:
                if position>=len(sequence):raise ValueError('Scored position outside productive human alignment')
                wt=sequence[position]
            elif e.get('coding_bases') in CODE and not any(p.startswith('chrM:') for p in e['positions']):wt=CODE[e['coding_bases']]
            nuclear = not any(p.startswith('chrM:') for p in e['positions'])
            if wt is not None and wt not in AA:
                t['status']='unmapped_reference_gap' if wt in '-?' else 'excluded_ambiguous_residue'
            elif wt is not None and wt not in anc|der:t['status']='reference_state_mismatch'
            elif nuclear and wt is not None and e.get('coding_bases') in CODE and CODE[e['coding_bases']]!=wt:
                t['status']='reference_state_mismatch'
            else:
                t['status']='ready'
                if wt:
                    t['hg38_ref_aa']=wt
                    alts=(der-{wt} if wt in anc else set()) | (anc-{wt} if wt in der else set())
                    t['caas_alt_aas']=','.join(sorted(alts))
                    if not alts:t['status']='excluded_no_change'
                    if ref and t['status']=='ready':t['variants']=variants_for(ref,aa,alts);t['candidate_count_known']=True
                    elif nuclear and len(e['positions'])==3 and e.get('coding_bases') in CODE and e['strand']:
                        small=dict(protein_sequence='A'*(aa-1)+wt,strand=e['strand'],codons={str(aa):dict(positions=e['positions'],coding_bases=e['coding_bases'])})
                        t['variants']=variants_for(small,aa,alts);t['candidate_count_known']=True
                    if t['candidate_count_known']:
                        unreachable=alts-{v['alt_aa'] for v in t['variants']}
                        t['unqueryable_alt_aas']=','.join(sorted(unreachable))
                        if not t['variants'] and alts:t['status']='excluded_requires_multiple_nucleotides'
        for variant in t['variants']:
            variant['orientation']=candidate_orientation(t,variant['ref_aa'],variant['alt_aa'])
        targets.append(t)
    return targets


def candidate_orientation(target,wt,alt):
    anc=set(target['ancestral_aas']);der=set(target['derived_aas'])
    direct=wt in anc and alt in der and wt!=alt
    reverse=wt in der and alt in anc and wt!=alt
    if direct and reverse:return 'direct_and_reverse_candidate'
    return 'direct_candidate' if direct else 'reverse_candidate' if reverse else None


def output_row(target, source, status=None, **extra):
    result = {k:target.get(k,'') for k in BASE}
    result.update(source=source,status=status or target['status'],candidate_count=len(target['variants']),annotated_count=0)
    result.update(extra)
    return result


def finish(output,header,rows,source,inputs,metadata=None):
    write_table(output,header,rows)
    statuses={}
    for row in rows:statuses.setdefault(row['status'],set()).add((row['Gene'],row['Position'],row['side']))
    write_json(str(output)+'.manifest.json',dict(schema='phylophere.vep_output.v2',source=source,
        complete=True,inputs=inputs,metadata=metadata or {},
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
