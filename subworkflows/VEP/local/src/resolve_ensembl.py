"""Resolve the declared human protein against Ensembl for legacy MAP targets.

Only used by the enabled Ensembl database source. No database or reference is
silently downloaded for PrimateAI/COSMIC. Responses are kept in the task audit.
"""
import json
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request,urlopen
from vep_common import AA,CODE,candidate_orientation,same_id,variants_for,write_json


class Ensembl:
    def __init__(self):self.cache={}
    def get(self,path):
        if path not in self.cache:
            try:
                with urlopen(Request('https://rest.ensembl.org/'+path,headers={'Accept':'application/json'}),timeout=45) as r:
                    self.cache[path]=json.load(r)
            except HTTPError as exc:
                if exc.code==404:self.cache[path]=None
                else:raise RuntimeError(f'Ensembl reference lookup failed: HTTP {exc.code}') from exc
        return self.cache[path]


def resolve(data,client=None):
    client=client or Ensembl()
    for t in data['targets']:
        if t['status']!='ready' or t['variants']:continue
        if any(p.startswith('chrM:') for p in t['genomic_positions']):
            t.update(status='unmapped_codon',detail='Ensembl SNV enumeration currently supports the nuclear genetic code only');continue
        pid=t.get('reference_protein_id','')
        if not pid:
            t.update(status='unmapped_protein_id',detail='Ensembl annotation needs a declared human protein ID');continue
        enc=quote(pid,safe='')
        protein=client.get(f'sequence/id/{enc}?type=protein')
        identity=client.get(f'lookup/id/{enc}')
        if not protein or not identity:
            t.update(status='unresolved_annotation',detail='Declared protein not found in Ensembl');continue
        seq=protein.get('seq','');aa=t['reference_aa_position']
        if aa>len(seq) or seq[aa-1] not in AA:
            t.update(status='unmapped_reference_gap');continue
        wt=seq[aa-1]
        if t.get('hg38_ref_aa') and t['hg38_ref_aa']!=wt:
            t.update(status='reference_state_mismatch',detail='Productive human alignment and declared Ensembl protein disagree');continue
        tx=identity.get('Parent','')
        if not tx:raise ValueError('Ensembl protein response lacks parent transcript')
        transcript=client.get('lookup/id/'+quote(tx,safe=''))
        cds=client.get('sequence/id/'+quote(tx,safe='')+'?type=cds')
        mapped=client.get(f'map/translation/{enc}/{aa}..{aa}')
        if not cds or not mapped or not transcript:
            t.update(status='unresolved_annotation',detail='Ensembl CDS or genomic mapping unavailable');continue
        coding=cds.get('seq','')[(aa-1)*3:aa*3].upper()
        if CODE.get(coding)!=wt:
            t.update(status='reference_state_mismatch',detail='CDS translation differs from declared protein');continue
        intervals=mapped.get('mappings',[])
        if not intervals:
            t.update(status='unmapped_codon');continue
        if any(x.get('assembly_name')!='GRCh38' for x in intervals):raise ValueError('Ensembl mapping is not GRCh38')
        strands={x['strand'] for x in intervals};chroms={x['seq_region_name'] for x in intervals}
        if len(strands)!=1 or len(chroms)!=1:raise ValueError('Ensembl codon mapping has mixed chromosomes/strands')
        strand='+' if strands.pop()==1 else '-'
        positions=sorted({p for x in intervals for p in range(x['start'],x['end']+1)},reverse=strand=='-')
        if len(positions)!=3:raise ValueError('Ensembl translation mapping does not span three bases')
        chrom=str(chroms.pop());chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
        if chrom=='chrMT':chrom='chrM'
        coords=[f'{chrom}:{p}' for p in positions]
        declared=t['genomic_positions']
        if declared and coords[:len(declared)]!=declared:
            t.update(status='reference_state_mismatch',detail='Ensembl protein uses different genomic coordinates from the MAP');continue
        alts={a for a in t['derived_aas']+t['ancestral_aas'] if candidate_orientation(t,wt,a)}
        ref=dict(protein_sequence=seq,strand=strand,codons={str(aa):dict(positions=coords,coding_bases=coding)})
        t['variants']=variants_for(ref,aa,alts)
        for v in t['variants']:v['orientation']=candidate_orientation(t,wt,v['alt_aa'])
        t.update(hg38_ref_aa=wt,genomic_positions=coords,strand=strand,candidate_count_known=True,
                 reference_transcript_id=transcript['id']+'.'+str(transcript['version']) if transcript.get('version') else tx,
                 reference_gene_id=transcript.get('Parent',''),reference_status='ensembl_sequence_and_MAP_locus_checked')
        t['unqueryable_alt_aas']=','.join(sorted(alts-{v['alt_aa'] for v in t['variants']}))
        if not t['variants']:t['status']='excluded_requires_multiple_nucleotides' if alts else 'excluded_no_change'
    if hasattr(client,'cache'):write_json('ensembl_reference_responses.json',client.cache)
    return data
