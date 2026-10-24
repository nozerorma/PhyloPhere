#!/usr/bin/env python3
"""Project prepared PP2A records using accession-verified, unique exact peptides.

Requires the supplied UniProt FASTA (including explicit isoforms) and a productive
reference-position TSV (gene,protein_position,position,amino_acid). Missing, ambiguous,
trimmed or mismatching peptides stay pending. Exports hits and measured background
separately; candidate records remain traceable and are never independent replicates.
"""
import argparse,csv,gzip,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path


def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()


def fasta(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt') as f:
  header=None;buf=[]
  for line in f:
   if line.startswith('>'):
    if header:yield header,''.join(buf)
    header=line[1:].strip();buf=[]
   else:buf.append(line.strip())
  if header:yield header,''.join(buf)


def peptide_from_raw(raw,study):
 pairs=json.loads(raw);values=dict(pairs)
 text=str(values.get('Peptide' if study=='Kauko' else 'Annotated Sequence',''))
 if study=='Kauko': text=text.split(',')[0]
 elif '.' in text:
  parts=text.split('.')
  if len(parts)>=3:text='.'.join(parts[1:-1])
 text=re.sub(r'\[[^]]*\]|\([^)]*\)','',text).upper().strip('_ ')
 return text if re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',text or '') else ''


def occurrences(seq,peptide):
 start=0;found=[]
 while True:
  pos=seq.find(peptide,start)
  if pos<0:return found
  found.append(pos);start=pos+1


def project(prepared,uniprot,reference,output):
 output.mkdir(parents=True,exist_ok=True)
 candidates=list(csv.DictReader((prepared/'site_candidates.tsv').open(),delimiter='\t'))
 needed={r['accession'] for r in candidates};proteins={}
 for header,seq in fasta(uniprot):
  parts=header.split('|');acc=parts[1] if len(parts)>2 else header.split()[0]
  if acc in needed:
   gene=re.search(r'\bGN=(\S+)',header)
   if gene:proteins[acc]=(gene[1],seq)
 genes={p[0] for p in proteins.values()};ref=defaultdict(dict)
 with reference.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['gene'] in genes:ref[r['gene']][int(r['protein_position'])]=(int(r['position']),r['amino_acid'])
 target={gene:''.join(data.get(i,(None,'?'))[1] for i in range(1,max(data)+1)) for gene,data in ref.items()}
 raw={}
 with gzip.open(prepared/'source_records.tsv.gz','rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):raw[r['record_id']]=r['source_row_json']
 fields=list(candidates[0])+['gene','source_peptide','reference_protein_position','alignment_position','mapping_status']
 counts=Counter();custom=set();coverage=set()
 with (output/'mapped_sites.tsv').open('w',newline='') as mf,(output/'pending_sites.tsv').open('w',newline='') as pf:
  mw=csv.DictWriter(mf,fieldnames=fields,delimiter='\t');mw.writeheader()
  pw=csv.DictWriter(pf,fieldnames=fields,delimiter='\t');pw.writeheader()
  for original in candidates:
   r=dict(original,gene='',source_peptide='',reference_protein_position='',alignment_position='',mapping_status='')
   protein=proteins.get(r['accession']);reason=''
   if r['site_status']!='reported_site':reason='ambiguous_source_mapping'
   elif protein is None:reason='accession_or_isoform_not_in_reference_fasta'
   else:
    gene,seq=protein;r['gene']=gene
    peptide=peptide_from_raw(raw[r['record_id']],r['study']);r['source_peptide']=peptide
    pos=int(r['protein_position'])-1
    if not peptide:reason='unparseable_peptide'
    elif pos<0 or pos>=len(seq) or seq[pos] not in 'STY':reason='source_site_not_phosphorylatable_or_outside_sequence'
    elif r['amino_acid'] and seq[pos]!=r['amino_acid']:reason='source_residue_mismatch'
    else:
     source_hits=[start for start in occurrences(seq,peptide) if start<=pos<start+len(peptide)]
     if len(source_hits)!=1:reason='peptide_does_not_uniquely_cover_source_site'
     elif gene not in target:reason='no_productive_human_reference'
     else:
      hits=occurrences(target[gene],peptide)
      if len(hits)!=1:reason='peptide_absent_or_nonunique_in_productive_reference'
      else:
       rp=hits[0]+pos-source_hits[0]+1;col,aa=ref[gene][rp]
       if aa!=seq[pos]:reason='target_residue_mismatch'
       else:r.update(reference_protein_position=rp,alignment_position=col,coordinate_space='CAAS_0based',mapping_status='validated_accession_and_unique_exact_peptide')
   key=r['study']+'/'+r['collection_role']
   if reason:r['mapping_status']=reason;pw.writerow(r);counts[key+'/'+reason]+=1
   else:
    mw.writerow(r);counts[key+'/validated']+=1
    term='PP2A:'+r['study']+':'+r['condition']
    desc=f"{r['study']} PMID {r['pmid']}; perturbation-associated phosphosite; source records retained in mapped_sites.tsv"
    row=(r['gene'],str(r['alignment_position']),term,desc)
    (custom if r['collection_role']=='hits' else coverage).add(row)
 for name,rows in [('pp2a_custom.tsv',custom),('pp2a_measured_background.tsv',coverage)]:
  with (output/name).open('w') as f:
   for row in sorted(rows):f.write('\t'.join(row)+'\n')
 manifest=dict(status='mapped_with_explicit_pending_records',counts=dict(counts),unique_hit_members=len(custom),
               unique_measured_members=len(coverage),coordinate_space='CAAS_0based',
               sources={str(p):sha(p) for p in [prepared/'site_candidates.tsv',prepared/'source_records.tsv.gz',uniprot,reference]},
               validation='Explicit accession/isoform; source peptide covers source residue; exact peptide occurs once in productive reference; ambiguous sources excluded',
               pending_policy='Pending records are unknown, never negative sites',
               coverage='Measured background exported separately; not included among hit terms',
               interpretation='Perturbation-associated sites, not proof of direct PP2A substrates. No new statistical cutoff or inference.')
 (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest['counts'],indent=2))


if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['prepared','uniprot','reference','output']:p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();project(a.prepared,a.uniprot,a.reference,a.output)
