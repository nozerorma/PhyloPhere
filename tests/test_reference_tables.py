"""Coordinate projections and biological source/peptide checks for real reference bundles."""
import csv,gzip,json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
sys.path.insert(0,str(ROOT/'subworkflows/ENRICHMENT/local/src'))
from build_position_gmt import collect_ucr_members,collect_fubar_members
from map_pp2a_to_caas import project
from build_reference_tables import project_ppga


def test_explicit_sources_work_without_map_and_reject_unknown_space(tmp_path):
 u=tmp_path/'ucr.tsv';u.write_text('gene\tposition\tregion_type\tmethod\tvariability\tcoordinate_space\nG\t0\tcore\tsliding\t0.1\tCAAS_0based\n')
 assert collect_ucr_members(u,{}, {'G'})[0]=={'G':{0}}
 f=tmp_path/'fubar.tsv';f.write_text('gene\thg38_aa_pos\tis_pos_hit\tis_neg_hit\talignment_position\tcoordinate_space\nG\t9\t1\t0\t0\tCAAS_0based\n')
 assert collect_fubar_members(f,{}, {'G'})[0]=={'G':{0}}
 u.write_text(u.read_text().replace('CAAS_0based','invented'))
 with pytest.raises(ValueError,match='coordinate_space'):collect_ucr_members(u,{}, {'G'})


def test_pp2a_requires_source_and_target_peptide_validation(tmp_path):
 prepared=tmp_path/'prepared';prepared.mkdir();out=tmp_path/'out'
 fields=['record_id','study','pmid','cell_line','sheet','source_row','collection_role','condition','accession','amino_acid','protein_position','localization_percent','site_status','coordinate_space']
 rows=[]
 for i,(acc,status,role) in enumerate([('P12345','reported_site','hits'),('P12345-2','reported_site','hits'),('P12345','ambiguous_mapping','hits'),('P12345','reported_site','measured_background')]):
  r=dict.fromkeys(fields,'');r.update(record_id=str(i),study='Stewart',pmid='123',accession=acc,amino_acid='S',protein_position=3,site_status=status,collection_role=role,condition='mutation');rows.append(r)
 with (prepared/'site_candidates.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 with gzip.open(prepared/'source_records.tsv.gz','wt') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['record_id','source_row_json'])
  for r in rows:w.writerow([r['record_id'],json.dumps([['Annotated Sequence','[R].ASDF.[K]']])])
 uni=tmp_path/'uni.fa';uni.write_text('>sp|P12345|TEST_HUMAN GN=G\nMASDFK\n')
 ref=tmp_path/'ref.tsv';ref.write_text('gene\tprotein_position\tposition\tamino_acid\nG\t1\t0\tM\nG\t2\t2\tA\nG\t3\t3\tS\nG\t4\t4\tD\nG\t5\t5\tF\nG\t6\t6\tK\n')
 project(prepared,uni,ref,out)
 mapped=list(csv.DictReader((out/'mapped_sites.tsv').open(),delimiter='\t'))
 assert len(mapped)==2 and all(r['alignment_position']=='3' for r in mapped)
 assert len((out/'pp2a_custom.tsv').read_text().splitlines())==1
 assert len((out/'pp2a_measured_background.tsv').read_text().splitlines())==1
 pending=list(csv.DictReader((out/'pending_sites.tsv').open(),delimiter='\t'))
 assert {r['mapping_status'] for r in pending}=={'accession_or_isoform_not_in_reference_fasta','ambiguous_source_mapping'}
 # A known site alone is insufficient if the peptide crosses missing trimmed residues.
 ref.write_text(ref.read_text().replace('G\t4\t4\tD\n',''))
 project(prepared,uni,ref,tmp_path/'out_missing')
 assert not (tmp_path/'out_missing/pp2a_custom.tsv').read_text()


def test_ppga_projection_detects_fubar_map_disagreement(tmp_path):
 base=tmp_path/'base';(base/'MAP').mkdir(parents=True)
 (base/'reports/variability').mkdir(parents=True);(base/'reports/positive_selection').mkdir()
 (base/'MAP/G.map.tsv').write_text('ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n5\tselected\t1\t9\tchr1:42\n')
 (base/'PROT/bmge').mkdir(parents=True);(base/'PROT/bmge/G.fa').write_text('>Homo_sapiens\nS\n')
 (base/'PROT_VAR_RAW').mkdir();(base/'PROT_VAR_RAW/G.entropy.tsv').write_text('gene\tposition\tvariability\tg\nG\t1\t0.1\t0\n')
 (base/'reports/variability/ucr_positions.tsv').write_text('gene\tmethod\tposition\nG\tsliding\t1\n')
 source=base/'reports/positive_selection/fubar_sites.tsv'
 source.write_text('gene\tsite\tori_codon_col\thg38_aa_pos\thg38_nt_coord\nG\t1\t5\t9\tchr1:42\n')
 project_ppga(base,tmp_path/'out')
 with gzip.open(tmp_path/'out/ucr_positions.tsv.gz','rt') as f:assert next(csv.DictReader(f,delimiter='\t'))['position']=='0'
 source.write_text(source.read_text().replace('\t9\t','\t10\t'))
 with pytest.raises(ValueError,match='hg38_aa_pos mismatch'):project_ppga(base,tmp_path/'bad')


def test_annotation_only_does_not_retain_normalized_members(tmp_path):
 u=tmp_path/'ucr.tsv';u.write_text('gene\tposition\tregion_type\tmethod\tvariability\tcoordinate_space\nG\t0\tcore\tsliding\t0.1\tCAAS_0based\n')
 standard=collect_ucr_members(u,{}, {'G'})
 compact=collect_ucr_members(u,{}, {'G'},include_members=False)
 assert standard[:-1]==compact[:-1] and standard[-1] and compact[-1]==[]


def test_reference_preflight_checks_alignment_hash_and_inventory(tmp_path):
 import shutil,subprocess,hashlib
 nf=shutil.which('nextflow')
 if not nf:pytest.skip('Nextflow unavailable')
 align=tmp_path/'align';align.mkdir();fa=align/'G.fa';fa.write_text('>Homo_sapiens\nAS\n')
 table=tmp_path/'pfam.tsv';table.write_text('gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tcoordinate_space\n')
 meta=dict(alignment_dir=str(align),alignments={'G.fa':hashlib.sha256(fa.read_bytes()).hexdigest()},table_sha256=hashlib.sha256(table.read_bytes()).hexdigest())
 Path(str(table)+'.manifest.json').write_text(json.dumps(meta))
 script=ROOT/'.smoke_reference_preflight.nf'
 script.write_text("include { validatePositionReferences } from './workflows/enrichment.nf'\nworkflow { validatePositionReferences(); log.info 'REFERENCE_PREFLIGHT_PASSED' }\n")
 cmd=[nf,'run',str(script),'--alignment',str(align),'--posenrich_domains','true','--domain_variability_file',str(table),'-ansi-log','false']
 try:
  result=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=90)
  assert result.returncode==0,result.stdout+result.stderr
  fa.write_text('>Homo_sapiens\nAT\n')
  result=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=90)
  assert result.returncode!=0 and 'Productive alignment changed' in result.stdout+result.stderr
  fa.write_text('>Homo_sapiens\nAS\n');(align/'NEW.fa').write_text('>Homo_sapiens\nAS\n')
  result=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=90)
  assert result.returncode!=0 and 'Alignment inventory changed' in result.stdout+result.stderr
 finally:script.unlink(missing_ok=True)


def test_pfam_preserves_sequence_and_domain_cutoffs(tmp_path,monkeypatch):
 import compute_domain_variability as p
 align=tmp_path/'align';align.mkdir();(align/'G.fa').write_text('>Homo_sapiens\nAST\n')
 hmm=tmp_path/'hmm';hmm.write_text('test');clans=tmp_path/'clans';clans.write_text('')
 monkeypatch.setattr(p,'ensure_pfam_cache',lambda _:(str(hmm),str(clans)))
 def scan(argv,**kwargs):
  assert argv[argv.index('-E')+1]=='0.001'
  rows=[]
  for acc,seq_e,domain_e in [('PF00001',1e-5,1e-5),('PF00002',.002,1e-5),('PF00003',1e-5,.02)]:
   fields=['x']*22
   for i,v in {0:'domain',1:acc+'.1',3:'G',6:seq_e,11:domain_e,17:1,18:2}.items():fields[i]=str(v)
   rows.append(' '.join(fields)+' desc')
  Path(argv[argv.index('--domtblout')+1]).write_text('\n'.join(rows)+'\n')
 monkeypatch.setattr(p.subprocess,'run',scan)
 out=tmp_path/'out';monkeypatch.setattr(sys,'argv',['compute_domain_variability.py','--alignment-dir',str(align),'--output-dir',str(out)])
 p.main()
 rows=list(csv.DictReader((out/'domain_variability.tsv').open(),delimiter='\t'))
 assert [r['pfam_id'] for r in rows]==['PF00001']


def test_direct_fubar_keeps_tested_sites_without_human_residue(tmp_path):
 f=tmp_path/'fubar.tsv';f.write_text('gene\thg38_aa_pos\tis_pos_hit\tis_neg_hit\talignment_position\tcoordinate_space\nG\tNA\t1\t0\t4\tCAAS_0based\n')
 result=collect_fubar_members(f,{}, {'G'})
 assert result[0]=={'G':{4}}
 assert result[2][('G',4)]=='positive'
 assert len(result[-1])==2 and 'hg38_aa_pos=NA' in result[-1][0]['detail']
