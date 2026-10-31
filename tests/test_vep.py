"""Scientific VEP transfer, reference identities, literal keys and snapshot contracts."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'subworkflows/VEP/local/src'
sys.path.insert(0,str(SRC))
import vep_common as vc
from map_database import map_database
from build_vep_variants import build
from join_vep_output import join
from run_ensembl import resource


def js(path,data):
    path.write_text(json.dumps(data));return path


def fixture(path,gene='001',strand='+',rows=None):
    path.mkdir(parents=True,exist_ok=True)
    maps=path/'maps';maps.mkdir(exist_ok=True)
    scores=path/'scores.tsv'
    if rows is None:rows=[(gene,0,'top',.8,'T/A','A','T'),(gene,0,'bottom',.7,'T/A','T','A')]
    scores.write_text('Gene\tPosition\tside\tCAAS_score\tcaas\tancestral_aa\tderived_aa\n'+''.join('\t'.join(map(str,r))+'\n' for r in rows))
    coords=['chr1:100','chr1:101','chr1:900'] if strand=='+' else ['chr1:900','chr1:101','chr1:100']
    mp=maps/(gene+'.map.tsv')
    mp.write_text('status\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\nselected\t1\t1\t'+coords[0]+'\n')
    ref=dict(schema='phylophere.vep_reference.v1',gene=gene,species='homo_sapiens',assembly='GRCh38',
             annotation_release=115,source_uri='fixture://reference/115',map_sha256=vc.digest(mp),
             gene_id='ENSG1.1',transcript_id='ENST1.2',protein_id='ENSP1.2',protein_sequence='A',
             protein_sha256=hashlib.sha256(b'A').hexdigest(),strand=strand,
             codons={'1':dict(positions=coords,coding_bases='GCT')})
    js(Path(str(mp)+'.vep.json'),ref)
    proteins=path/'proteins.tsv';proteins.write_text('gene\thuman_protein_id\n'+gene+'\tENSP1.2\n')
    return scores,maps,proteins,mp,ref


def prepared(path,**kwargs):
    scores,maps,proteins,mp,ref=fixture(path,**kwargs)
    targets=path/'targets.json'
    js(targets,vc.pack_targets(vc.prepare(scores,maps,proteins),{'scores':vc.digest(scores)}))
    return targets,mp,ref


def database(path,source='primateai',rows=None):
    header=['chr','pos','non_flipped_ref','non_flipped_alt','ref_aa','alt_aa','change_position_1based',
            'gene_name','ensembl_gene_id','score_PAI3D','percentile_PAI3D','prediction'] if source=='primateai' else [
            'CHROMOSOME','GENOME_START','GENOMIC_WT_ALLELE','GENOMIC_MUT_ALLELE','MUTATION_AA','TRANSCRIPT_ACCESSION','GENE_ID']
    if rows is None:
        rows=[['chr1','100','G','A','A','T','1','ENST1.2','ENSG1','0.8','0.9','pathogenic']] if source=='primateai' else [
            ['1','100','G','A','p.A1T','ENST1.2','ENSG1']]
    db=path/(source+'.tsv')
    vc.write_table(db,header,[dict(zip(header,r)) for r in rows])
    manifest=dict(schema='phylophere.vep_database.v1',source=source,species='homo_sapiens',assembly='GRCh38',
                  coordinates='genomic_1based_forward',version='fixture1',annotation_release=115,
                  source_uri='fixture://'+source,sha256=vc.digest(db),scope_gene_ids=['ENSG1'])
    meta=js(Path(str(db)+'.vep.json'),manifest)
    return db,meta


def mapped(path,source='primateai',rows=None):
    targets,_,_=prepared(path)
    db,meta=database(path,source,rows)
    out=path/'out.tsv'
    map_database(targets,db,meta,out,source)
    return list(vc.table(out)),json.loads(Path(str(out)+'.manifest.json').read_text())


def test_side_and_orientation_preserved(tmp_path):
    rows,manifest=mapped(tmp_path)
    assert [(r['Gene'],r['Position'],r['side'],r['orientation']) for r in rows]==[
        ('001','0','top','direct_candidate'),('001','0','bottom','reverse_candidate')]
    assert all(r['status']=='annotated' and r['coverage']=='complete' for r in rows)
    assert manifest['targets_by_status']=={'annotated':2}


@pytest.mark.parametrize('gene',['NA','001','G.internal.v2','"Q"',"a 'quote'",'G[1]'])
def test_literal_ids(tmp_path,gene):
    scores,maps,proteins,_,_=fixture(tmp_path,gene)
    targets=vc.prepare(scores,maps,proteins)
    assert [t['Gene'] for t in targets]==[gene,gene]
    out=tmp_path/'literal.tsv';vc.write_table(out,vc.BASE,[vc.output_row(t,'ensembl') for t in targets])
    assert [r['Gene'] for r in vc.table(out)]==[gene,gene]


@pytest.mark.parametrize('change,match',[
    (lambda r:r.update(Position='-1'),'integer'),
    (lambda r:r.update(Position='0.5'),'integer'),
    (lambda r:r.update(side='none'),'side'),
    (lambda r:r.update(CAAS_score='nan'),'finite'),
    (lambda r:r.update(CAAS_score='1.1'),'finite'),
    (lambda r:r.update(CAAS_score='0'),'finite'),
    (lambda r:r.update(ancestral_aa='ALA'),'residue'),
])
def test_scores_reject_invalid(tmp_path,change,match):
    scores,*_=fixture(tmp_path)
    rows=list(vc.table(scores));change(rows[0]);vc.write_table(scores,vc.columns(scores),rows)
    with pytest.raises(ValueError,match=match):vc.load_scores(scores)


def test_duplicate_scores_rejected_before_join(tmp_path):
    scores,*_=fixture(tmp_path)
    with scores.open('a') as f:f.write(scores.read_text().splitlines()[1]+'\n')
    with pytest.raises(ValueError,match='Duplicate'):vc.load_scores(scores)


@pytest.mark.parametrize('anc,der,status',[('','T','unresolved_asr'),('A','','unresolved_asr'),('A','X','excluded_ambiguous_residue'),('A','-','excluded_gap'),('A','W','excluded_requires_multiple_nucleotides'),('V','T','reference_state_mismatch'),('A','A','excluded_no_change')])
def test_no_fabricated_ancestor_or_unmappable_transition(tmp_path,anc,der,status):
    scores,maps,_,_,_=fixture(tmp_path,rows=[('001',0,'top',.8,'T/PS',anc,der)])
    t=vc.prepare(scores,maps)[0]
    assert t['status']==status and not t['variants']
    assert vc.anc_der_from_row('','','T/PS','top')==(set(),set())


@pytest.mark.parametrize('strand',['+','-'])
def test_split_codon_genomic_bases_and_strand(tmp_path,strand):
    targets,_,_=prepared(tmp_path,strand=strand)
    data,variants=build(targets,tmp_path/'out.vcf')
    v=list(variants.values())[0]
    assert v['pos']==(100 if strand=='+' else 900)
    assert (v['ref'],v['alt'])==(('G','A') if strand=='+' else ('C','T'))
    assert len(variants)==1 and v['alt_aa']=='T'
    assert len(data['targets'])==2


@pytest.mark.parametrize('field,value,match',[
    ('assembly','GRCh37','GRCh38'),('species','pan_troglodytes','homo_sapiens'),
    ('strand','','strand'),('map_sha256','wrong','checksum'),('protein_sha256','wrong','checksum'),
    ('protein_id','ENSP1','Versioned'),('transcript_id','ENST1','Versioned')])
def test_reference_identity_and_checksums(tmp_path,field,value,match):
    scores,maps,_,mp,ref=fixture(tmp_path)
    ref[field]=value;js(Path(str(mp)+'.vep.json'),ref)
    with pytest.raises(ValueError,match=match):vc.prepare(scores,maps)


def test_codon_mismatch_rejected(tmp_path):
    scores,maps,_,mp,ref=fixture(tmp_path)
    ref['codons']['1']['coding_bases']='ACT';js(Path(str(mp)+'.vep.json'),ref)
    with pytest.raises(ValueError,match='translation'):vc.prepare(scores,maps)


def test_old_map_and_missing_map_have_explicit_status(tmp_path):
    scores,maps,_,mp,_=fixture(tmp_path)
    Path(str(mp)+'.vep.json').unlink()
    assert {t['status'] for t in vc.prepare(scores,maps)}=={'unverified_reference'}
    mp.unlink()
    assert {t['status'] for t in vc.prepare(scores,maps)}=={'unmapped_no_map'}


def test_duplicate_or_ambiguous_map_rejected(tmp_path):
    scores,maps,_,mp,_=fixture(tmp_path)
    other=maps/'001.isoform.map.tsv';other.write_text(mp.read_text())
    with pytest.raises(ValueError,match='Multiple MAP'):vc.prepare(scores,maps)
    other.unlink()
    with mp.open('a') as f:f.write('selected\t1\t1\tchr1:100\n')
    with pytest.raises(ValueError,match='Repeated MAP'):vc.prepare(scores,maps)


def test_duplicate_and_incompatible_proteins_fail(tmp_path):
    scores,maps,proteins,*_=fixture(tmp_path)
    proteins.write_text('gene\thuman_protein_id\n001\tENSP1.2\n001\tENSP2.1\n')
    with pytest.raises(ValueError,match='duplicate'):vc.prepare(scores,maps,proteins)
    proteins.write_text('gene\thuman_protein_id\n001\tENSP1.3\n')
    with pytest.raises(ValueError,match='disagrees'):vc.prepare(scores,maps,proteins)


@pytest.mark.parametrize('source',['primateai','cosmic'])
def test_empty_snapshot_not_a_biological_negative(tmp_path,source):
    rows,audit=mapped(tmp_path,source,rows=[])
    assert {r['status'] for r in rows}==({'no_hit'} if source=='cosmic' else {'missing_annotation'})
    assert audit['complete']


def test_isoform_mismatch_does_not_match(tmp_path):
    rows,_=mapped(tmp_path,rows=[['chr1','100','G','A','A','T','999','ENST_OTHER','ENSG1','.8','.9','pathogenic']])
    assert {r['status'] for r in rows}=={'missing_annotation'}


@pytest.mark.parametrize('source',['primateai','cosmic'])
def test_matching_transcript_wrong_aa_position_fails(tmp_path,source):
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path,source)
    text=db.read_text().replace('\t1\tENST1','\t999\tENST1') if source=='primateai' else db.read_text().replace('p.A1T','p.A999T')
    db.write_text(text);m=vc.read_json(meta);m['sha256']=vc.digest(db);js(meta,m)
    with pytest.raises(ValueError,match='consequence'):map_database(targets,db,meta,tmp_path/'out.tsv',source)


def test_score_range_and_database_checksum(tmp_path):
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path)
    db.write_text(db.read_text().replace('0.8','1.8'))
    with pytest.raises(ValueError,match='checksum'):map_database(targets,db,meta,tmp_path/'out','primateai')
    m=vc.read_json(meta);m['sha256']=vc.digest(db);js(meta,m)
    with pytest.raises(ValueError,match='score'):map_database(targets,db,meta,tmp_path/'out','primateai')


def test_snapshot_duplicates_coalesce(tmp_path):
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path,'cosmic')
    with db.open('a') as f:f.write(db.read_text().splitlines()[1]+'\n')
    m=vc.read_json(meta);m['sha256']=vc.digest(db);js(meta,m)
    out=tmp_path/'out';map_database(targets,db,meta,out,'cosmic')
    assert len(list(vc.table(out)))==2


def test_census_scope_not_inferred_from_hits(tmp_path):
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path,'cosmic',rows=[])
    m=vc.read_json(meta);m['scope_gene_ids']=['ENSG2'];js(meta,m)
    out=tmp_path/'out';map_database(targets,db,meta,out,'cosmic')
    assert {r['status'] for r in vc.table(out)}=={'outside_source_scope'}


def fake_output(path,vid='chr1:100:G:A',position='1',feature='ENST1.2',extra='CANONICAL=YES;ENSP=ENSP1.2'):
    path.write_text('## version 115\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tProtein_position\tAmino_acids\tExtra\n'+
                    f'{vid}\t1:100\tA\tENSG1\t{feature}\tTranscript\tmissense_variant\t{position}\tA/T\t{extra}\n')


def test_ensembl_join_retains_side_without_duplicating_variants(tmp_path):
    targets,_,_=prepared(tmp_path);data,variants=build(targets,tmp_path/'snvs.vcf')
    raw=tmp_path/'raw';fake_output(raw)
    out=tmp_path/'out';join(raw,data,variants,out,{}, {})
    assert len(list(vc.table(out)))==2
    fake_output(raw,feature='ENST2.1');join(raw,data,variants,out,{}, {})
    assert {r['status'] for r in vc.table(out)}=={'unresolved_annotation'}


@pytest.mark.parametrize('kind',['foreign_id','wrong_position','ragged','no_header'])
def test_ensembl_output_corruption_fails(tmp_path,kind):
    targets,_,_=prepared(tmp_path);data,variants=build(targets,tmp_path/'snvs.vcf')
    raw=tmp_path/'raw';fake_output(raw,vid='foreign' if kind=='foreign_id' else 'chr1:100:G:A',position='999' if kind=='wrong_position' else '1')
    if kind=='ragged':raw.write_text(raw.read_text()+'bad\n')
    if kind=='no_header':raw.write_text('')
    with pytest.raises(ValueError):join(raw,data,variants,tmp_path/'out',{}, {})


def cache_fixture(path):
    path.mkdir(exist_ok=True)
    (path/'homo_sapiens/115_GRCh38').mkdir(parents=True,exist_ok=True)
    (path/'homo_sapiens/115_GRCh38/info.txt').write_text('cache fixture')
    (path/'genome.fa').write_text('>1\nGCT\n');(path/'genome.fa.fai').write_text('index fixture')
    files={str(p.relative_to(path)):vc.digest(p) for p in path.rglob('*') if p.is_file() and p.name!='vep_resource.json'}
    js(path/'vep_resource.json',dict(schema='phylophere.vep_cache.v1',species='homo_sapiens',assembly='GRCh38',release=115,
        source_uri='fixture://cache',fasta='genome.fa',files=files))


def test_offline_resources_check_every_file_and_release(tmp_path):
    cache_fixture(tmp_path)
    assert resource(tmp_path,115)[0]==tmp_path/'genome.fa'
    with pytest.raises(ValueError,match='Incompatible'):resource(tmp_path,114)
    (tmp_path/'genome.fa').write_text('>1\nACT\n')
    with pytest.raises(ValueError,match='checksum'):resource(tmp_path,115)


def test_empty_scores_are_valid_not_missing(tmp_path):
    scores,maps,_,_,_=fixture(tmp_path,rows=[])
    assert vc.prepare(scores,maps)==[]
    scores.write_text('')
    with pytest.raises(ValueError):vc.prepare(scores,maps)


def test_report_literal_keys_sides_and_coverage(tmp_path):
    r=subprocess.run(['Rscript',str(ROOT/'tests/test_vep_report.R'),str(ROOT),str(tmp_path)],capture_output=True,text=True,timeout=60)
    assert r.returncode==0,r.stdout+r.stderr
    assert 'VEP_REPORT_PASS' in r.stdout


def test_full_comparison_render_and_regeneration(tmp_path):
    targets,_,_=prepared(tmp_path,gene='G.internal.2')
    for source,name in [('primateai','pai.tsv'),('cosmic','cosmic_out.tsv')]:
        db,meta=database(tmp_path,source)
        map_database(targets,db,meta,tmp_path/name,source)
    data,variants=build(targets,tmp_path/'variants.vcf')
    raw=tmp_path/'raw';fake_output(raw)
    join(raw,data,variants,tmp_path/'ens.tsv',{}, {})
    scores=tmp_path/'scores.tsv';rows=list(vc.table(scores))
    for r in rows:r.update({'p.emp':'.01','p.adj_bh':'.02','caas_nominal':'TRUE','caas_adjusted':'TRUE'})
    vc.write_table(scores,list(rows[0]),rows)
    before={name:vc.digest(tmp_path/name) for name in ['pai.tsv','cosmic_out.tsv','ens.tsv']}
    r=subprocess.run(['Rscript',str(ROOT/'tests/test_vep_render.R'),str(ROOT),str(tmp_path)],capture_output=True,text=True,timeout=120)
    assert r.returncode==0,r.stdout+r.stderr
    assert 'VEP_RENDER_PASS' in r.stdout
    assert {name:vc.digest(tmp_path/name) for name in before}==before


def test_mixed_snv_and_multinucleotide_candidates_remain_partial(tmp_path):
    targets,_,_=prepared(tmp_path,rows=[('001',0,'top',.8,'TW/A','A','T,W')])
    db,meta=database(tmp_path)
    out=tmp_path/'out';map_database(targets,db,meta,out,'primateai')
    r=next(vc.table(out))
    assert r['status']=='annotated' and r['coverage']=='partial' and r['unqueryable_alt_aas']=='W'


def test_conflicting_snapshot_predictions_rejected(tmp_path):
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path)
    with db.open('a') as f:f.write(db.read_text().splitlines()[1].replace('0.8','0.2').replace('pathogenic','benign')+'\n')
    m=vc.read_json(meta);m['sha256']=vc.digest(db);js(meta,m)
    with pytest.raises(ValueError,match='Conflicting'):map_database(targets,db,meta,tmp_path/'out','primateai')


def test_prepared_reference_registry_avoids_quadratic_repetition(tmp_path):
    targets,_,_=prepared(tmp_path)
    compact=vc.read_json(targets)
    assert list(compact['references'])==['001']
    assert all('reference' not in t for t in compact['targets'])
    loaded=vc.read_targets(targets)['targets']
    assert loaded[0]['reference'] is loaded[1]['reference']


@pytest.mark.parametrize('chrom',['chrM','chr23','chr0'])
def test_unsupported_chromosomes_cannot_use_nuclear_translation(tmp_path,chrom):
    scores,maps,_,mp,ref=fixture(tmp_path)
    mp.write_text(mp.read_text().replace('chr1:',chrom+':'))
    ref['map_sha256']=vc.digest(mp)
    ref['codons']['1']['positions']=[x.replace('chr1:',chrom+':') for x in ref['codons']['1']['positions']]
    js(Path(str(mp)+'.vep.json'),ref)
    with pytest.raises(ValueError,match='coordinate'):vc.prepare(scores,maps)


def test_gzip_snapshot_has_the_same_literal_contract(tmp_path):
    import gzip
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path)
    compressed=tmp_path/'pai.tsv.gz'
    with gzip.open(compressed,'wb') as f:f.write(db.read_bytes())
    m=vc.read_json(meta);m['sha256']=vc.digest(compressed);js(meta,m)
    out=tmp_path/'out';map_database(targets,compressed,meta,out,'primateai')
    assert len(list(vc.table(out)))==2


def test_minus_strand_exact_forward_alleles_match(tmp_path):
    targets,_,_=prepared(tmp_path,strand='-')
    db,meta=database(tmp_path,'cosmic',rows=[['1','900','C','T','p.A1T','ENST1.2','ENSG1']])
    out=tmp_path/'out';map_database(targets,db,meta,out,'cosmic')
    assert {r['status'] for r in vc.table(out)}=={'annotated'}


def test_source_release_cannot_silently_transfer(tmp_path):
    targets,_,_=prepared(tmp_path)
    db,meta=database(tmp_path);m=vc.read_json(meta);m['annotation_release']=114;js(meta,m)
    with pytest.raises(ValueError,match='releases differ'):map_database(targets,db,meta,tmp_path/'out','primateai')
