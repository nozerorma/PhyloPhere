"""Positional annotation contracts: trimmed columns, retained sources and report inputs."""
import csv
import importlib.util
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT/'subworkflows/ENRICHMENT/local/src/build_position_gmt.py'


def write(path, content):
    path.write_text(content)
    return str(path)


def test_coordinate_mapping_layers_and_retired_collections(tmp_path):
    maps=tmp_path/'maps';maps.mkdir()
    # First protein residue is 10, one column removed, one human gap. Neither
    # counting selected rows nor equating raw/trimmed columns can pass this fixture.
    write(maps/'G1.map.tsv', 'ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
          '1\tselected\t1\t10\tchr1:1000000\n'
          '2\tremoved\tNA\t11\tchr1:1000002\n'
          '3\tselected\t2\tNA\tNA\n'
          '4\tselected\t3\t12\tchr1:1000005\n')
    bg=write(tmp_path/'bg','G1\n')
    domain=write(tmp_path/'domains','gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tcoordinate_space\n'
          'G1\tPF1\tD1\tdesc\tCL1\tC1\t11\t12\treference_protein_1based\n'
          'G1\tPF2\tD2\tdesc\tNA\tNA\t12\t12\treference_protein_1based\n')
    ucr=write(tmp_path/'ucr','gene\tposition\tregion_type\tmethod\tvariability\nG1\t4\tcore\tsliding\t0.1\nG1\t2\tcore\tsliding\t0.2\n')
    fubar=write(tmp_path/'fubar','gene\thg38_aa_pos\tis_pos_hit\tis_neg_hit\nG1\t10\t1\t0\nG1\t12\t0\t1\n')
    fade=write(tmp_path/'fade.csv','gene,position,max_bf,target_aa\nG1,0,100,S\nOUT,1,100,T\n')
    custom=write(tmp_path/'custom','G1\t0\tPP2A\tvalidated mapping fixture\n')
    old=write(tmp_path/'retired','invalid legacy data is never read\n')
    out=tmp_path/'out'
    result=subprocess.run([sys.executable,str(BUILDER),'--map_dir',str(maps),'--cleaned_background',bg,
        '--domain_variability_file',domain,'--ucr_positions_file',ucr,'--fubar_sites_file',fubar,
        '--fade_sites_top_file',fade,'--custom_marker_file',custom,'--egg_members_file',old,
        '--cosmic_db',old,'--output_dir',str(out)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert 'retired' in result.stderr
    assert {p.name for p in out.glob('*.gmt')} == {'pfam_domains.gmt','pfam_clans.gmt','genomic_locations.gmt','custom_features.gmt'}
    layers={r[0]:r[2:] for r in csv.reader((out/'characterization_layers.tsv').open(),delimiter='\t')}
    assert layers['UCR_core']==['G1:2']
    assert layers['FUBAR_positive']==['G1:0']
    assert layers['FUBAR_purifying']==['G1:2']
    assert layers['FADE_top_sig']==['G1:0']
    rows={int(r['Position']):r for r in csv.DictReader((out/'position_characterization.tsv').open(),delimiter='\t')}
    assert rows[2]['pfam_domain']=='D1;D2' and rows[2]['pfam_clan']=='C1'
    assert rows[0]['fade_top_sig']=='1'
    domains=(out/'pfam_domains.gmt').read_text()
    assert 'G1:1' not in domains # reference gap is not part of either Pfam hit
    bins=(out/'genomic_locations.gmt').read_text()
    assert 'chr1_0M' in bins and 'chr1_1M' in bins


def test_all_optional_inputs_absent(tmp_path):
    result=subprocess.run([sys.executable,str(BUILDER),'--output_dir',str(tmp_path)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert (tmp_path/'position_characterization.tsv').read_text().count('\n')==1


def test_report_regeneration_uses_sites_and_annotations(tmp_path):
    from gui.generation.report_registry import Listing, _compare_find_slots, _posenrich_find_slots
    files=['selection/fade/top/fade_sites_top.csv','selection/fade/bottom/fade_sites_bottom.csv',
           'selection/fade/top/fade_site_bf_top.tsv','vep/ensembl_vep_mapped.tsv',
           'posenrich/gmts/position_characterization.tsv']
    for name in files:
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.touch()
    for find in (_compare_find_slots,_posenrich_find_slots):
        slots={s.key:s.path for s in find(Listing(str(tmp_path)))}
        assert str(slots['fade_sites_top']).endswith('fade_sites_top.csv')
        assert str(slots['vep_ensembl']).endswith('ensembl_vep_mapped.tsv')
    assert slots['vep_ensembl'] is not None
    slots={s.key:s.path for s in _compare_find_slots(Listing(str(tmp_path)))}
    assert slots['position_characterization'] is not None


def test_pp2a_parser_keeps_ambiguity():
    spec=importlib.util.spec_from_file_location('pp2a',ROOT/'bin/prepare_pp2a.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    sites,problems=mod.parse_modifications('Q13620 2xPhospho [T46(100); S53(99.4)]')
    assert [(s[1],s[2],s[3]) for s in sites]==[('T',46,'100'),('S',53,'99.4')]
    assert not problems and all(s[-1]=='reported_site' for s in sites)
    sites,_=mod.parse_modifications('Q09666 1xPhospho [S3054(100)];1xPhospho [S3836(100)]')
    assert len(sites)==2 and all(s[-1]=='ambiguous_mapping' for s in sites)
    assert mod.parse_modifications('Q8NB15 1xPhospho [S/T]')[1]


def test_pfam_producer_preserves_columns_through_gaps_and_unknown_residues(tmp_path, monkeypatch):
    sys.path.insert(0,str(ROOT/'bin'))
    import compute_domain_variability as producer
    align=tmp_path/'align';align.mkdir()
    write(align/'G1.fa','>Homo_sapiens\nAX-S\n')
    clans=tmp_path/'clans';write(clans,'PF00001\tCL1\tClan\tDomain\tDescription\n')
    hmm=tmp_path/'hmm';write(hmm,'test HMM')
    monkeypatch.setattr(producer,'ensure_pfam_cache',lambda _:(str(hmm),str(clans)))
    def hmmscan(argv,**kwargs):
        assert argv[argv.index('--cpu')+1] == '3'
        fields=['x']*22
        for index,value in {1:'PF00001.1',3:'G1',6:'1e-8',11:'1e-10',17:'1',18:'2'}.items():fields[index]=value
        Path(argv[argv.index('--domtblout')+1]).write_text(' '.join(fields)+'\n')
    monkeypatch.setattr(producer.subprocess,'run',hmmscan)
    out=tmp_path/'domains'
    monkeypatch.setattr(sys,'argv',['compute_domain_variability.py','--alignment-dir',str(align),'--output-dir',str(out),'--cpus','3'])
    producer.main()
    rows=list(csv.DictReader((out/'domain_variability.tsv').open(),delimiter='\t'))
    assert rows[0]['alignment_positions']=='0,3' and rows[0]['coordinate_space']=='CAAS_0based'


def test_retired_toggle_migrates_and_custom_sets_render():
    from gui.models.serialization import from_dict, to_dict
    from gui.generation.render import render_batch
    config=from_dict({'modules':{'enrichment':{'posenrich_eggnog':True,'posenrich_custom_file':'/data/positions.tsv'}}})
    # Existing project schema nests module configurations under modules.
    assert 'posenrich_eggnog' not in to_dict(config)['modules']['enrichment']
    assert 'export POSENRICH_CUSTOM_FILE="/data/positions.tsv"' in render_batch(config)


def test_regeneration_prefers_modular_characterization(tmp_path):
    from gui.generation.report_registry import Listing, _compare_find_slots
    for rel in ('posenrich/gmts/position_characterization.tsv',
                'posenrich/characterization/position_characterization.tsv'):
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    slots = {slot.key: slot.path for slot in _compare_find_slots(Listing(str(tmp_path)))}
    assert str(slots['position_characterization']).endswith('posenrich/characterization/position_characterization.tsv')
