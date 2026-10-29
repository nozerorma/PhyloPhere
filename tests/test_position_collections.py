"""Positional annotation contracts: trimmed columns, retained sources and report inputs."""
import csv
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'subworkflows/POSENRICH/local/src'))
import position_sources as ps  # noqa: E402


def write(path, content):
    path.write_text(content)
    return str(path)


def test_sources_are_projected_to_caas_columns(tmp_path):
    maps=tmp_path/'maps';maps.mkdir()
    # First protein residue is 10, one column removed, one human gap. Neither
    # counting selected rows nor equating raw/trimmed columns can pass this fixture.
    write(maps/'G1.map.tsv', 'ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
          '1\tselected\t1\t10\tchr1:1000000\n'
          '2\tremoved\tNA\t11\tchr1:1000002\n'
          '3\tselected\t2\tNA\tNA\n'
          '4\tselected\t3\t12\tchr1:1000005\n')
    genes={'G1'}
    domain=write(tmp_path/'domains','gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tcoordinate_space\n'
          'G1\tPF1\tD1\tdesc\tCL1\tC1\t11\t12\treference_protein_1based\n'
          'G1\tPF2\tD2\tdesc\tNA\tNA\t12\t12\treference_protein_1based\n')
    ucr=write(tmp_path/'ucr','gene\tposition\tregion_type\tmethod\tvariability\tcoordinate_space\nG1\t4\tcore\tsliding\t0.1\traw_alignment_1based\nG1\t2\tcore\tsliding\t0.2\traw_alignment_1based\n')
    fubar=write(tmp_path/'fubar','gene\thg38_aa_pos\tis_pos_hit\tis_neg_hit\nG1\t10\t1\t0\nG1\t12\t0\t1\n')
    fade=write(tmp_path/'fade.csv','gene,position,max_bf,target_aa\nG1,0,100,S\nOUT,1,100,T\n')
    custom=write(tmp_path/'custom','G1\t0\tPP2A\tvalidated mapping fixture\n')
    cache=ps.build_map_cache(str(maps),genes)
    pending=[]
    regions=[]
    domains,clans=ps.collect_pfam_members(domain,cache,genes,pending=pending,regions=regions)
    assert domains['PF1'][1]==['G1:2'] and domains['PF2'][1]==['G1:2']   # the reference gap (column 1) is in neither hit
    assert clans['CL1'][1]==['G1:2'] and 'NA' not in clans
    assert {(r['collection_id'],r['term_id'],r['start'],r['end']) for r in regions}=={
        ('pfam_domains','PF1',2,2),('pfam_domains','PF2',2,2),('pfam_clans','CL1',2,2)}
    core,flank,region,variability=ps.collect_ucr_members(ucr,cache,genes,pending=pending)
    assert core=={'G1':{2}} and not flank and region=={('G1',2):'core'}   # raw column 4 is trimmed column 3; raw column 2 was removed
    positive,purifying,calls=ps.collect_fubar_members(fubar,cache,genes,pending=pending)
    assert positive=={'G1':{0}} and purifying=={'G1':{2}} and calls=={('G1',0):'positive',('G1',2):'purifying'}
    assert ps.collect_fade_sites(fade,cache,genes,pending=pending)=={'G1':{0}}   # OUT is not in the gene list
    windows=ps.collect_genomic_windows(cache)
    assert windows['chr1_0M'][1]==['G1:0'] and windows['chr1_0M'][0].startswith('Genomic window')
    terms,members=ps.collect_custom_members(custom,cache,genes)
    assert terms['PP2A'][1]==['G1:0'] and members[0]['gene']=='G1'


def test_absent_sources_give_empty_structures():
    assert ps.collect_pfam_members(None,{},set())==({},{})
    assert ps.collect_ucr_members(None,{},set())==({},{},{},{})
    assert ps.collect_fubar_members(None,{},set())==({},{},{})
    assert ps.collect_fade_sites(None,{},set())=={}
    assert ps.collect_custom_members(None,{},set())==({},[])
    assert ps.build_map_cache(None,None)=={}


def test_report_regeneration_uses_sites_and_annotations(tmp_path):
    from gui.generation.report_registry import Listing, _compare_find_slots
    files=['selection/fade/top/fade_sites_top.csv','selection/fade/bottom/fade_sites_bottom.csv',
           'selection/fade/top/fade_site_bf_top.tsv','vep/ensembl_vep_mapped.tsv',
           'posenrich/characterization/position_characterization.tsv']
    for name in files:
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.touch()
    for find in (_compare_find_slots,):
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


def test_saved_projects_with_retired_orthogroup_fields_load():
    from gui.models.serialization import from_dict, to_dict
    retired = {'eggnog_taxid': '9443', 'auto_fetch_eggnog': True, 'egg_members_file': '/x/m.tsv', 'egg_annotations_file': '/x/a.tsv'}
    for section in ('enrichment', 'posenrich'):
        config = from_dict({'modules': {section: dict(retired)}})
        dumped = to_dict(config)['modules']
        assert not set(retired) & set(dumped['posenrich']) and not set(retired) & set(dumped['enrichment'])


def test_custom_position_sets_live_in_the_posenrich_section():
    from gui.models.serialization import from_dict, to_dict
    from gui.generation.render import render_batch
    # A project saved with the field in the enrichment section moves it to the POSENRICH section.
    old = from_dict({'modules': {'enrichment': {'posenrich_custom_file': '/data/old.tsv'}}})
    dumped = to_dict(old)['modules']
    assert dumped['posenrich']['posenrich_custom_file'] == '/data/old.tsv'
    assert 'posenrich_custom_file' not in dumped['enrichment']
    # The value in the POSENRICH section wins over a stale one in the enrichment section.
    both = from_dict({'modules': {'posenrich': {'posenrich_custom_file': '/data/new.tsv'},
                                  'enrichment': {'posenrich_custom_file': '/data/old.tsv'}}})
    assert to_dict(both)['modules']['posenrich']['posenrich_custom_file'] == '/data/new.tsv'
    assert 'export POSENRICH_CUSTOM_FILE="/data/new.tsv"' in render_batch(both)


def test_custom_sets_field_is_on_the_posenrich_tab_only():
    import pytest
    pytest.importorskip('PySide6')
    from gui.widgets.tabs import enrichment_tab, posenrich_tab

    def names(spec):
        fields = spec.essential_fields + spec.advanced_fields + spec.fallback_fields
        return [getattr(item, 'name', None) for item in fields]
    assert 'posenrich_custom_file' in names(posenrich_tab.SPEC)
    assert 'posenrich_custom_file' not in names(enrichment_tab.SPEC)
