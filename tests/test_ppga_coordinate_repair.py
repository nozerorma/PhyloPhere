"""End-to-end PPGA repair: trimmed MSA, query residues and original MAP differ."""
import csv
import gzip
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'bin'), str(ROOT/'subworkflows/ENRICHMENT/local/src')]
from ppga_position_references import project_ppga
from build_position_gmt import build_map_cache, collect_pfam_members, collect_ucr_members


def fixture(tmp_path):
    base = tmp_path/'base'
    for name in ['PROT/bmge', 'MAP', 'PROT_VAR_RAW', 'HMMSCAN',
                 'reports/variability', 'reports/positive_selection']:
        (base/name).mkdir(parents=True, exist_ok=True)
    (base/'PROT/bmge/G.fa').write_text('>Homo_sapiens\nAX-S\n>Other\nAT-S\n')
    (base/'MAP/G.map.tsv').write_text(
        'ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
        '1\tremoved\tNA\t8\tchr1:18\n'
        '5\tselected\t1\t9\tchr1:21\n'
        '6\tselected\t2\t10\tchr1:24\n'
        '7\tselected\t3\tNA\tNA\n'
        '8\tselected\t4\t11\tchr1:27\n')
    (base/'PROT_VAR_RAW/G.entropy.tsv').write_text(
        'gene\tposition\tvariability\tg\nG\t1\t0.1\t0\nG\t2\t0.2\t0\n'
        'G\t3\t0.3\t1\nG\t4\t0.4\t0\n')
    (base/'reports/variability/ucr_positions.tsv').write_text(
        'gene\tposition\tregion_type\tmethod\tvariability\n'
        'G\t1\tcore\tsliding\t0.1\nG\t4\tflank_up\tsliding\t0.4\n')
    (base/'reports/variability/domain_variability.tsv').write_text(
        'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tn_positions\tmean_variability\n'
        'G\tPF1\tDomain\tdesc\tCL1\tClan\t2\t3\t999\t999\n')
    v = ['x']*22
    for i, value in {1:'PF1.1',3:'G',5:'3',17:'2',18:'3'}.items():
        v[i] = value
    (base/'HMMSCAN/hmmscan.domtblout').write_text(' '.join(v)+'\n')
    (base/'reports/positive_selection/fubar_sites.tsv').write_text(
        'gene\tsite\tori_codon_col\thg38_aa_pos\thg38_nt_coord\nG\t4\t8\t11\tchr1:27\n')
    return base


def read(path):
    with gzip.open(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def test_repair_preserves_sites_instead_of_double_projecting(tmp_path):
    base = fixture(tmp_path)
    source = (base/'reports/variability/ucr_positions.tsv').read_bytes()
    out = tmp_path/'repaired'
    manifest = project_ppga(base, out)
    ucr = read(out/'ucr_positions.tsv.gz')
    assert [r['position'] for r in ucr] == ['0', '3']
    assert [r['source_position'] for r in ucr] == ['1', '4']
    domain = read(out/'domain_variability.tsv.gz')[0]
    # X consumes a query position; the alignment gap does not. No interval fill.
    assert domain['alignment_positions'] == '1,3'
    assert domain['n_positions'] == '2'
    assert float(domain['mean_variability']) == pytest.approx(.3)
    cache = build_map_cache(str(base/'MAP'), {'G'})
    assert collect_ucr_members(out/'ucr_positions.tsv.gz', cache, {'G'})[0] == {'G': {0}}
    terms = collect_pfam_members(out/'domain_variability.tsv.gz', cache, {'G'})[0]
    assert terms['PF1'][1] == ['G:1', 'G:3']
    assert (base/'reports/variability/ucr_positions.tsv').read_bytes() == source
    assert manifest['status'] == 'complete'
    meta = json.loads(Path(str(out/'ucr_positions.tsv.gz')+'.manifest.json').read_text())
    assert meta['producer_coordinate_contract'] == 'productive_alignment_1based'
    assert meta['alignments'] and meta['maps']
    with pytest.raises(ValueError, match='already exists'):
        project_ppga(base, out)


@pytest.mark.parametrize('bad', ['query_length', 'endpoint', 'map', 'entropy', 'ucr_bounds'])
def test_incompatible_inputs_never_produce_complete_bundle(tmp_path, bad):
    base = fixture(tmp_path)
    if bad == 'query_length':
        p = base/'HMMSCAN/hmmscan.domtblout';v=p.read_text().split();v[5]='4';p.write_text(' '.join(v)+'\n')
    elif bad == 'endpoint':
        p=base/'reports/variability/domain_variability.tsv';p.write_text(p.read_text().replace('\t2\t3\t999', '\t2\t4\t999'))
    elif bad == 'map':
        p=base/'MAP/G.map.tsv';p.write_text(p.read_text().replace('6\tselected\t2', '6\tselected\t1'))
    elif bad == 'entropy':
        p=base/'PROT_VAR_RAW/G.entropy.tsv';p.write_text(p.read_text().replace('G\t4\t0.4\t0\n', ''))
    else:
        p=base/'reports/variability/ucr_positions.tsv';p.write_text(p.read_text().replace('G\t4\tflank', 'G\t5\tflank'))
    out=tmp_path/'bad'
    with pytest.raises(ValueError):
        project_ppga(base, out)
    assert not (out/'reference_bundle_manifest.json').exists()


def test_ambiguous_ucr_is_rejected_before_projection(tmp_path):
    base=fixture(tmp_path)
    with pytest.raises(ValueError, match='has no coordinate_space'):
        collect_ucr_members(base/'reports/variability/ucr_positions.tsv', {}, {'G'})
