"""Paired checkout integration: native alignment producers -> PhyloPhere, no converter."""
import csv
import importlib.util
import os
from pathlib import Path
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
UPSTREAM=Path(os.environ.get('PHYLOPHERE_ALIGNMENT_PIPELINE',ROOT.parent/'ortholog_characterizator'))
SRC=UPSTREAM/'subworkflows/variability/local'
pytestmark=pytest.mark.skipif(not (SRC/'map_domain_variability.py').exists(),reason='paired alignment pipeline checkout unavailable')


def load(name,path):
 sys.path.insert(0,str(path.parent))
 spec=importlib.util.spec_from_file_location(name,path)
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 return m


def test_native_producer_outputs_are_consumed_without_coordinate_conversion(tmp_path,monkeypatch):
 pfam=load('upstream_pfam_contract',SRC/'map_domain_variability.py')
 ucr=load('upstream_ucr_contract',SRC/'aggregate_ucr.py')
 sys.path.insert(0,str(ROOT/'subworkflows/POSENRICH/local/src'))
 import position_sources as consumer
 prot=tmp_path/'prot';prot.mkdir();(prot/'G.fa').write_text('>Homo_sapiens\nAX-S\n')
 ent=tmp_path/'entropy';ent.mkdir();(ent/'G.entropy.tsv').write_text(
 'gene\tposition\tt\tr\tg\tC_trident\tvariability\n'
 'G\t1\t0\t0\t0\t.9\t.1\nG\t2\t0\t0\t0\t.8\t.2\nG\t3\t0\t0\t1\t.7\t.3\nG\t4\t0\t0\t0\t.6\t.4\n')
 scan=tmp_path/'scan';v=['x']*22
 for i,val in {0:'Domain',1:'PF1.1',3:'G',5:'3',11:'1e-8',17:'2',18:'3'}.items():v[i]=val
 scan.write_text(' '.join(v)+' description\n');domain=tmp_path/'domain.tsv'
 monkeypatch.setattr(sys,'argv',['map_domain_variability.py','--prot_raw_dir',str(prot),'--var_dir',str(ent),
 '--hmmscan_domtblout',str(scan),'--out_tsv',str(domain),'--alignment_coordinate_space','trimmed_alignment_1based'])
 pfam.main()
 window=dict(gene='G',ucr_id='u',method='sliding',start_pos=4,end_pos=4,flank_start=4,flank_end=4,
 n_positions=1,mean_C_trident=.6,min_C_trident=.6,mean_variability=.4,mean_gap=0)
 monkeypatch.setattr(ucr,'load_ucr_tsvs',lambda _: {'G':[window]})
 monkeypatch.setattr(ucr,'load_taxonomy',lambda *a: {})
 monkeypatch.setattr(ucr,'compute_ucr_clade_variability',lambda *a,**k: [])
 out=tmp_path/'ucr'
 monkeypatch.setattr(sys,'argv',['aggregate_ucr.py','--ucr_dir','unused','--entropy_dir',str(ent),'--prot_dir',str(prot),
 '--taxid_tsv','unused','--domain_tsv',str(domain),'--out_dir',str(out),'--alignment_coordinate_space','trimmed_alignment_1based'])
 ucr.main()
 maps=tmp_path/'maps';maps.mkdir();(maps/'G.map.tsv').write_text(
 'ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
 '1\tremoved\tNA\t8\tchr1:18\n5\tselected\t1\t10\tchr1:21\n6\tselected\t2\t11\tchr1:24\n'
 '7\tselected\t3\tNA\tNA\n8\tselected\t4\t12\tchr1:27\n')
 cache=consumer.build_map_cache(str(maps),{'G'})
 assert consumer.collect_pfam_members(domain,cache,{'G'})[0]['PF1'][1]==['G:1','G:3']
 assert consumer.collect_ucr_members(out/'ucr_positions.tsv',cache,{'G'})[0]=={'G':{3}}
 with (out/'ucr_positions.tsv').open() as f:r=next(csv.DictReader(f,delimiter='\t'))
 assert r['position']=='4' and r['coordinate_space']=='trimmed_alignment_1based'
 assert r['pfam_domain']=='PF1'
