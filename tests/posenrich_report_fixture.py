"""Build real POSENRICH component tables for rendering tests, in isolated scratch."""
import csv
import sys
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from posenrich_fixtures import canonical, reference, annotation, calibrate_scores

root=Path(sys.argv[1]);p=Path(sys.argv[2]);p.mkdir(parents=True,exist_ok=True)
src=root/'subworkflows/POSENRICH/local/src'
(p/'maps').mkdir(exist_ok=True)
for g in ['G1','G2']:
    path=p/f'maps/{g}.map.tsv'
    path.write_text('ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'+''.join(f'{i}\tselected\t{i}\t{i}\tchr1:{1000*i}\n' for i in range(1,7)))
    reference(path,g,'A'*6)
(p/'bg.txt').write_text('G1\nG2\n')
(p/'background.output').write_text('G1\t0,1,2,3,4,5\nG2\t0,1,2,3,4,5\n')
obs=pd.DataFrame([(g,i,'bottom' if g=='G2' and i==5 else 'top', .9-.1*i) for g in ('G1','G2') for i in range(6)],columns=['Gene','Position','side','asr_path_score'])
obs.to_csv(p/'observed.tsv',sep='\t',index=False)
scores=obs.rename(columns={'asr_path_score':'CAAS_score'}).assign(**{'p.emp':.2,'p.adj_bh':.5,'nd':2,'null_unseen':'FALSE','caas_nominal':'FALSE','caas_adjusted':'FALSE'})
scores.to_csv(p/'position_scores.tsv',sep='\t',index=False)
(p/'custom.tsv').write_text('G1\t0\tTERM_A\tlocalized\nG1\t1\tTERM_A\tlocalized\nG2\t3\tTERM_A\tlocalized\nG2\t4\tTERM_A\tlocalized\n'+''.join(f'G1\t{i}\tTERM_FULL\twhole protein\n' for i in range(6)))
annotation(p/'custom.tsv',['G1','G2'],tested=p/'background.output',background=p/'bg.txt')
rng=np.random.default_rng(1)
block=canonical(p/'prior',p/'observed.tsv',['G1','G2'],[(g,f'b_{c}',i,float(rng.uniform(.05,1)),np.nan) for c in range(1,41) for g in ['G1','G2'] for i in range(6) if rng.random()<.4],[f'b_{c}' for c in range(1,41)],tested=p/'background.output')
calibrate_scores(p/'position_scores.tsv',p/'prior',p/'background.output',{'G1','G2'})
def run(script,*args):subprocess.run([sys.executable,str(src/script),*map(str,args)],check=True)
run('posenrich_universe.py','--background-genes',p/'bg.txt','--background-positions',p/'background.output','--null-positions',p/'prior/null_positions.parquet','--observed',p/'observed.tsv','--observed-contract',p/'observed.tsv.scoring.json','--counts',p/'prior/counts.json','--roster',p/'prior/cycle_roster.txt','--blocks',block,'--output-dir',p/'universe')
u=p/'universe/universe.tsv.gz';m=p/'universe/universe_manifest.json'
run('posenrich_registry.py','--map_dir',p/'maps','--cleaned_background',p/'bg.txt','--position_scores',p/'position_scores.tsv','--postproc',p/'observed.tsv','--universe',u,'--universe-manifest',m,'--custom_marker_file',p/'custom.tsv','--output_dir',p/'registry')
run('posenrich_characterize.py','--registry_dir',p/'registry','--output_dir',p/'characterization')
for direction in ['all','top','bottom']:
    run('accumulation_test.py','--universe',u,'--universe-manifest',m,'--blocks',block,'--change-side','both' if direction=='all' else direction,'--output',p/f'accumulation/accumulation_{direction}_results.csv')
for name,extra in [('within_gene',[]),('within_gene_holdout',['--holdout','10'])]:
    run('posenrich_within_gene.py','--universe',u,'--observed',p/'position_scores.tsv','--registry-dir',p/'registry','--blocks',block,'--output-dir',p/name,'--min-genes','1',*extra)
run('posenrich_proximity.py','--registry_dir',p/'registry','--output_dir',p/'proximity')
