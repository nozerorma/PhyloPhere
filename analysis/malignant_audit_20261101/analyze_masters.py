"""Descriptive comparison of the 121 archived paired observed masters; no p-values."""
import csv,gzip,json,re
from pathlib import Path
import numpy as np
import pandas as pd

root=Path('.tmp/asr_deep_audit');snapshot=root/'snapshot/.tmp/malignant_audit/new_batch1/b0_observed'
out=Path('analysis/malignant_audit_20261101');old={};new={}
for r in csv.DictReader((root/'inputs/old_master.csv').open()):
 old[(r['gene'],int(r['msa_pos']),r['side'])]=r
oldgenes={k[0] for k in old}
def canonical(gene):
 candidates=[gene,re.sub(r'\.[A-Z][a-z]+_[a-z]+$','',gene)]
 candidates += [re.sub(r'\.\d+$','',g) for g in candidates[:]]
 hits=list(dict.fromkeys(g for g in candidates if g in oldgenes))
 if not hits:raise ValueError('No old gene match: '+gene)
 return hits[0]
for f in snapshot.glob('*.master.csv.gz'):
 with gzip.open(f,'rt') as fh:
  for r in csv.DictReader(fh):new[(canonical(r['gene']),int(r['msa_pos']),r['side'])]=r
genes=sorted({k[0] for k in new});old={k:v for k,v in old.items() if k[0] in genes}
rows=[]
for k in sorted(old.keys()|new.keys()):
 a,b=old.get(k),new.get(k)
 rows.append(dict(gene=k[0],position=k[1],side=k[2],old_score=float(a['asr_path_score']) if a else np.nan,new_score=float(b['asr_path_score']) if b else np.nan,old_type=a['convergence_type'] if a else '',new_type=b['convergence_type'] if b else '',old_present=a is not None,new_present=b is not None))
df=pd.DataFrame(rows);df['delta']=df.new_score-df.old_score;df.to_csv(out/'master_comparison_121.tsv',sep='\t',index=False)
gp=[]
for g in genes:
 a=df[df.gene==g];shared=a[a.old_present&a.new_present];gp.append(dict(gene=g,old_rows=int(a.old_present.sum()),new_rows=int(a.new_present.sum()),changed_scores=int((shared.delta.abs()>1e-9).sum()),changed_types=int((shared.old_type!=shared.new_type).sum()),max_abs_delta=float(shared.delta.abs().max()),old_max=float(a.old_score.max()),new_max=float(a.new_score.max()),old_positive=int((a.old_score>0).sum()),new_positive=int((a.new_score>0).sum())))
gd=pd.DataFrame(gp);gd['max_delta']=gd.new_max-gd.old_max
gd['old_rank']=gd.old_max.rank(method='average',ascending=False);gd['new_rank']=gd.new_max.rank(method='average',ascending=False);gd.to_csv(out/'gene_comparison_121.tsv',sep='\t',index=False)
shared=df[df.old_present&df.new_present];oldpos={k[:2] for k in old};newpos={k[:2] for k in new}
summary=dict(genes=len(genes),old_rows=len(old),new_rows=len(new),shared_rows=len(shared),old_only_rows=len(old.keys()-new.keys()),new_only_rows=len(new.keys()-old.keys()),old_only_positions=len(oldpos-newpos),new_only_positions=len(newpos-oldpos),changed_scores=int((shared.delta.abs()>1e-9).sum()),changed_types=int((shared.old_type!=shared.new_type).sum()),genes_with_score_changes=int((gd.changed_scores>0).sum()),genes_with_max_changes=int((gd.max_delta.abs()>1e-9).sum()),max_abs_row_delta=float(shared.delta.abs().max()),row_abs_delta_quantiles={str(q):float(shared.delta.abs().quantile(q)) for q in [.5,.9,.95,.99]},positive_to_zero=int(((shared.old_score>0)&(shared.new_score==0)).sum()),zero_to_positive=int(((shared.old_score==0)&(shared.new_score>0)).sum()),spearman_gene_max=float(gd.old_rank.corr(gd.new_rank)),largest_changes=gd.reindex(gd.max_abs_delta.sort_values(ascending=False).index).head(8).to_dict('records'))
for n in [10,20,50]:
 a=set(gd.sort_values(['old_max','gene'],ascending=[False,True]).head(n).gene);b=set(gd.sort_values(['new_max','gene'],ascending=[False,True]).head(n).gene);summary[f'top{n}']=dict(shared=len(a&b),old_only=sorted(a-b),new_only=sorted(b-a))
(out/'master_summary_121.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
