"""Compare full observed scoring and axes replay with fixed local inputs."""
import csv,json,logging,pickle,re,sys
from pathlib import Path
from collections import defaultdict
sys.path.insert(0,str(Path('subworkflows/CT_DISAMBIGUATION/local').resolve()))
logging.basicConfig(level=logging.ERROR)
from src.core.driver import score_labelings,pool_labelings
from src.core.observed import observed_entries,analyze_observed
from src.core.labelings import read_trait_pairs,observed_pss
from src.convergence.disambiguate_single import analyze_gene_disambiguation
p=Path('.tmp/asr_deep_audit/inputs').resolve();g=sys.argv[1] if len(sys.argv)>1 else 'A3GALT2'
ctxfile=p/('new_context.pkl' if g=='A1BG' else f'{g}_new_context.pkl')
ctx=pickle.loads(ctxfile.read_bytes());gene=g+'.Homo_sapiens'
raw=list(csv.DictReader((p/f'{g}.discovery.tsv').open(),delimiter='\t'))
pairs=read_trait_pairs(p/'traitfiles_ok_dir');weights=observed_pss(p/'traitfiles_ok_dir/contrast_hypotheses_pairs.tsv')
bytag=defaultdict(list)
for r in raw:
 h=int(re.search(r'H(\d+)',r['trait'])[1]);tag=f'b_0~H{h}';e=observed_entries(gene,[r])[0];e.trait='';bytag[tag].append(e)
labels={f'b_0~H{h}':(list(zip(*ps))[0],list(zip(*ps))[1]) for h,ps in pairs.items()}
scored=score_labelings(ctx,gene,sorted(bytag),labels,bytag,0.1)
pooled=pool_labelings(scored,{'b_0':weights},labels)
full,diag=analyze_observed(ctx,gene,observed_entries(gene,raw),pairs,weights,0.1,keep_unpooled=True)
def tab(rows):return {(r.position,r.side):r.asr_path_score for r in rows}
a=tab(full);b=tab([r for c,rs in pooled for r in rs]);d=[(k,a.get(k),b.get(k)) for k in a.keys()|b.keys() if k not in a or k not in b or abs(a[k]-b[k])>1e-10]
print(g,'full/replay fixed discovery',len(a),len(b),d)
Path('.tmp/asr_deep_audit/'+g+'_replay_debug.json').write_text(json.dumps({'full':[{ 'position':k[0],'side':k[1],'score':v} for k,v in a.items()],'axes':[{ 'position':k[0],'side':k[1],'score':v} for k,v in b.items()]},indent=2))
from src.utils.perm_stream import ReplaySettings,ReplayContext,stream_gene
block=p/'audit_b0';block.mkdir(exist_ok=True)
with (block/'resample_perms.tab').open('w') as f:
 for tag,(fg,bg) in labels.items():f.write(tag+'\t'+','.join(fg)+'\t'+','.join(bg)+'\n')
with (block/'fop_pairs.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['cycle','hypothesis_id','pair','species1','species2','pss_score'])
 for (h,d),(fg,bg,pss) in weights.items():w.writerow(['b_0',h,d,fg,bg,pss])
settings=ReplaySettings.from_args(['--patterns','1,2,3','--miss_pair','--max_conserved','3','--min_non_overlapping','2','--max_fg_gaps','0','--max_bg_gaps','0','--max_gaps','0','--max_fg_miss','0','--max_bg_miss','0','--max_miss','0'],str(Path('subworkflows/CT/local').resolve()),str(p/'traitfiles_ok_dir'),'fasta')
replay=ReplayContext(settings,str(block),fop_pairs_file=str(block/'fop_pairs.tsv'))
pooled,_=stream_gene(gene,str(p/(gene+'.fa')),replay,lambda:ctx,0.1)
c=tab([r for cy,rs in pooled for r in rs]);print(g,'actual kernel replay vs full',[(k,a.get(k),c.get(k)) for k in a.keys()|c.keys() if k not in a or k not in c or abs(a[k]-c[k])>1e-10])
