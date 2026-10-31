"""Canonical positional fixtures with literal IDs, explicit references and null provenance."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def seal_block(block):
    (block/'scoring_artifacts.json').write_text(json.dumps(dict(version=1, files={
        n: hashlib.md5((block/n).read_bytes()).hexdigest()
        for n in ['matrices.rds','positions.parquet','cycles.txt','counts.json']})))


def universe_provenance(tested):
    """The `universe` entry of a null's provenance: the digest of every line of a background.output (SCORING checks it)."""
    out = {}
    for line in Path(tested).read_text().splitlines():
        gene, _, text = line.partition('\t')
        positions = [] if text in ('', 'NULL') else [int(p) for p in text.split(',')]
        out[gene] = hashlib.sha256(f"{gene}\t{','.join(map(str, sorted(positions))) if positions else 'NULL'}".encode()).hexdigest()
    return out


def canonical(root, observed, genes, rows=(), cycles=None, real_matrix=False, tested=None):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    cycles=cycles or [f'b_{i}' for i in range(1,6)]
    prov=dict(version=1, scientific=dict(model='LG',tree_sha256='f'*64,posterior_threshold=.8,taxid_sha256=None,weights_sha256=None),
              cluster=dict(enabled=True,remove=True,minlen=3,maxcaas=.7,min_positions=None),
              gene_filter=dict(mode='none'),strict=False,maps={},lengths={})
    Path(str(observed)+'.scoring.json').write_text(json.dumps(dict(prov, observed_sha256=digest(observed))))
    if tested is not None:prov['universe']=universe_provenance(tested)
    block=root/'blocks/blk_000001';block.mkdir(parents=True,exist_ok=True)
    (root/'cycle_roster.txt').write_text('\n'.join(cycles)+'\n')
    (block/'cycles.txt').write_text('\n'.join(cycles)+'\n')
    df=pd.DataFrame(rows,columns=['gene','cycle','position','top','bottom'])
    df=df.astype({'gene':'str','cycle':'str','position':'int32','top':'float64','bottom':'float64'})
    df['all']=np.fmax(df.top,df.bottom)
    pq.write_table(pa.Table.from_pandas(df,preserve_index=False).cast(pa.schema([('gene',pa.string()),('cycle',pa.string()),('position',pa.int32()),('top',pa.float64()),('bottom',pa.float64()),('all',pa.float64())])),block/'positions.parquet')
    aggregate=df.groupby(['gene','position'],sort=False).size().rename('nd').reset_index()
    first={}
    for r in df.itertuples(): first[(r.gene,r.position)]=min(first.get((r.gene,r.position),len(cycles)),cycles.index(r.cycle))
    aggregate['first_rank']=[first[(r.gene,r.position)] for r in aggregate.itertuples()]
    aggregate=aggregate.astype({'gene':'str','position':'int32','nd':'int32','first_rank':'int32'})
    pq.write_table(pa.Table.from_pandas(aggregate,preserve_index=False).cast(pa.schema([('gene',pa.string()),('position',pa.int32()),('nd',pa.int32()),('first_rank',pa.int32())])),root/'null_positions.parquet')
    counts=dict(cycles=len(cycles),cycles_present=sorted(set(df.cycle)),scored_positions=len(df),score_rule='caas_pos_species_pss_v4',scoring_provenance=prov)
    (block/'counts.json').write_text(json.dumps(counts))
    (root/'counts.json').write_text(json.dumps(dict(counts,cycles_present=len(set(df.cycle)))))
    if real_matrix:
        (root/'genes.txt').write_text('\n'.join(genes)+'\n')
        df[['gene','cycle','top','bottom','all']].to_csv(root/'matrix_input.tsv',sep='\t',index=False,na_rep='NA',quoting=csv.QUOTE_NONE)
        code="""a<-commandArgs(TRUE);g<-readLines(file.path(a[1],'genes.txt'));cy<-readLines(file.path(a[1],'cycle_roster.txt'));x<-read.delim(file.path(a[1],'matrix_input.tsv'),colClasses='character',quote='',na.strings=NULL);m<-matrix(0,length(g),length(cy),dimnames=list(g,cy));ss<-list();for(d in c('all','top','bottom')){v<-m;for(i in seq_len(nrow(x))){z<-suppressWarnings(as.numeric(x[[d]][i]));if(!is.na(z))v[x$gene[i],x$cycle[i]]<-max(v[x$gene[i],x$cycle[i]],z)};ss[[if(d=='all')'global' else d]]<-v};saveRDS(list(genes=g,cycles=cy,gene_stat='max',stat=ss),file.path(a[1],'blocks/blk_000001/matrices.rds'))"""
        result=subprocess.run(['Rscript','-e',code,str(root)],capture_output=True,text=True)
        assert result.returncode==0,result.stderr
    else:
        (block/'matrices.rds').write_bytes(b'checksum-only matrix fixture; real workflow fixtures use RDS')
    seal_block(block)
    return block


def reference(path,gene,sequence='A'*12,strand='+'):
    rows=list(csv.DictReader(Path(path).open(),delimiter='\t'))
    codons={}
    for row in rows:
        if row['status']=='selected' and row['hg38_aa_pos'] not in ('NA','-','') and row['hg38_nt_coord'] not in ('NA','-',''):
            chrom,nt=row['hg38_nt_coord'].split(':');nt=int(nt);step=1 if strand=='+' else -1
            codons[row['hg38_aa_pos']]=dict(positions=[f'{chrom}:{nt+step*i}' for i in range(3)],coding_bases='GCT')
    meta=dict(schema='phylophere.vep_reference.v1',gene=gene,species='homo_sapiens',assembly='GRCh38',annotation_release=115,
              source_uri='fixture:pinned-reference',map_sha256=digest(path),gene_id='ENSG000001.1',transcript_id='ENST000001.1',protein_id='ENSP000001.1',
              protein_sequence=sequence,protein_sha256=hashlib.sha256(sequence.encode()).hexdigest(),strand=strand,codons=codons)
    Path(str(path)+'.vep.json').write_text(json.dumps(meta))
    return {k:meta[k] for k in ('protein_id','transcript_id','gene_id','annotation_release','assembly','species')}


def annotation(path, genes, identities=None, tested=None, background=None, coverage='tested_positions_in_declared_genes'):
    meta=dict(table_sha256=digest(path),coordinate_space='CAAS_0based',source_uri='fixture:prepared',version='1',coverage=coverage,coverage_genes=list(genes),reference_identities=identities or {})
    if tested and background:
        meta.update(tested_sha256=digest(tested),background_sha256=digest(background))
    else:
        # Isolated descriptive-reader fixtures do not claim validated inference.
        meta['coverage']=''
    Path(str(path)+'.manifest.json').write_text(json.dumps(meta))
    return meta


def calibrate_scores(scores, null_root, tested, genes):
    """Small independent fixture oracle: max-side score, inclusive rank, fixed BH family."""
    null_root=Path(null_root)
    blocks=[pd.read_parquet(b/'positions.parquet') for b in (null_root/'blocks').iterdir()]
    df=pd.concat(blocks,ignore_index=True)
    table=pd.read_csv(scores,sep='\t',dtype=str,keep_default_na=False,quoting=csv.QUOTE_NONE)
    best=table.assign(CAAS_score=table.CAAS_score.astype(float)).groupby(['Gene','Position']).CAAS_score.max()
    N=len((null_root/'cycle_roster.txt').read_text().splitlines())
    M=0
    for line in Path(tested).read_text().splitlines():
        g,ps=line.split('\t')
        if g in genes:M+=len(ps.split(',')) if ps else 0
    ps,nds={},{}
    for (g,p),score in best.items():
        vals=df[(df.gene==g)&(df.position==int(p))]['all']
        ps[g,p]=(1+int((vals>=score-1e-12).sum()))/(N+1)
        nds[g,p]=len(vals)
    keys=sorted(ps,key=lambda k:ps[k]);q={};running=1.
    for i in range(len(keys)-1,-1,-1):
        k=keys[i];running=min(running,ps[k]*M/(i+1));q[k]=running
    table['p.emp']=[ps[g,p] for g,p in zip(table.Gene,table.Position)]
    table['p.adj_bh']=[q[g,p] for g,p in zip(table.Gene,table.Position)]
    table['nd']=[nds[g,p] for g,p in zip(table.Gene,table.Position)]
    table['null_unseen']=['TRUE' if nd==0 else 'FALSE' for nd in table.nd]
    table['caas_nominal']=['TRUE' if v<=.05 else 'FALSE' for v in table['p.emp']]
    table['caas_adjusted']=['TRUE' if v<=.05 else 'FALSE' for v in table['p.adj_bh']]
    table.to_csv(scores,sep='\t',index=False,quoting=csv.QUOTE_NONE)
