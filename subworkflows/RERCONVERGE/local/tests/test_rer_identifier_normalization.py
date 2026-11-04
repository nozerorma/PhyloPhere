"""RER renaming preserves numerical values, null identity and source artifacts."""
import csv
import os
import subprocess
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[4]
LOCAL = ROOT/'subworkflows/RERCONVERGE/local'
RESOLVER = ROOT/'subworkflows/SCORING/local/src/scoring_io.R'


def run(wd, script, *args, env=None):
    return subprocess.run(['Rscript',str(script),*map(str,args)],cwd=wd,env=env,capture_output=True,text=True)


def make_fixture(path, checksum):
    env=os.environ.copy();env['PHYLOPHERE_VERIFY_CHECKSUMS']=str(checksum)
    code=f'''source("{LOCAL}/rer_batches.R")
r<-data.frame(Rho=c(.2,-.3,.1),N=c(5,5,5),P=c(.2,.3,.4),p.adj=c(.4,.4,.4),row.names=c("KMT2D","AC005551.1","RERONLY"))
attr(r,"rer_settings")<-list(matrix_contract=list(genes=rownames(r)),method="test")
r<-rer_write_batches(r,2,2,0,function(n) list(corRho=matrix(c(.1,.2,.3,.4,.5,.6),3,dimnames=list(rownames(r),NULL)),corStat=matrix(1:6,3,dimnames=list(rownames(r),NULL))),"null")
saveRDS(r,"result.rds")'''
    r=run(path,'-e',code,env=env);assert r.returncode==0,r.stdout+r.stderr
    (path/'reference.txt').write_text('KMT2D.Homo_sapiens\nAC005551.1.Homo_sapiens\n')
    return env


@pytest.mark.parametrize('checksum',[0,1])
def test_renames_result_and_all_null_blocks_without_mutating_source(tmp_path,checksum):
    env=make_fixture(tmp_path,checksum)
    r=run(tmp_path,LOCAL/'rer_tables.R','result.rds','null','trait','.05','.05','reference.txt',RESOLVER,env=env)
    assert r.returncode==0,r.stdout+r.stderr
    with (tmp_path/'rer_gene_identifiers.tsv').open() as f:
        mapping=list(csv.DictReader(f,delimiter='\t'))
    assert [(r['original_gene'],r['gene']) for r in mapping]==[
        ('KMT2D','KMT2D.Homo_sapiens'),('AC005551.1','AC005551.1.Homo_sapiens'),('RERONLY','RERONLY')]
    code='''orig<-readRDS("result.rds"); mapped<-readRDS("validated_rer.output");
stopifnot(identical(rownames(orig),c("KMT2D","AC005551.1","RERONLY")),
          identical(rownames(mapped),c("KMT2D.Homo_sapiens","AC005551.1.Homo_sapiens","RERONLY")),
          identical(orig$Rho,mapped$Rho),identical(orig$P,mapped$P),identical(orig$p.perm,mapped$p.perm),
          identical(attr(mapped,"rer_settings")$matrix_contract$genes,rownames(mapped)))
for (name in c("blk_000001.rds","blk_000002.rds")) {
 x<-readRDS(file.path("null",name))$corStat;y<-readRDS(file.path("normalized_null",name))$corStat
 stopifnot(identical(unname(x),unname(y)),identical(colnames(x),colnames(y)),identical(rownames(y),rownames(mapped)))
}'''
    check=run(tmp_path,'-e',code,env=env);assert check.returncode==0,check.stdout+check.stderr
    regen=tmp_path/'regenerate';regen.mkdir()
    r=run(regen,LOCAL/'rer_tables.R',tmp_path/'validated_rer.output',tmp_path/'normalized_null','trait','.05','.05',env=env)
    assert r.returncode==0,r.stdout+r.stderr


def test_normalized_result_without_null_does_not_claim_old_manifest(tmp_path):
    env=make_fixture(tmp_path,0)
    r=run(tmp_path,LOCAL/'rer_tables.R','result.rds','NO_FILE','trait','.05','.05','reference.txt',RESOLVER,env=env)
    assert r.returncode==0,r.stdout+r.stderr
    check=run(tmp_path,'-e','x<-readRDS("validated_rer.output");stopifnot(is.null(attr(x,"rer_null_manifest")))',env=env)
    assert check.returncode==0,check.stdout+check.stderr
    assert not (tmp_path/'normalized_null').exists()


def test_ambiguous_reference_is_rejected(tmp_path):
    env=make_fixture(tmp_path,0)
    (tmp_path/'reference.txt').write_text('KMT2D.Homo_sapiens\nKMT2D.Papio_anubis\n')
    r=run(tmp_path,LOCAL/'rer_tables.R','result.rds','null','trait','.05','.05','reference.txt',RESOLVER,env=env)
    assert r.returncode!=0 and 'Ambiguous producer gene alias' in r.stderr


def test_nextflow_publishes_matching_normalized_null_and_summary(tmp_path):
    import shutil
    if not shutil.which('nextflow'):
        pytest.skip('Nextflow unavailable')
    env=make_fixture(tmp_path,1)
    (tmp_path/'subworkflows').symlink_to(ROOT/'subworkflows',target_is_directory=True)
    (tmp_path/'main.nf').write_text('''params.outdir = 'out'
params.traitname = 'test'
params.evidence_p_threshold = 0.05
params.evidence_fdr_threshold = 0.05
params.use_singularity = false
params.use_apptainer = false
include { RER_TABLES } from './subworkflows/RERCONVERGE/rer_tables.nf'
workflow {
  RER_TABLES(Channel.value(file('result.rds')),Channel.value(file('null')),
             Channel.value(file('reference.txt')),Channel.value('test-sources'))
}
''')
    r=subprocess.run(['nextflow','run','main.nf','-ansi-log','false'],cwd=tmp_path,env=env,capture_output=True,text=True,timeout=180)
    assert r.returncode==0,r.stdout+r.stderr
    out=tmp_path/'out/rerconverge/rer_results'
    assert 'KMT2D.Homo_sapiens' in (out/'rerconverge_summary_test.tsv').read_text()
    assert (out/'normalized_null/rer_null_manifest.json').is_file()
    check=run(tmp_path,'-e','''x<-readRDS("out/rerconverge/rer_results/normalized_null/blk_000001.rds");stopifnot(rownames(x$corStat)[1]=="KMT2D.Homo_sapiens")''',env=env)
    assert check.returncode==0,check.stdout+check.stderr
