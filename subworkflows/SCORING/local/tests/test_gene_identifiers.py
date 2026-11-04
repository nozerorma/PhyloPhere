"""Exercise integration of precomputed producer IDs without rewriting CAAS nulls."""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / 'subworkflows/SCORING/local/src'


def run_scoring(tmp_path, fade_extra=''):
    for name in ('scoring_compute.R','scoring_io.R','scoring_input.py','aa_grouping.R','null_position_p.R','observed_core_scores.py'):
        shutil.copy(SRC / name, tmp_path / name)
    core = tmp_path/'src/core'; core.mkdir(parents=True)
    shutil.copy(ROOT/'subworkflows/CT_DISAMBIGUATION/local/src/core/scores.py', core/'scores.py')
    (core/'__init__.py').touch();(core.parent/'__init__.py').touch()
    (tmp_path/'observed.tsv').write_text('Gene\tPosition\tside\tasr_path_score\nKMT2D.Homo_sapiens\t0\ttop\t0.8\nAC005551.1.Homo_sapiens\t0\ttop\t0.7\n')
    (tmp_path/'background.txt').write_text('KMT2D.Homo_sapiens\nAC005551.1.Homo_sapiens\nOTHER.Homo_sapiens\n')
    (tmp_path/'fade.tsv').write_text('gene\tmax_bf\nKMT2D\t120\nOTHER\t110\n'+fade_extra)
    (tmp_path/'rer.tsv').write_text('gene\tRho\tP\tp.perm\tp.perm.adj\nKMT2D\t0.5\t0.001\t0.03\t0.09\nOTHER\t0.2\t0.01\t0.1\t0.15\nAC005551.1\t0.1\t0.1\t0.2\t0.2\n')
    env = os.environ.copy();env['OPENBLAS_NUM_THREADS']='1'
    r=subprocess.run([sys.executable,'observed_core_scores.py','--input','observed.tsv','--positions-out','pos.tsv','--genes-out','genes.tsv'],cwd=tmp_path,env=env,capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    r=subprocess.run(['Rscript','scoring_compute.R','--postproc','observed.tsv','--core_positions','pos.tsv','--core_genes','genes.tsv','--background','background.txt','--fade_top','fade.tsv','--rer','rer.tsv'],cwd=tmp_path,env=env,capture_output=True,text=True)
    return r


def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def test_unifies_producer_ids_and_keeps_rer_null_namespace(tmp_path):
    r=run_scoring(tmp_path);assert r.returncode==0,r.stdout+r.stderr
    scores={x['Gene']:x for x in rows(tmp_path/'gene_scores.tsv')}
    assert set(scores)=={'KMT2D.Homo_sapiens','AC005551.1.Homo_sapiens','OTHER.Homo_sapiens'}
    assert scores['KMT2D.Homo_sapiens']['fade_significant_top']=='TRUE'
    assert float(scores['KMT2D.Homo_sapiens']['rer_min_pval'])==0.03
    assert scores['OTHER.Homo_sapiens']['gene_caas_score']=='NA'
    assert scores['OTHER.Homo_sapiens']['fade_significant_top']=='TRUE'
    assert {x['gene'] for x in rows(tmp_path/'fcs_stats_rer.tsv')}=={'KMT2D','OTHER','AC005551.1'}
    # The cleaned CAAS universe includes OTHER without an observed hit: calibration
    # must add its CAAS zero to the existing FADE/RER row, never create a duplicate.
    (tmp_path/'blocks').mkdir()
    evidence=ROOT/'subworkflows/SCORING/local/src/scoring_evidence.R'
    r=subprocess.run(['Rscript',str(evidence),'--genes','gene_scores.tsv','--positions','position_scores.tsv',
                      '--stats','fcs_stats.tsv','--universe','background.txt','--blocks','blocks'],
                     cwd=tmp_path,capture_output=True,text=True)
    assert r.returncode==0,r.stdout+r.stderr
    final={x['Gene']:x for x in rows(tmp_path/'gene_scores.tsv')}
    assert len(final)==3 and float(final['OTHER.Homo_sapiens']['gene_caas_score'])==0
    assert final['OTHER.Homo_sapiens']['fade_significant_top']=='TRUE'


def test_colliding_producer_aliases_rejected(tmp_path):
    r=run_scoring(tmp_path,'KMT2D.Homo_sapiens\t140\n')
    assert r.returncode!=0 and 'Ambiguous producer gene alias' in r.stderr


def test_numeric_part_is_not_an_alias(tmp_path):
    expression = (f"source('{SRC / 'scoring_io.R'}'); "
                  "stopifnot(is.na(scoring_match_genes('AC005551','AC005551.1.Homo_sapiens')), "
                  "scoring_match_genes('AC005551.1','AC005551.1.Homo_sapiens') == 1L, is.na(scoring_match_genes('X.Homo_sapiens','X.Papio_anubis')))")
    r=subprocess.run(['Rscript','-e',expression],cwd=tmp_path,capture_output=True,text=True)
    assert r.returncode==0,r.stdout+r.stderr
