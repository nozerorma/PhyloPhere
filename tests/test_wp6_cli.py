"""Input resolution and failure contracts for the productive block path."""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
R = os.environ.get('RSCRIPT', shutil.which('Rscript') or '')
PREPARE = ROOT / 'subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R'


def run_r(tmp_path, *args):
    if not R:
        pytest.skip('Rscript required')
    return subprocess.run([R, *map(str, args)], cwd=tmp_path,
                          env=dict(os.environ, TMPDIR=str(tmp_path)), capture_output=True, text=True)


@pytest.mark.parametrize('source', ['caas', 'rer'])
def test_no_null_has_explicit_output(tmp_path, source):
    got = run_r(tmp_path, PREPARE, source, 'size_adj_max', 7, 'enabled')
    assert got.returncode == 0, got.stderr
    if source == 'caas':
        assert (tmp_path / 'empty.tsv').read_text().startswith('ranking\tdatabase\tpathway')
    else:
        assert (tmp_path / 'blocks/no_null.rds').exists()


def test_legacy_import_splits_columns(tmp_path):
    generated = run_r(tmp_path, '-e', '''
        m<-matrix(1,4,15,dimnames=list(paste0('g',1:4),paste0('b_',1:15)))
        saveRDS(list(caas_corStat_byrank=list(global=m,top=m,bottom=m),gene_stat='max'),'null_input')
    ''')
    assert generated.returncode == 0, generated.stderr
    got = run_r(tmp_path, PREPARE, 'caas', 'max', 7, 'enabled')
    assert got.returncode == 0, got.stderr
    assert len(list((tmp_path / 'blocks').glob('*.rds'))) == 3


def test_mismatched_gene_stat_is_rejected(tmp_path):
    run_r(tmp_path, '-e', "saveRDS(list(gene_stat='max'),'null_input')")
    got = run_r(tmp_path, PREPARE, 'caas', 'fmax', 7, 'enabled')
    assert got.returncode != 0
    assert 'gene_stat mismatch' in got.stderr


def test_rer_overlapping_directory_is_rejected(tmp_path):
    run_r(tmp_path, '-e', '''
      dir.create('null_input')
      m<-matrix(1,4,2,dimnames=list(paste0('g',1:4),c('p1','p2')))
      saveRDS(list(corStat=m),'null_input/a.rds');saveRDS(list(corStat=m),'null_input/b.rds')
    ''')
    got = run_r(tmp_path, PREPARE, 'rer', 'max', 7, 'enabled')
    assert got.returncode != 0
    assert 'overlapping' in got.stderr
