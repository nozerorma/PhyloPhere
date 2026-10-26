"""The sole CAAS maximum: direction semantics, CLI, project migration and stale-null rejection."""
import csv
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import pytest

ROOT=Path(__file__).parents[1]
LOCAL=ROOT/'subworkflows/CT_DISAMBIGUATION/local'
spec=importlib.util.spec_from_file_location('max_only_scores',LOCAL/'src/core/scores.py')
scores=importlib.util.module_from_spec(spec);spec.loader.exec_module(scores)
OBS=ROOT/'subworkflows/SCORING/local/src/observed_core_scores.py'


def test_maximum_is_independent_of_pool_and_position_multiplicity():
    assert scores.GENE_STATS==('max',)
    assert scores.gene_stat_value([.2,.7,.7],[])==.7
    assert scores.gene_stat_value([.2,.7],[.99]*200)==.7
    assert scores.gene_stat_value([None,0,float('nan')],[.5]) is None
    got=scores.gene_scores([{'top':.7,'bottom':.2},{'bottom':.9}],{'all':[],'top':[],'bottom':[]})
    assert got=={'all':.9,'top':.7,'bottom':.9}
    for retired in ('fmax','size_adj_max'):
        with pytest.raises(ValueError,match='Only'):scores.gene_stat_value([.7],[.7],retired)
    assert not hasattr(scores,'size_adj_max')


def test_observed_cli_uses_raw_maximum_and_rejects_retired_options(tmp_path):
    inp=tmp_path/'positions.tsv'
    inp.write_text('Gene\tPosition\tside\tasr_path_score\ng\t0\ttop\t.7\ng\t0\tbottom\t.2\ng\t1\tbottom\t.9\nh\t0\ttop\t0\n')
    args=[sys.executable,str(OBS),'--input',str(inp),'--positions-out',str(tmp_path/'pos.tsv'),'--genes-out',str(tmp_path/'genes.tsv')]
    got=subprocess.run(args,text=True,capture_output=True)
    assert got.returncode==0,got.stderr
    with (tmp_path/'genes.tsv').open() as f:rows={x['Gene']:x for x in csv.DictReader(f,delimiter='\t')}
    assert [float(rows['g'][k]) for k in ('gene_caas_score','gene_caas_score_top','gene_caas_score_bottom')]==[.9,.7,.9]
    assert rows['h']['gene_caas_score']=='NA'
    for retired in ('fmax','size_adj_max'):
        bad=subprocess.run(args+['--gene-stat',retired],text=True,capture_output=True)
        assert bad.returncode!=0 and 'invalid choice' in bad.stderr


def test_saved_projects_migrate_but_in_memory_overrides_are_rejected(caplog):
    from gui.models.serialization import from_dict,to_dict
    from gui.generation.context import build_context
    from gui.generation.validate import validate
    from gui.generation.render import render_single,render_batch
    for retired in ('fmax','size_adj_max'):
        p=from_dict({'modules':{'caas':{'gene_stat':retired}}})
        assert p.modules.caas.gene_stat=='max'
        assert to_dict(p)['modules']['caas']['gene_stat']=='max'
        assert 'legacy gene scores/null matrices must be rebuilt' in caplog.text
        assert 'export GENE_STAT="max"' in render_batch(p)
        assert '${GENE_STAT:-max}' in render_single(p)
        p.modules.caas.gene_stat=retired
        assert any('only gene_stat=max' in x for x in validate(p))
        with pytest.raises(ValueError,match='Only'):build_context(p)


@pytest.mark.parametrize('stamp', ['size_adj_max','fmax',None])
def test_legacy_null_is_rejected_even_when_runtime_requests_max(tmp_path,stamp):
    r=shutil.which('Rscript')
    if not r:pytest.skip('Rscript required')
    label='NULL' if stamp is None else repr(stamp)
    code=f"saveRDS(list(gene_stat={label},caas_corStat_byrank=list(global=matrix(.7,1,3))), 'null_input')"
    gen=subprocess.run([r,'-e',code],cwd=tmp_path,capture_output=True,text=True)
    assert gen.returncode==0,gen.stderr
    run=subprocess.run([r,str(ROOT/'subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R'),'caas','max','7','enabled'],cwd=tmp_path,capture_output=True,text=True)
    assert run.returncode!=0 and ('gene_stat mismatch' in run.stderr or 'lacks gene_stat stamp' in run.stderr)
