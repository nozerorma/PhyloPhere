"""WP8 migration, saved templates, isolated cleanup and empty DOMINO contracts."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
ENV = dict(os.environ, NXF_OFFLINE='true')

@pytest.mark.parametrize('policy,expected', [('none', False), ('prune_intermediates', False), ('remove_work', True)])
def test_cleanup_migration_preserves_work_choice(policy, expected):
    from gui.models.serialization import from_dict, to_dict
    p = from_dict({'runtime': {'cleanup_policy': policy, 'remove_intermediates': True}})
    assert p.runtime.remove_work is expected
    assert 'remove_intermediates' not in to_dict(p)['runtime']
    p = from_dict({'runtime': {'cleanup_policy': policy, 'remove_work': not expected}})
    assert p.runtime.remove_work is not expected


def test_all_templates_load_save_and_render(tmp_path):
    from gui.project_io import load_project, save_project
    from gui.generation.render import render_single, render_batch
    paths = sorted((ROOT/'gui/templates').glob('*.json')) + sorted((ROOT/'validation').glob('tier*/templates/*.json'))
    for i, path in enumerate(paths):
        p = load_project(path)
        target = tmp_path/f'{i}.json'
        save_project(target, p)
        assert load_project(target) == p
        for render in (render_single, render_batch):
            text = render(p)
            assert 'REMOVE_INTERMEDIATES' not in text and '"remove_intermediates"' not in text
            result = subprocess.run(['bash', '-n'], input=text, text=True, capture_output=True)
            assert result.returncode == 0, (path, result.stderr)


def test_bool_exports_are_json_safe_lowercase():
    """Regression: a raw Python bool in a template export (RUN_POSENRICH in the
    sbatch template) reached params.json as True/False and broke Nextflow's JSON
    parser. Every boolean export must render as the lowercase string, and no
    capitalized True/False may reach any generated script."""
    from gui.project_io import load_project
    from gui.generation.render import render_single, render_batch
    project = load_project(ROOT / 'gui/templates/neoplasia_toy.json')
    assert project.modules.posenrich.enabled is True
    batch = render_batch(project)
    assert 'export RUN_POSENRICH=true' in batch
    single = render_single(project)
    assert '"posenrich": $_RUN_POSENRICH,' in single
    for text, name in ((batch, 'sbatch_array'), (single, 'run_single')):
        assert '="True"' not in text, f'{name}: capitalized boolean leaked'
        assert '="False"' not in text, f'{name}: capitalized boolean leaked'


def runner_tail():
    from gui.project_io import load_project
    from gui.generation.render import render_single
    text = render_single(load_project(ROOT/'gui/templates/tier1_pepc_c4.json'))
    return text[text.index('nextflow run "${REPO_DIR}/main.nf"'):]

@pytest.mark.parametrize('success,remove', [(True, True), (False, True), (True, False)])
def test_real_nextflow_cleanup_isolated(tmp_path, success, remove):
    """Execute the rendered launcher's completion path with a real one-process Nextflow run."""
    if not shutil.which('nextflow'): pytest.skip('Nextflow required')
    code = tmp_path/'code';code.mkdir()
    base = tmp_path/'scratch';base.mkdir()
    work = base/'work';run = base/'run';home = base/'home'
    for d in (work, run, home): d.mkdir()
    out = tmp_path/'published';out.mkdir()
    (out/'canonical.rds').write_bytes(b'canonical null fixture\n')
    before = hashlib.sha256((out/'canonical.rds').read_bytes()).hexdigest()
    (work/'payload').write_bytes(b'x'*1048576)
    # Use the installed runtime/plugins cache through a symlink; cleanup must retain the cache.
    default_home = Path(os.environ.get('NXF_HOME', str(Path.home()/'.nextflow')))
    if (default_home/'framework').exists(): (home/'framework').symlink_to(default_home/'framework', target_is_directory=True)
    if (default_home/'plugins').exists(): (home/'plugins').symlink_to(default_home/'plugins', target_is_directory=True)
    command = 'echo published > result.txt' if success else 'exit 17'
    (code/'main.nf').write_text('''nextflow.enable.dsl=2
process CHECK {
 publishDir params.outdir, mode: 'copy'
 output:
 path 'result.txt'
 script:
 """COMMAND"""
}
workflow { CHECK() }
'''.replace('COMMAND', command))
    env = dict(ENV, REPO_DIR=str(code), WORK_BASE=str(base), WORK_DIR=str(work), RUN_DIR=str(run), NXF_HOME=str(home),
               MEASURE_DIR=str(tmp_path), RESULTS_BASE=str(out), TRAIT='fixture', RUN_TAG='', _REMOVE_WORK=str(remove).lower(), TMPDIR=str(tmp_path))
    tail = runner_tail().replace('    # Validate every target',
        '    du -s -B1 "$WORK_DIR" "$RUN_DIR" "$NXF_HOME" > "$MEASURE_DIR/du_before.tsv"\n'
        '    find "$RESULTS_BASE" -type f -exec sha256sum {} + > "$MEASURE_DIR/published_before.sha256"\n'
        '    # Validate every target')
    script='set -Eeuo pipefail\nNF_ENGINE=(-work-dir "$WORK_DIR" --outdir "$RESULTS_BASE" -ansi-log false)\n'+tail
    result = subprocess.run(['bash', '-c', script], cwd=run, env=env, text=True, capture_output=True, timeout=100)
    (tmp_path/'launcher.log').write_text(result.stdout+result.stderr)
    assert (result.returncode == 0) is success, result.stdout+result.stderr
    assert hashlib.sha256((out/'canonical.rds').read_bytes()).hexdigest() == before
    assert work.exists() is (not (success and remove))
    assert run.exists() is work.exists() and home.exists() is work.exists()
    if success: assert (out/'result.txt').read_text().strip() == 'published'
    allocated = None
    if success and remove:
        allocated = sum(int(line.split()[0]) for line in (tmp_path/'du_before.tsv').read_text().splitlines())
        assert allocated >= 1048576
        verified = subprocess.run(['sha256sum','-c',str(tmp_path/'published_before.sha256')],capture_output=True,text=True)
        assert verified.returncode == 0, verified.stdout+verified.stderr
    (tmp_path/'measurement.json').write_text(json.dumps({'success':success,'remove_work':remove,'fixture_payload_bytes':1048576,
        'allocated_bytes_recovered':allocated,'payload_bytes_recovered':1048576 if not work.exists() else 0,'published_sha256_before':before,
        'published_sha256_after':hashlib.sha256((out/'canonical.rds').read_bytes()).hexdigest()}))


def test_cleanup_rejects_results_inside_work(tmp_path):
    base=tmp_path/'base';work=base/'work';out=work/'results';out.mkdir(parents=True)
    run=base/'run';run.mkdir();home=base/'home';home.mkdir()
    env=dict(os.environ,WORK_BASE=str(base),WORK_DIR=str(work),RUN_DIR=str(run),NXF_HOME=str(home),RESULTS_BASE=str(out),
             REPO_DIR=str(ROOT),TRAIT='guard',RUN_TAG='',_REMOVE_WORK='true')
    tail=runner_tail();tail=tail[tail.index('cleanup_work()'):]
    result=subprocess.run(['bash','-c','set -euo pipefail\n'+tail],env=env,capture_output=True,text=True)
    assert result.returncode != 0 and 'Refusing' in result.stderr
    assert out.exists() and run.exists() and home.exists()


def test_domino_zero_edges_and_empty_report_graph(tmp_path):
    db=tmp_path/'db';db.mkdir()
    for kind, content in [('links','protein1 protein2 combined_score\np1 p2 600\n'),
                          ('info','#string_protein_id\tpreferred_name\tprotein_size\tannotation\np1\tG1\t1\tx\np2\tG2\t1\tx\n')]:
        with gzip.open(db/f'9606.protein.{kind}.v12.0.txt.gz','wt') as f:f.write(content)
    bg=tmp_path/'background';bg.write_text('G1\nG2\n')
    out=tmp_path/'out'
    result=subprocess.run(['python',str(ROOT/'subworkflows/ENRICHMENT/local/src/build_domino_network.py'),
        '--cache-dir',str(tmp_path/'cache'),'--string-db-dir',str(db),'--cleaned-background',str(bg),'--output-dir',str(out)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    for name in ('network.sif','network_edge_scores.tsv','slices.txt'):assert (out/name).stat().st_size==0
    direct=subprocess.run(['python',str(ROOT/'subworkflows/ENRICHMENT/local/src/run_domino_modules.py'),
        '--network',str(out/'network.sif'),'--slices',str(out/'slices.txt'),
        '--gene-lists-dir',str(tmp_path/'absent_lists'),'--output-dir',str(tmp_path/'modules')],capture_output=True,text=True)
    assert direct.returncode==0,direct.stderr
    assert list((tmp_path/'modules').iterdir())==[]
    if not shutil.which('Rscript'):pytest.skip('Rscript required')
    report=(ROOT/'subworkflows/ENRICHMENT/local/13.AMI_analysis.Rmd').read_text()
    function=report[report.index('build_tool_network <-'):report.index('# Per-gene flag table')]
    result=subprocess.run(['Rscript','-e',function+f'\nstopifnot(igraph::vcount(build_tool_network("{out}/network.sif","{out}/network_edge_scores.tsv"))==0)'],capture_output=True,text=True)
    assert result.returncode==0,result.stderr


def test_retired_selectors_migrate_without_changing_live_resources():
    from gui.models.serialization import from_dict
    rows=[{'selector':'CAAS_PERMS_REBUILD','cpus':'99'},{'selector':'CAAS_CORE_MERGE','cpus':'7'}]
    p=from_dict({'resources':{'process_overrides':rows}})
    assert len(p.resources.process_overrides)==1
    assert p.resources.process_overrides[0].selector=='CAAS_CORE_MERGE'
    assert p.resources.process_overrides[0].cpus=='7'


def test_acceptance_renderer_preserves_seed_and_collects_resources(tmp_path):
    from gui.project_io import load_project
    template=ROOT/'gui/templates/neoplasia_def.json'
    saved=load_project(template)
    result=subprocess.run(['python',str(ROOT/'validation/tier1/scripts/render_tier1_scripts.py'),
        '--template',str(template),'--outdir',str(tmp_path),'--caas-full-perms','1000',
        '--rer-perm-batches','10','--asr-cache',str(tmp_path/'frozen_asr'),'--keep-work','--collect-resources'],
        capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    resolved=load_project(tmp_path/'resolved_project.json')
    assert resolved.modules.rer.rer_batch_seed == saved.modules.rer.rer_batch_seed
    assert resolved.modules.caas.caas_full_perms=='1000' and not resolved.runtime.remove_work
    assert resolved.modules.disambiguation.ct_disambig_asr_cache_dir==str(tmp_path/'frozen_asr')
    single=next(tmp_path.glob('*single*.sh')).read_text()
    assert '-with-trace ' in single and '-with-report ' in single and '-with-timeline ' in single
    for p in tmp_path.glob('*.sh'):
        assert subprocess.run(['bash','-n',str(p)]).returncode==0
