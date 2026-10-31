"""Nextflow regression tests for REPORTING staging, ownership and prune-list caching.

Minimal Rmd bodies exercise the actual process scripts without expensive figures.
The real R helpers and report chunks are covered by test_reporting.R.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / 'subworkflows/TRAIT_ANALYSIS/local'


@pytest.fixture
def reporting_project(tmp_path):
    if not shutil.which('nextflow') or not shutil.which('Rscript'):
        pytest.skip('Run with Nextflow and Rscript from the phylophere environment')
    local = tmp_path / 'subworkflows/TRAIT_ANALYSIS/local'
    local.mkdir(parents=True)
    shutil.copytree(LOCAL / 'src', local / 'src')
    bodies = {
        '0.Data_pruning.Rmd': '''
stopifnot(file.exists(params$prune_list))
if(nzchar(params$prune_list_secondary))stopifnot(file.exists(params$prune_list_secondary))
dir.create('data_exploration/0.Data-pruning',recursive=TRUE)
writeLines(readLines(params$prune_list),'data_exploration/prune_marker.txt')
file.copy(params$tree_file,'data_exploration/0.Data-pruning/pruned_tree_file.nwk')
write.table(read.csv(params$trait_file),'data_exploration/0.Data-pruning/pruned_trait_file.tsv',sep='\\t',row.names=FALSE)
write.csv(data.frame(species='S1',score=1),'data_exploration/0.Data-pruning/pruned_trait_stats.csv',row.names=FALSE)
''',
        '1.Dataset_exploration.Rmd': '''
source('src/io_utils.R')
dir.create('data_exploration/1.Data-exploration/1.Species_distribution',recursive=TRUE,showWarnings=FALSE)
write_trait_archive(params$trait_file,'data_exploration/1.Data-exploration/1.Species_distribution/original_trait_file.tsv')
write.csv(data.frame(species='S1',score=1),'data_exploration/1.Data-exploration/1.Species_distribution/trait_stats.csv',row.names=FALSE)
writeLines('dataset','data_exploration/dataset_marker.txt')
''',
        '2.Phenotype_exploration.Rmd': '''
stopifnot(file.exists('data_exploration/dataset_marker.txt'))
writeLines('phenotype','data_exploration/phenotype_marker.txt')
''',
    }
    for name, body in bodies.items():
        source = (LOCAL / name).read_text()
        header = source[:source.index('\n---', 4) + 4]
        (local / name).write_text(header + '\n```{r}\n' + body + '\n```\n')
    (tmp_path / 'traits.csv').write_text('species,family,score\nS1,F,1\nS2,F,2\n')
    (tmp_path / 'tree.nwk').write_text('(S1:1,S2:1);\n')
    (tmp_path / 'prune.txt').write_text('S1\n')
    (tmp_path / 'secondary').mkdir()
    (tmp_path / 'secondary/prune.txt').write_text('S2\n')
    (tmp_path / 'main.nf').write_text(
        f"include {{ REPORTING }} from '{ROOT / 'workflows/reporting.nf'}'\n"
        'workflow { REPORTING(Channel.value(file(params.trait_file)),Channel.value(file(params.tree_file))) }\n'
    )
    (tmp_path / 'audit.config').write_text("""
process.executor='local'
process.shell=['/bin/bash','-euo','pipefail']
params {
 trait_file='traits.csv'; tree_file='tree.nwk'; outdir='out'; prune_data=false
 use_singularity=false; use_apptainer=false; seed='1998'; clade_name='primates'
 taxon_of_interest='family'; sp_colname='species'; traitname='score'; n_trait=''; c_trait=''
 branch_trait=''; secondary_trait=''; pss_top_pct='0.05'; perm_strategy='best_model'
 trait_type='continuous'; max_contrasts='0'; prune_list='prune.txt'; prune_list_secondary=''
}
""")
    return tmp_path


def run_project(project, *, prune=False, secondary=False, resume=False):
    args = ['nextflow', '-log', 'nextflow.log', 'run', 'main.nf', '-c', 'audit.config',
            '-w', 'work', '-ansi-log', 'false', '--prune_data', str(prune).lower()]
    if secondary:
        args += ['--prune_list_secondary', 'secondary/prune.txt']
    if resume:
        args += ['-resume']
    result = subprocess.run(args, cwd=project, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    tasks = {}
    for command in (project / 'work').glob('*/*/.command.sh'):
        script = command.read_text()
        for stage, filename in [('prune', '0.Data_pruning.Rmd'),
                                ('dataset', '1.Dataset_exploration.Rmd'),
                                ('phenotype', '2.Phenotype_exploration.Rmd')]:
            if filename in script:
                tasks.setdefault(stage, []).append(command.parent)
    return tasks


@pytest.mark.parametrize(('prune', 'secondary'), [(False, False), (True, False), (True, True)])
def test_result_directories_are_owned_by_their_tasks(reporting_project, prune, secondary):
    tasks = run_project(reporting_project, prune=prune, secondary=secondary)
    for stage in ['dataset', 'phenotype'] + (['prune'] if prune else []):
        output = tasks[stage][0] / 'data_exploration'
        assert output.is_dir() and not output.is_symlink()
        assert (output / 'phenotype_marker.txt').exists() == (stage == 'phenotype')
    published = reporting_project / 'out/data_exploration'
    assert (published / 'phenotype_marker.txt').read_text().strip() == 'phenotype'
    archive = published / '1.Data-exploration/1.Species_distribution/original_trait_file.tsv'
    assert '\t' in archive.read_text().splitlines()[0]
    if prune:
        assert (published / 'prune_marker.txt').read_text().strip() == 'S1'
        script = (tasks['prune'][0] / '.command.sh').read_text()
        assert 'prune_primary/prune.txt' in script
        if secondary:
            assert 'prune_secondary/prune.txt' in script


def test_resume_invalidates_changed_prune_list(reporting_project):
    first = run_project(reporting_project, prune=True)
    (reporting_project / 'prune.txt').write_text('S2\n')
    second = run_project(reporting_project, prune=True, resume=True)
    assert len(first['prune']) == 1 and len(second['prune']) == 2
    assert (reporting_project / 'out/data_exploration/prune_marker.txt').read_text().strip() == 'S2'
