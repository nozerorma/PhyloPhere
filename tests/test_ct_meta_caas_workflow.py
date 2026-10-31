"""CT_META_CAAS process contracts, safe parameter transport and resume caching."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def project(tmp_path):
    if not shutil.which('nextflow') or not shutil.which('Rscript'):
        pytest.skip('Run with conda run -n phylophere')
    p = tmp_path / 'project with spaces'
    local = p / 'subworkflows/CT_META_CAAS/local'
    local.mkdir(parents=True)
    for rel in ['workflows/ct_meta_caas.nf', 'subworkflows/CT_META_CAAS/ctpp_meta_caas.nf',
                'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R',
                'subworkflows/ENRICHMENT/local/src/report_presentation.R']:
        dest = p / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, dest)
    for name in ['7.CAAS_pattern_annotation.Rmd', '15.CAAS_significance_report.Rmd']:
        source = (ROOT / 'subworkflows/CT_META_CAAS/local' / name).read_text()
        header = source[:source.index('\n---', 4) + 4]
        body = "source('src/report_style.R');set_report_style()\n"
        body += "stopifnot(file.exists(params$" + ('discovery_input' if name.startswith('7.') else 'position_scores_input') + "))\n"
        if name.startswith('15.'):
            body += "dir.create('significance');writeLines(params$traitname,'significance/meta_caas_significance.tsv')\n"
        (local / name).write_text(header + '\n```{r}\n' + body + 'cat(params$traitname)\n```\n')
    for name, text in [("discovery's input.tab", 'gene\tposition\tpattern\ng1\t0\t1\n'),
                       ('empty.tab', 'gene\tposition\tpattern\n'), ('background.txt', 'g1\n')]:
        (p / name).write_text(text)
    (p / 'main.nf').write_text("""
include { CT_META_CAAS; } from './workflows/ct_meta_caas.nf'
include { CAAS_SIGNIFICANCE_REPORT } from './subworkflows/CT_META_CAAS/ctpp_meta_caas.nf'
process OTHER {
 output:
 path 'other.txt'
 script:
 "echo completed > other.txt"
}
workflow {
 def d = file(params.discovery_from)
 if (params.significance) {
  CAAS_SIGNIFICANCE_REPORT(d,d,d,d,d,d)
 } else {
  CT_META_CAAS(params.live ? Channel.value(d) : null, params.live ? Channel.value(file(params.background_input)) : null)
 }
 OTHER()
}
""")
    (p / 'audit.config').write_text("""
process.executor='local'
params {
 discovery_from="discovery's input.tab"; background_input='background.txt'
 outdir='out'; use_singularity=false; use_apptainer=false
 traitname='test'; seed=1998; live=false; significance=false
}
""")
    return p


def run(p, *extra, success=True):
    result = subprocess.run(['nextflow', '-log', 'nextflow.log', 'run', 'main.nf',
                             '-c', 'audit.config', '-w', 'work', '-ansi-log', 'false', *extra],
                            cwd=p, capture_output=True, text=True, timeout=120)
    output = result.stdout + result.stderr
    assert (result.returncode == 0) == success, output
    return output


def tasks(p, marker):
    return [f.parent for f in (p / 'work').glob('*/*/.command.sh') if marker in f.read_text()]


@pytest.mark.parametrize('live', [False, True])
def test_empty_discovery_renders_and_other_work_completes(project, live):
    run(project, '--discovery_from', 'empty.tab', '--live', str(live).lower())
    assert (project / 'out/html_reports/7.CAAS_pattern_annotation.html').exists()
    assert len(tasks(project, 'echo completed')) == 1
    assert (tasks(project, 'echo completed')[0] / 'other.txt').exists()


@pytest.mark.parametrize('n', [0, 2])
def test_background_directory_requires_one_file(project, n):
    bg = project / 'backgrounds'; bg.mkdir()
    for k in range(n):
        (bg / f'{k}_background_genes.output').write_text('g1\n')
    (bg / 'directory_background_genes').mkdir()  # directories are not candidates
    output = run(project, '--background_input', 'backgrounds', success=False)
    assert f'exactly one background_genes file' in output and f'found {n}' in output


@pytest.mark.parametrize('significance', [False, True])
def test_single_background_directory_and_source_resume_invalidation(project, significance):
    bg = project / 'backgrounds'; bg.mkdir()
    (bg / 'background_genes.output').write_text('g1\n')
    args = ['--background_input', 'backgrounds', '--significance', str(significance).lower()]
    run(project, *args)
    stem = '15.CAAS_significance_report' if significance else '7.CAAS_pattern_annotation'
    marker = f"rmarkdown::render('{stem}"
    rmd = project / f'subworkflows/CT_META_CAAS/local/{stem}.Rmd'
    rmd.write_text(rmd.read_text() + '\n<!-- changed report -->\n')
    run(project, *args, '-resume')
    assert len(tasks(project, marker)) == 2
    style = project / 'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'
    style.write_text(style.read_text() + '\n# changed style\n')
    run(project, *args, '-resume')
    assert len(tasks(project, marker)) == 3
    if significance:
        helper = project / 'subworkflows/ENRICHMENT/local/src/report_presentation.R'
        helper.write_text(helper.read_text() + '\n# changed table helper\n')
        run(project, *args, '-resume')
        assert len(tasks(project, marker)) == 4


@pytest.mark.parametrize('significance', [False, True])
def test_parameters_are_literal_data_and_outputs_publish(project, significance):
    trait = "O'Brien $(touch INJECTED) `touch INJECTED2`\nPHYLOPHERE_REPORT_PARAMS"
    run(project, '--traitname', trait, '--significance', str(significance).lower())
    rendered = tasks(project, "rmarkdown::render('")
    assert len(rendered) == 1
    assert json.loads((rendered[0] / 'report_params.json').read_text())['traitname'] == trait
    assert not (rendered[0] / 'INJECTED').exists()
    assert not (rendered[0] / 'INJECTED2').exists()
    if significance:
        assert (project / 'out/meta_caas/significance/meta_caas_significance.tsv').read_text().rstrip('\n') == trait
        assert (project / 'out/html_reports/15.CAAS_significance_report.html').exists()
