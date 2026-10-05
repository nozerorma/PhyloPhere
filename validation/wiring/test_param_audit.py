"""Every surface of a pipeline parameter agrees (param_audit.py), and the audit itself detects the defects it claims to.

The first group runs the audit on the tree under test: the only findings left are the ones the allowlist justifies.
The second group copies the tracked files to a scratch directory, injects one known defect at a time and checks that the
audit reports it, so a silent audit cannot pass for a clean tree.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
sys.path.insert(0, str(HERE))
import param_audit as pa  # noqa: E402

needs_nextflow = pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")
ALLOW = HERE / "param_audit_allow.tsv"


@pytest.fixture(scope="module")
def findings():
    return pa.audit(ROOT)


@needs_nextflow
def test_every_finding_is_justified_in_the_allowlist(findings):
    allowed = pa.load_allowlist(ALLOW)
    open_ = [f for f in findings if (f[0], f[1]) not in allowed]
    assert not open_, "\n".join("\t".join(f) for f in open_)


@needs_nextflow
def test_the_allowlist_has_no_entry_the_audit_no_longer_reports(findings):
    seen = {(f[0], f[1]) for f in findings}
    assert not [k for k in pa.load_allowlist(ALLOW) if k not in seen]


def test_every_allowlist_entry_says_why():
    assert all(reason.strip() for reason in pa.load_allowlist(ALLOW).values())


# ── the audit detects what it claims to ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def copy(tmp_path_factory):
    """The tracked files of the tree under test, minus the large ones, in a scratch directory."""
    dst = tmp_path_factory.mktemp("audit_copy")
    for name in pa._tracked(ROOT):
        src = ROOT / name
        if name.startswith(("validation/tier1/input", ".claude")) or not src.is_file() or src.stat().st_size > 2_000_000:
            continue
        (dst / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst / name)
    subprocess.run(["git", "init", "-q"], cwd=dst, check=True)
    subprocess.run(["git", "add", "-A"], cwd=dst, check=True, capture_output=True)
    return dst


def _edited(copy, path, old, new, count=1):
    """Context-free helper: apply one replacement, return a restore function."""
    f = copy / path
    text = f.read_text()
    assert text.count(old) == count, (path, old, text.count(old))
    f.write_text(text.replace(old, new))
    return lambda: f.write_text(text)


def _codes(copy):
    return {(c, n) for c, n, _ in pa.audit(copy)}


MUTATIONS = [
    ("a parameter declared and read by nobody", "conf/scoring.config", "    scoring_p_emp_thr ", "    zz_dead_param = 1\n    scoring_p_emp_thr ",
     ("DECLARED_UNREAD", "zz_dead_param")),
    ("a parameter read and declared nowhere", "main.nf", "def run_caas_permulation ", "def zz_probe = params.zz_undeclared_param\n        def run_caas_permulation ",
     ("READ_UNDECLARED", "zz_undeclared_param")),
    ("a parameter the GUI writes and nothing declares", "gui/generation/templates/run_single.sh.j2", '  "caas_evidence_top_n":',
     '  "zz_gui_only": "x",\n  "caas_evidence_top_n":', ("GUI_UNDECLARED", "zz_gui_only")),
    ("a model field no template renders", "gui/models/modules.py", '    caas_evidence_top_n: str = "0"',
     '    zz_unrendered: str = ""\n    caas_evidence_top_n: str = "0"', ("FIELD_NOT_RENDERED", "modules.scoring.zz_unrendered")),
    ("a model field without a tab field", "gui/models/modules.py", '    caas_evidence_top_n: str = "0"',
     '    zz_hidden: str = ""\n    caas_evidence_top_n: str = "0"', ("FIELD_NOT_IN_TAB", "modules.scoring.zz_hidden")),
    ("a default that differs between the model and conf", "gui/models/modules.py", '    scoring_p_emp_thr: str = "0.05"',
     '    scoring_p_emp_thr: str = "0.07"', ("DEFAULT_MISMATCH", "scoring_p_emp_thr")),
    ("a tab field naming nothing in the model", "gui/widgets/tabs/scoring_tab.py", '        Section("Evidence of the best positions"),',
     '        FieldSpec(name="zz_ghost", label="ghost"),\n        Section("Evidence of the best positions"),', ("TAB_FIELD_NOT_IN_MODEL", "modules.scoring.zz_ghost")),
    ("a blank that replaces a non-empty conf default", "gui/generation/templates/run_single.sh.j2",
     '"lg_dat_path": "${LG_DAT_PATH:-$REPO_DIR/subworkflows/SELECTION/local/dat/${_FADE_MODEL_LOWER}.dat}",', '"lg_dat_path": "${LG_DAT_PATH:-}",',
     ("BLANK_BURNS_DEFAULT", "lg_dat_path")),
    ("a parameter the GUI stops writing", "gui/generation/templates/run_single.sh.j2", '  "scoring_gene_perm_pooled": {{ run_defaults.scoring_gene_perm_pooled }},\n', "",
     ("NOT_IN_GUI", "scoring_gene_perm_pooled")),
    ("a template key the model no longer has", "gui/templates/cancer_no_prune_multi.json", '"scoring_gene_top_pct"', '"zz_old_key": 1, "scoring_gene_top_pct"',
     ("TEMPLATE_STALE_KEY", "cancer_no_prune_multi.json: modules.scoring.zz_old_key")),
]


@needs_nextflow
@pytest.mark.parametrize("what,path,old,new,expected", MUTATIONS, ids=[m[0] for m in MUTATIONS])
def test_the_audit_reports_an_injected_defect(copy, what, path, old, new, expected):
    base = _codes(copy)
    assert expected not in base
    restore = _edited(copy, path, old, new)
    try:
        assert expected in _codes(copy)
    finally:
        restore()
