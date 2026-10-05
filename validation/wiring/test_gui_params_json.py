"""The params.json a generated run writes: valid JSON for every saved template, carrying what the GUI fields hold.

The batch script exports the model values as shell variables and the single-run script expands them inside a heredoc. The
test runs both pieces of the rendered scripts with bash, so a Jinja condition, a quote or a missing comma shows up as
invalid JSON, and a field that never reaches the file as a missing or wrong value.
"""
import glob
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.project_io import load_project  # noqa: E402

TEMPLATES = sorted(glob.glob(str(ROOT / "gui/templates/*.json")))


def params_json(project):
    """The params.json that the rendered scripts write: the batch exports, then the single-run heredoc, run by bash."""
    batch, single = render_batch(project), render_single(project)
    exports = [l.strip() for l in batch.splitlines() if re.match(r'^\s*export \w+=', l)]
    start = single.index('cat > "$PARAMS_JSON" <<PARAMS_EOF')
    body = single[single.index("\n", start) + 1:single.index("\nPARAMS_EOF", start)]
    unquoted = set(re.findall(r'^\s*(?:\{%.*?%\})*\s*"\w+"\s*:\s*\$\{?(\w+)\}?\s*,?\s*$', body, re.M))  # boolean switches the script computes
    script = "\n".join([f"{v}=false" for v in sorted(unquoted)] + exports + ["cat <<PARAMS_EOF", body, "PARAMS_EOF"])
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True).stdout
    return json.loads(out)


@pytest.mark.parametrize("template", TEMPLATES, ids=[Path(t).name for t in TEMPLATES])
def test_a_saved_template_writes_valid_json(template):
    d = params_json(load_project(Path(template)))
    assert d["caas_evidence_top_n"] == "0" and d["fcs_fdr_permsum"] == "0.05"
    assert "fcs_fdr_wilcoxon" not in d and "fcs_fdr_lachenbruch" not in d      # blank follows fcs_fdr: the key is not written


def _project():
    return load_project(Path(ROOT / "gui/templates/cancer_no_prune_multi.json"))


def test_the_values_of_the_new_fields_reach_the_file():
    p = _project()
    p.modules.caas.caas_perms_postproc = False
    p.modules.scoring.scoring_gene_perm_pooled = True
    p.modules.scoring.scoring_hypotheses_pairs = "/x/pairs_scoring.tsv"
    p.modules.disambiguation.ct_disambig_hypotheses_pairs = "/x/pairs_obs.tsv"
    p.modules.enrichment.fcs_fdr_wilcoxon = "0.2"
    p.modules.enrichment.fcs_fdr_lachenbruch = "0.3"
    p.modules.enrichment.fcs_fdr_permsum = "0.07"
    p.modules.enrichment.pfam_cache_dir = "/x/pfam"
    d = params_json(p)
    assert d["caas_perms_postproc"] is False and d["scoring_gene_perm_pooled"] is True
    assert (d["scoring_hypotheses_pairs"], d["ct_disambig_hypotheses_pairs"]) == ("/x/pairs_scoring.tsv", "/x/pairs_obs.tsv")
    assert (d["fcs_fdr_wilcoxon"], d["fcs_fdr_lachenbruch"], d["fcs_fdr_permsum"]) == ("0.2", "0.3", "0.07")
    assert d["pfam_cache_dir"] == "/x/pfam"


def test_the_defaults_of_the_new_fields_are_the_conf_defaults():
    d = params_json(_project())
    assert d["caas_perms_postproc"] is True and d["scoring_gene_perm_pooled"] is False and d["caas_permulation_enrichment"] is True
    assert d["scoring_hypotheses_pairs"] == "" and d["ct_disambig_hypotheses_pairs"] == "" and d["pfam_cache_dir"] == ""


def test_one_blank_fdr_gate_does_not_hide_the_other():
    p = _project()
    p.modules.enrichment.fcs_fdr_lachenbruch = "0.3"
    d = params_json(p)
    assert d["fcs_fdr_lachenbruch"] == "0.3" and "fcs_fdr_wilcoxon" not in d


def test_every_key_the_audit_sees_in_the_template_is_in_the_file():
    sys.path.insert(0, str(ROOT / "validation/wiring"))
    import param_audit
    keys = set(param_audit.gui_params(ROOT))
    written = set(params_json(_project()))
    # keys behind a Jinja condition are absent when their field is blank
    assert keys - written <= {"fcs_fdr_wilcoxon", "fcs_fdr_lachenbruch"}, sorted(keys - written)
    assert not written - keys
