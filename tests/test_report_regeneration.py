"""Independent report regeneration: scripts the GUI builds must match the migrated Rmd contracts.

detect_reports()/build_script() are exercised on a synthetic outdir that carries the
HTML of every migrated report (11/5/12/14/POSENRICH modular/13/15) plus the published
artefacts they read. Each generated script must (a) parse as bash, (b) pass only
params the target Rmd declares in its YAML frontmatter, and (c) stage the nominal
active gene lists (active_*.tsv) for AMI, never the retired slice_*.tsv.
"""
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

RMD_OF = {
    "scoring": "subworkflows/SCORING/local/11.Scoring_report.Rmd",
    "rer": "subworkflows/RERCONVERGE/local/5.RERconverge_report.Rmd",
    "fcs_scoring": "subworkflows/ENRICHMENT/local/12.FCS_general_report.Rmd",
    "posenrich": "subworkflows/ENRICHMENT/local/14.Position_enrichment_report.Rmd",
    "posenrich_modular": "subworkflows/ENRICHMENT/local/posenrich_modular_report.Rmd",
    "ami": "subworkflows/ENRICHMENT/local/13.AMI_analysis.Rmd",
    "compare": "subworkflows/ENRICHMENT/local/15.Comparison_report.Rmd",
}


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("fixture\n")
    return path


def _make_outdir(base: Path) -> Path:
    out = base / "outdir"
    reports = out / "html_reports"
    _touch(reports / "11.Scoring_report_trait.html")
    _touch(reports / "5.RERconverge_report.html")
    _touch(reports / "12.FCS_scoring_trait.html")
    _touch(reports / "14.Position_enrichment_report_trait.html")
    _touch(reports / "POSENRICH_modular_report_trait.html")
    _touch(reports / "13.AMI_analysis_trait.html")
    _touch(reports / "15.Comparison_report_trait.html")
    _touch(out / "scoring" / "position_scores.tsv")
    _touch(out / "scoring" / "gene_scores.tsv")
    _touch(out / "scoring" / "gene_correlations.tsv")
    _touch(out / "scoring" / "gene_lists" / "active_global.tsv")
    _touch(out / "scoring" / "gene_lists" / "active_top.tsv")
    _touch(out / "scoring" / "gene_lists" / "active_bottom.tsv")
    _touch(out / "cleaned_background.txt")
    _touch(out / "rerconverge" / "rer_results" / "trait.continuous.output")
    _touch(out / "fcs" / "fcs_results" / "fcs_enrich_input.tsv")
    _touch(out / "fcs" / "fcs_results" / "fcs_report_metadata.json")
    _touch(out / "fcs_stats.tsv")
    _touch(out / "posenrich" / "posenrich_characterization.tsv")
    _touch(out / "posenrich" / "posenrich_leading_edge.tsv")
    for sub in ("registry", "characterization", "localization", "global", "proximity"):
        (out / "posenrich" / sub).mkdir(parents=True, exist_ok=True)
    _touch(out / "ami" / "network.sif")
    (out / "ami" / "domino_modules").mkdir(parents=True, exist_ok=True)
    _touch(out / "ami" / "network_edge_scores.tsv")
    return out


def _declared_params(rmd: Path) -> set[str]:
    declared = set()
    in_params = False
    for line in rmd.read_text().splitlines():
        if re.match(r"^params:\s*$", line):
            in_params = True
            continue
        if in_params:
            m = re.match(r"^  ([A-Za-z_][A-Za-z0-9_.]*):", line)
            if m:
                declared.add(m.group(1))
            elif line and not line.startswith((" ", "\t", "#")):
                break
    return declared


def _rendered_params(script: str) -> set[str]:
    m = re.search(r"params = list\((.*?)\n\s*\)", script, re.DOTALL)
    assert m, "no params = list(...) block in the generated script"
    return set(re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_.]*)\s*=", m.group(1), re.MULTILINE))


def test_regeneration_scripts_match_rmd_contracts(tmp_path):
    from gui.generation.report_registry import detect_reports, build_script

    out = _make_outdir(tmp_path)
    found = {r.id: r for r in detect_reports(str(out), ROOT)}
    assert set(RMD_OF).issubset(found), f"missing detected reports: {set(RMD_OF) - set(found)}"

    for spec_id, rmd_rel in RMD_OF.items():
        report = found[spec_id]
        assert report.runnable, f"{spec_id}: {report.missing_required}"
        script = build_script(report, ROOT, use_singularity=False)
        check = subprocess.run(["bash", "-n"], input=script, text=True, capture_output=True)
        assert check.returncode == 0, f"{spec_id}: {check.stderr}"
        declared = _declared_params(ROOT / rmd_rel)
        assert declared, f"{spec_id}: no params parsed from {rmd_rel}"
        unknown = _rendered_params(script) - declared
        assert not unknown, f"{spec_id}: params not declared in {rmd_rel}: {sorted(unknown)}"

    ami_script = build_script(found["ami"], ROOT, use_singularity=False)
    assert "active_*.tsv" in ami_script and 'name="${b#active_}"' in ami_script
    assert "slice_" not in ami_script
    assert "string_cache_dir" in ami_script
    assert "percentile" not in ami_script.lower()
    for spec_id in ("posenrich", "compare"):
        script = build_script(found[spec_id], ROOT, use_singularity=False)
        assert "position_lists_dir" not in script and "gene_lists_dir" not in script
    assert "gene_lists_dir" not in build_script(found["fcs_scoring"], ROOT, use_singularity=False)


def test_fcs_regeneration_requires_original_settings(tmp_path):
    from gui.generation.report_registry import detect_reports

    out = _make_outdir(tmp_path)
    (out / "fcs" / "fcs_results" / "fcs_report_metadata.json").unlink()
    found = {r.id: r for r in detect_reports(str(out), ROOT)}
    assert not found["fcs_scoring"].runnable
    assert any("Saved FCS settings" in label for label in found["fcs_scoring"].missing_required)
