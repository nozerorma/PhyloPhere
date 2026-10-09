#!/usr/bin/env python3
# report_registry.py — Detection of the HTML reports of a run and standalone scripts that re-render them.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Detects which of the pipeline's Rmd reports exist in an output directory
(detect_reports) and builds, per report, a shell script that re-renders it from the
files already on disk there. Each script is modeled on the rmarkdown::render() call
of the Nextflow process that produced the report (subworkflows/*/*.nf), with outdir
paths in place of the staged inputs and with the render parameters fixed in this
module. The module never runs R itself.

Each report is a ReportSpec: a filename regex that identifies it, a find_slots
function that locates its input files, and a build_script function. The inputs of a
report are searched with a short ordered list of globs (most specific first) rooted
at outdir, which follow the subdirectory each upstream process publishes into (the
publishDir of its .nf file). A glob miss leaves the slot unresolved and never raises.
An unresolved optional slot is passed to R as NULL, as the Nextflow processes do for
their NO_* sentinels, so the report skips that section. An unresolved required slot
makes the report non-runnable until the user selects the file by hand; some inputs
(the --trait_file and --tree_file of TRAIT_ANALYSIS) are pipeline inputs rather than
outputs and are rarely found automatically.

Called by:  gui/widgets/common/regenerate_dialog.py (detect_reports, build_script);
            the generated scripts go to MainWindow._show_preview, the same
            preview/save window as "Generate Scripts"
Inputs:     the output directory of a run (local, or a remote listing from gui.remote)
Outputs:    one bash script per report (regenerate_<report id>.sh) that renders the
            report into a scratch directory and copies the HTML to the directories its
            Nextflow process publishes to
"""

from __future__ import annotations

import fnmatch
import re
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Callable, Optional


# ── Listing: local-or-remote filesystem view ──────────────────────────────────


class Listing:
    """Thin directory view (glob, is_dir, "/" chaining, as in pathlib) that the
    find_slots functions are written against, so they work unchanged on a local
    outdir and on a remote one reached over SSH. Remote mode makes no network call
    here: it is given a pre-fetched {relpath: is_dir} map (gui.remote.list_all_remote,
    one round trip made by the caller). Results are PurePosixPath, correct for a
    pipeline that runs on a Linux cluster whatever the host OS of the GUI.
    """

    def __init__(self, root: str, remote_entries: Optional[dict[str, bool]] = None, prefix: str = ""):
        self.root = PurePosixPath(root)
        self._remote = remote_entries  # None => local filesystem; else {relpath: is_dir}
        self._prefix = prefix  # relative to root

    @property
    def path(self) -> PurePosixPath:
        return self.root / self._prefix if self._prefix else self.root

    def __truediv__(self, other: str) -> "Listing":
        new_prefix = f"{self._prefix}/{other}" if self._prefix else str(other)
        return Listing(str(self.root), self._remote, new_prefix)

    def __str__(self) -> str:
        return str(self.path)

    def glob(self, pattern: str) -> list[PurePosixPath]:
        full_pattern = f"{self._prefix}/{pattern}" if self._prefix else pattern
        if self._remote is None:
            return sorted(PurePosixPath(p) for p in Path(str(self.root)).glob(full_pattern))
        # fnmatch matches the whole relative path and its '*' crosses '/', so
        # '**/x.tsv' patterns behave as in the local glob() branch above.
        hits = sorted(rel for rel in self._remote if fnmatch.fnmatch(rel, full_pattern))
        return [self.root / rel for rel in hits]

    def is_dir(self) -> bool:
        if self._remote is None:
            return Path(str(self.path)).is_dir()
        return self._remote.get(self._prefix) is True


# ── Slot / detection data model ───────────────────────────────────────────────


@dataclass
class InputSlot:
    key: str  # unique within one report's slot list; used as the override dict key
    label: str  # shown in the regenerate dialog
    required: bool
    path: Optional[PurePosixPath] = None  # resolved (glob hit) or user-overridden; None = unresolved
    note: str = ""  # short UI hint, e.g. "optional: passed as NULL if not set"


@dataclass
class DetectedReport:
    id: str
    display_name: str
    html_path: PurePosixPath  # the html_reports/*.html file that triggered detection
    rmd_filename: str
    local_dir: str  # relative to repo_dir, e.g. "subworkflows/SCORING/local"
    traitname: str
    slots: list[InputSlot]
    spec: "ReportSpec" = field(repr=False)
    script_filename: str = ""

    @property
    def runnable(self) -> bool:
        return all(s.path is not None for s in self.slots if s.required)

    @property
    def missing_required(self) -> list[str]:
        return [s.label for s in self.slots if s.required and s.path is None]


@dataclass
class ReportSpec:
    id: str
    display_name: str
    html_regex: re.Pattern
    find_slots: Callable[[Listing], list[InputSlot]]
    build_script: Callable[["DetectedReport", Path, bool], str]
    guess_traitname: Callable[[str], str] = lambda basename: "unknown_trait"


# ── Path-search helpers ───────────────────────────────────────────────────────


def _first_match(listing: Listing, *rel_globs: str) -> Optional[PurePosixPath]:
    """Return the first hit of the globs tried in order, most specific first, or None."""
    for pattern in rel_globs:
        hits = listing.glob(pattern)
        if hits:
            return hits[0]
    return None


def _slot(key: str, label: str, required: bool, path: Optional[PurePosixPath], note: str = "") -> InputSlot:
    if not note and not required:
        note = "optional — passed as NULL if not set"
    return InputSlot(key=key, label=label, required=required, path=path, note=note)


_TRAIT_SUFFIX_RE = re.compile(r"_(?P<trait>[A-Za-z0-9][A-Za-z0-9_.-]*)\.html$")


def _guess_trait_from_suffix(basename: str, prefix_to_strip: str) -> str:
    """Trait name from a report filename: '11.Scoring_report_um_myotis.html' → 'um_myotis'.

    Falls back to 'unknown_trait', the default several Nextflow report processes use
    when --traitname is empty, if the name carries no trait suffix.
    """
    stem = basename
    if stem.startswith(prefix_to_strip):
        stem = stem[len(prefix_to_strip):]
    m = _TRAIT_SUFFIX_RE.search("_" + stem) if not stem.startswith("_") else _TRAIT_SUFFIX_RE.search(stem)
    return m.group("trait") if m else "unknown_trait"


# ── Script-assembly helpers ───────────────────────────────────────────────────


def _r_arg(path: Optional[Path]) -> str:
    return f"'{path}'" if path is not None else "NULL"


def _entrypoint_prefix(use_singularity: bool) -> str:
    return "/usr/local/bin/_entrypoint.sh " if use_singularity else ""


def _render_call(rmd_filename: str, param_lines: list[str], output_file: str, extra_arg: str = "") -> str:
    """Build the rmarkdown::render() call of one report as R source text.

    extra_arg is an additional top-level render() argument (e.g. "envir = new.env()")
    and must be passed here, so that it is written inside the call, next to
    output_file. Appending it to the returned string would put it after the closing
    parenthesis, a stray top-level comma that R's parser rejects.
    """
    params = ",\n                ".join(param_lines)
    extra = f",\n            {extra_arg}" if extra_arg else ""
    return (
        f"rmarkdown::render(\n"
        f"            '{rmd_filename}',\n"
        f"            params = list(\n"
        f"                {params}\n"
        f"            ),\n"
        f"            output_file = '{output_file}'{extra}\n"
        f"        )"
    )


def _wrap_script(
    *,
    header_comment: str,
    repo_dir: Path,
    local_dir: str,
    pre_lines: list[str],
    render_block: str,
    output_file: str,
    publish_targets: list[Path],
    use_singularity: bool,
    post_lines: Optional[list[str]] = None,
) -> str:
    """Assemble one standalone bash script that re-renders a report.

    The script stages the report's local/ directory into a scratch work directory,
    runs pre_lines (input staging), renders, then copies the HTML to every
    publish_targets directory, keeping a .bak of any file it overwrites.
    post_lines run after a successful render and before the publish copy; some
    upstream processes also copy files after rendering (the independent-contrasts
    report is one).
    """
    entry = _entrypoint_prefix(use_singularity)
    pre = "\n".join(pre_lines)
    post = "\n".join(post_lines or [])
    copy_lines = []
    for target in publish_targets:
        copy_lines.append(f'mkdir -p "{target}"')
        copy_lines.append(
            f'[ -f "{target}/{output_file}" ] && cp "{target}/{output_file}" "{target}/{output_file}.bak" || true'
        )
        copy_lines.append(f'cp "{output_file}" "{target}/"')
    copy_block = "\n".join(copy_lines)

    return f"""#!/usr/bin/env bash
# {header_comment}
# Standalone regeneration — re-renders this report from files already produced
# in the chosen output directory. Does not rerun any pipeline computation.
set -euo pipefail

REPO_DIR="{repo_dir}"
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

cp -R "$REPO_DIR/{local_dir}"/* .
{pre}

{entry}Rscript -e "
    {render_block}
"

if [ ! -f "{output_file}" ]; then
    echo "ERROR: {output_file} was not produced." >&2
    exit 1
fi
{post}

{copy_block}
echo "Regenerated {output_file}"
"""


# ── SCORING ───────────────────────────────────────────────────────────────────


def _scoring_find_slots(outdir: Listing) -> list[InputSlot]:
    scoring = outdir / "scoring"
    fade_top = outdir / "selection" / "fade" / "top"
    fade_bot = outdir / "selection" / "fade" / "bottom"
    return [
        _slot("position_scores", "Position scores TSV", True, _first_match(scoring, "position_scores.tsv")),
        _slot("gene_scores", "Gene scores TSV", True, _first_match(scoring, "gene_scores.tsv")),
        _slot("gene_correlations", "Gene correlations TSV", True, _first_match(scoring, "gene_correlations.tsv")),
        _slot("stress_summary", "Stress summary TSV", False,
              _first_match(scoring, "position_score_stress_summary.tsv")),
        _slot("stress_correlations", "Stress correlations TSV", False,
              _first_match(scoring, "position_score_stress_correlations.tsv")),
        _slot("stress_rank", "Stress rank agreement TSV", False,
              _first_match(scoring, "position_score_stress_rank_agreement.tsv")),
        _slot("stress_overlap", "Stress top-overlap TSV", False,
              _first_match(scoring, "position_score_stress_top_overlap.tsv")),
        _slot("stress_variants", "Stress variants TSV", False,
              _first_match(scoring, "position_score_stress_variants.tsv")),
        _slot("fade_site_top", "FADE per-site BF (top)", False,
              _first_match(fade_top, "fade_site_bf_top.tsv")),
        _slot("fade_site_bot", "FADE per-site BF (bottom)", False,
              _first_match(fade_bot, "fade_site_bf_bottom.tsv")),
        # The genomic-coordinates file is params.gene_ensembl_file (resolved from
        # there in scoring.nf), the same file as the gene Ensembl slot of CT_POSTPROC.
        _slot("genomic_info", "Gene genomic coordinates TSV", False,
              _first_match(outdir, "postproc/postproc_inputs/*", "**/*genomic_info*.tsv", "**/*genomic*coords*.tsv")),
        _slot("caas_perms", "CAAS permulation RDS", False,
              _first_match(outdir, "**/*caas*perm*.rds")),
        _slot("caas_pos_sample", "CAAS position sample TSV", False,
              _first_match(outdir, "caas_permulation/perm_pos_sample.tsv", "**/*pos*sample*.tsv")),
        _slot("caas_pos_quantiles", "CAAS position quantiles TSV", False,
              _first_match(outdir, "caas_permulation/perm_pos_quantiles.tsv", "**/*pos*quantiles*.tsv")),
        _slot("filtered_discovery", "Filtered discovery TSV", False,
              _first_match(outdir, "**/filtered_discovery.tsv")),
        _slot("background_file", "Background gene list", False,
              _first_match(outdir, "**/cleaned_background*.txt", "**/*background*.txt")),
    ]


def _scoring_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = f"11.Scoring_report_{report.traitname}.html"
    render = _render_call(
        "11.Scoring_report.Rmd",
        [
            f"position_scores_file = '{s['position_scores']}'",
            f"gene_scores_file     = '{s['gene_scores']}'",
            f"gene_corr_file       = '{s['gene_correlations']}'",
            f"traitname            = '{report.traitname}'",
            "output_dir           = '.'",
            "top_pct              = 0.10",
            "gene_top_pct         = 0.10",
            f"stress_summary_file  = {_r_arg(s['stress_summary'])}",
            f"stress_corr_file     = {_r_arg(s['stress_correlations'])}",
            f"stress_rank_file     = {_r_arg(s['stress_rank'])}",
            f"stress_overlap_file  = {_r_arg(s['stress_overlap'])}",
            f"stress_variants_file = {_r_arg(s['stress_variants'])}",
            f"fade_site_top_file   = {_r_arg(s['fade_site_top'])}",
            f"fade_site_bot_file   = {_r_arg(s['fade_site_bot'])}",
            f"genomic_info_file    = {_r_arg(s['genomic_info'])}",
            f"caas_perms_file      = {_r_arg(s['caas_perms'])}",
            f"caas_pos_sample_file = {_r_arg(s['caas_pos_sample'])}",
            f"caas_pos_quantiles_file = {_r_arg(s['caas_pos_quantiles'])}",
            f"filtered_discovery_file = {_r_arg(s['filtered_discovery'])}",
            f"background_file      = {_r_arg(s['background_file'])}",
            "window_size_bp       = 1000000",
            "direction            = 'combined'",
            "seed                 = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="SCORING — 11.Scoring_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/SCORING/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "scoring", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── RERconverge ───────────────────────────────────────────────────────────────


def _rer_find_slots(outdir: Listing) -> list[InputSlot]:
    rer_dir = outdir / "rerconverge" / "rer_results"
    return [
        _slot("continuous_rds", "RERconverge continuous output RDS", True,
              _first_match(rer_dir, "*.continuous.output")),
    ]


def _rer_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "5.RERconverge_report.html"
    render = _render_call(
        "5.RERconverge_report.Rmd",
        [
            f"continuous_rds    = '{s['continuous_rds']}'",
            f"traitname         = '{report.traitname}'",
            "pval_threshold    = 0.05",
            "top_n_labels      = 15",
            "output_dir        = '.'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="RERCONVERGE — 5.RERconverge_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/RERCONVERGE/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "rerconverge" / "rer_results", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── FADE (one spec instance per direction) ────────────────────────────────────


def _make_fade_spec(direction: str) -> ReportSpec:
    def find_slots(outdir: Listing) -> list[InputSlot]:
        # fade_run.nf publishes the *.FADE.json files in selection/fade/<direction>/json/,
        # one level below fade_dir.
        fade_dir = outdir / "selection" / "fade" / direction
        json_dir = fade_dir / "json"
        json_files = sorted(json_dir.glob("*.FADE.json"))
        return [
            _slot("json_dir", f"FADE JSON directory ({direction})", True,
                  json_dir if json_files else None,
                  note=f"{len(json_files)} *.FADE.json file(s) found" if json_files else "no *.FADE.json files found"),
            # selection_utils.nf publishes the <direction>_species.txt files in the
            # shared selection/species_sets/, not inside fade_dir.
            _slot("fg_list_file", "Foreground species list", False,
                  _first_match(outdir / "selection" / "species_sets", f"{direction}_species.txt")),
        ]

    def build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
        s = {slot.key: slot.path for slot in report.slots}
        output_file = f"6.FADE_report_{direction}.html"
        # Symlink instead of copy: the *.FADE.json files of a full run can total tens
        # of GB and the scratch directory is removed right after the render reads
        # them once, so copying would only double the disk use.
        pre_lines = [f'ln -sf "{s["json_dir"]}"/*.FADE.json . 2>/dev/null || true']
        render = _render_call(
            "6.FADE_report.Rmd",
            [
                "json_dir        = '.'",
                f"direction       = '{direction}'",
                f"traitname       = '{report.traitname}'",
                "bf_threshold    = 100",
                "min_genes_hmap  = 3",
                "output_dir      = '.'",
                f"fg_list_file    = {_r_arg(s['fg_list_file'])}",
            ],
            output_file,
        )
        return _wrap_script(
            header_comment=f"FADE ({direction}) — 6.FADE_report.Rmd",
            repo_dir=repo_dir,
            local_dir="subworkflows/FADE/local",
            pre_lines=pre_lines,
            render_block=render,
            output_file=output_file,
            publish_targets=[report.html_path.parent.parent / "selection" / "fade" / direction, report.html_path.parent],
            use_singularity=use_singularity,
        )

    return ReportSpec(
        id=f"fade_{direction}",
        display_name=f"FADE report ({direction})",
        html_regex=re.compile(rf"^6\.FADE_report_{direction}\.html$"),
        find_slots=find_slots,
        build_script=build_script,
    )


# ── CT_ACCUMULATION ───────────────────────────────────────────────────────────


def _accumulation_find_slots(outdir: Listing) -> list[InputSlot]:
    accum = outdir / "accumulation"
    dirs_present = [
        d for d in ("top", "bottom", "all")
        if next(iter((accum / d / "randomization").glob("*_aggregated_results.csv")), None) is not None
    ]
    return [
        _slot("randomization", "Per-direction randomization CSVs", True,
              accum if dirs_present else None,
              note=f"directions found: {', '.join(dirs_present) or 'none'}"),
        _slot("global_csv", "Global aggregation CSV", True,
              _first_match(accum / "aggregation", "*_global.csv")),
    ]


def _accumulation_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    accum = report.html_path.parent.parent / "accumulation"
    output_file = "10.Accumulation_report.html"
    pre_lines = [
        "mkdir -p accum_root/top/randomization accum_root/bottom/randomization "
        "accum_root/all/randomization accum_root/aggregation",
    ]
    for direction in ("top", "bottom", "all"):
        pre_lines.append(
            f'for f in "{accum}/{direction}/randomization"/*_aggregated_results.csv; do '
            f'[ -f "$f" ] && cp "$f" accum_root/{direction}/randomization/; done'
        )
    pre_lines.append(f'for f in "{accum}/aggregation"/*_global.csv; do [ -f "$f" ] && cp "$f" accum_root/aggregation/; done')
    render = _render_call(
        "10.Accumulation_report.Rmd",
        [
            "accum_dir      = 'accum_root'",
            f"traitname      = '{report.traitname}'",
            "pval_threshold = 0.05",
            "output_dir     = '.'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="CT_ACCUMULATION — 10.Accumulation_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/CT_ACCUMULATION/local",
        pre_lines=pre_lines,
        render_block=render,
        output_file=output_file,
        publish_targets=[accum / "aggregation", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── CT_POSTPROC ───────────────────────────────────────────────────────────────


def _postproc_find_slots(outdir: Listing) -> list[InputSlot]:
    postproc = outdir / "postproc"
    return [
        _slot("discovery_file", "Discovery CSV (post-disambiguation)", True,
              _first_match(outdir, "ct_disambiguation/caas_convergence_master.csv", "**/*discovery*master*.csv")),
        _slot("filter_summary", "Filter summary TSV", True,
              _first_match(postproc, "summary_statistics/filter_summary.tsv")),
        _slot("gene_ensembl_file", "Gene Ensembl length file", True,
              _first_match(outdir, "postproc/postproc_inputs/*", "**/*gene_ensembl*", "**/*ensembl*length*")),
        _slot("gene_stats_file", "Gene stats file (discarded summary)", True,
              _first_match(postproc, "gene_filtering/discarded_summary.tsv")),
    ]


def _postproc_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    postproc = report.html_path.parent.parent / "postproc"
    output_file = "8.Characterization_report.html"
    render = _render_call(
        "8.Characterization_report.Rmd",
        [
            f"discovery_file = '{s['discovery_file']}'",
            f"filter_summary_file = '{s['filter_summary']}'",
            f"filter_dir = '{postproc}'",
            f"gene_ensembl_file = '{s['gene_ensembl_file']}'",
            f"gene_stats_file = '{s['gene_stats_file']}'",
            "processing_mode = 'filter'",
            "output_dir = '.'",
            "extreme_threshold = 3",
            "iqr_multiplier = 1.5",
            "alpha_threshold = 0.05",
            "gene_filter_mode = 'strict'",
            "seed = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="CT_POSTPROC — 8.Characterization_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/CT_POSTPROC/local",
        pre_lines=[
            'find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true',
        ],
        render_block=render,
        output_file=output_file,
        publish_targets=[postproc, report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── CT_META_CAAS ──────────────────────────────────────────────────────────────


def _meta_caas_find_slots(outdir: Listing) -> list[InputSlot]:
    # CAAS_META_CAAS_REPORT (workflows/ct_meta_caas.nf) takes the discovery.tab of CT
    # directly, not caas_convergence_master.csv (the input of the post-processing
    # characterization report).
    return [
        _slot("discovery_input", "Discovery input TSV", True,
              _first_match(outdir, "caastools/discovery.tab", "**/discovery.tab", "**/*discovery*.tab")),
        _slot("background_input", "Background genes file", True,
              _first_match(outdir, "**/*background_genes*")),
    ]


def _meta_caas_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "7.CAAS_pattern_annotation.html"
    render = _render_call(
        "7.CAAS_pattern_annotation.Rmd",
        [
            f"discovery_input = '{s['discovery_input']}'",
            f"background_input = '{s['background_input']}'",
            "output_dir = '.'",
            "caap_mode = FALSE",
            "seed = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="CT_META_CAAS — 7.CAAS_pattern_annotation.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/CT_META_CAAS/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        # Publishes under meta_caas/ only; a signification/ tree already present in
        # the outdir is left untouched.
        publish_targets=[report.html_path.parent.parent / "meta_caas", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── CAAS_SIGNIFICANCE (after SCORING; distinct from CT_META_CAAS) ─────────────
# Joins the meta_caas/global_meta_caas.tsv of CT_META_CAAS with the
# position_scores.tsv and gene_scores.tsv of SCORING, so it runs after scoring.


def _signif_significance_find_slots(outdir: Listing) -> list[InputSlot]:
    scoring = outdir / "scoring"
    slots = list(_meta_caas_find_slots(outdir))
    slots.append(
        _slot("global_meta_input", "Filtered discovery / CAAS TSV", True,
              _first_match(outdir, "postproc/gene_filtering/filtered_discovery.tsv",
                            "**/*filtered_discovery.tsv",
                            "meta_caas/meta_caas/global_meta_caas.tsv",
                            "signification/meta_caas/global_meta_caas.tsv",
                            "**/*global_meta_caas.tsv")
              or _first_match(outdir, "meta_caas/meta_caas/meta_caas.tsv",
                               "signification/meta_caas/meta_caas.tsv",
                               "**/*meta_caas.tsv"))
    )
    slots.append(_slot("position_scores_input", "Position scores TSV", True,
                        _first_match(scoring, "position_scores.tsv")))
    return slots


def _signif_significance_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "16.CAAS_significance_report.html"
    render = _render_call(
        "16.CAAS_significance_report.Rmd",
        [
            f"global_meta_input = '{s['global_meta_input']}'",
            f"position_scores_input = '{s['position_scores_input']}'",
            "output_dir = '.'",
            "seed = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="CAAS_SIGNIFICANCE — 16.CAAS_significance_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/CT_META_CAAS/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "meta_caas" / "significance", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── ASR_DIAGNOSTICS ────────────────────────────────────────────────────────────


def _asr_find_slots(outdir: Listing) -> list[InputSlot]:
    return [
        _slot("disambig_dir", "ct_disambiguation/ output directory", True,
              (outdir / "ct_disambiguation").path if (outdir / "ct_disambiguation").is_dir() else None),
    ]


def _asr_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "9.ASR_diagnostics.html"
    render = _render_call(
        "9.ASR_diagnostics.Rmd",
        [
            f"disambig_dir        = '{s['disambig_dir']}'",
            "output_dir          = '.'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="ASR_DIAGNOSTICS — 9.ASR_diagnostics.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/ASR_DIAGNOSTICS/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "asr_diagnostics", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── ENRICHMENT: FCS (CAAS report only) ────────────────────────────────────────
# RER_FCS_REPORT is not covered: its output filename comes from a caller-supplied
# report_label, so it cannot be recognized from the outdir (see UNSUPPORTED_NOTE).


def _fcs_scoring_find_slots(outdir: Listing) -> list[InputSlot]:
    fcs = outdir / "fcs"
    return [
        _slot("fcs_stats", "FCS stats TSV", True, _first_match(outdir, "**/fcs_stats.tsv")),
        _slot("universe", "Background/universe gene list", True,
              _first_match(outdir, "**/cleaned_background*.txt")),
        _slot("perms_file", "Permulation null file", False, _first_match(outdir, "**/*perm*.rds", "**/*perm*.tsv")),
        _slot("gene_lists", "Gene percentile lists directory", False,
              (outdir / "scoring" / "gene_lists").path if (outdir / "scoring" / "gene_lists").is_dir() else None),
    ]


def _fcs_scoring_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = f"12.FCS_scoring_{report.traitname}.html"
    render = _render_call(
        "12.FCS_general_report.Rmd",
        [
            f"stats_file    = '{s['fcs_stats']}'",
            f"universe_file = '{s['universe']}'",
            "gmt_dir       = paste0(REPO_DIR, '/subworkflows/ENRICHMENT/dat')",
            f"project_name  = 'Scoring_FCS_{report.traitname}'",
            "num_g         = 5",
            "fdr_thr       = 0.15",
            "pperm_thr     = 0.025",
            "top_n         = 25",
            f"traitname     = '{report.traitname}'",
            f"perms_file    = {_r_arg(s['perms_file'])}",
            f"gene_lists_dir = {_r_arg(s['gene_lists'])}",
            "seed           = '1998'",
        ],
        output_file,
    )
    # gmt_dir refers to REPO_DIR from inside the R params list, so REPO_DIR is set
    # before the render() call. It is single-quoted because the block sits inside a
    # double-quoted `Rscript -e "..."` bash string, which a double-quoted R literal
    # would close early.
    render_prefix = f"REPO_DIR <- '{repo_dir}'\n    "
    return _wrap_script(
        header_comment="ENRICHMENT — 12.FCS_general_report.Rmd (Scoring/CAAS)",
        repo_dir=repo_dir,
        local_dir="subworkflows/ENRICHMENT/local",
        pre_lines=[],
        render_block=render_prefix + render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "fcs", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── ENRICHMENT: POSENRICH ─────────────────────────────────────────────────────


def _posenrich_find_slots(outdir: Listing) -> list[InputSlot]:
    posenrich = outdir / "posenrich"
    return [
        # posenrich.nf emits this file (emit: results) as
        # posenrich_characterization.tsv, a name without "results".
        _slot("results", "Posenrich results TSV", True,
              _first_match(posenrich, "posenrich_characterization.tsv", "*results*.tsv")),
        _slot("leading_edge", "Posenrich leading-edge TSV", True, _first_match(posenrich, "*leading_edge*.tsv")),
        _slot("position_scores", "Position scores TSV", False, _first_match(outdir, "scoring/position_scores.tsv")),
        _slot("gene_scores", "Gene scores TSV", False, _first_match(outdir, "scoring/gene_scores.tsv")),
        _slot("vep_primateai", "VEP PrimateAI TSV", False, _first_match(outdir, "**/*primateai*.tsv")),
        _slot("vep_cosmic", "VEP COSMIC TSV", False, _first_match(outdir, "**/*cosmic*.tsv")),
        _slot("genomic_info", "Gene genomic coordinates TSV", False,
              _first_match(outdir, "postproc/postproc_inputs/*", "**/*genomic_info*.tsv")),
        _slot("fade_sites_top", "FADE per-site BF (top)", False,
              _first_match(outdir, "selection/fade/top/fade_site_bf_top.tsv")),
        _slot("fade_sites_bottom", "FADE per-site BF (bottom)", False,
              _first_match(outdir, "selection/fade/bottom/fade_site_bf_bottom.tsv")),
        _slot("position_lists", "Position percentile lists directory", False,
              (outdir / "scoring" / "position_lists").path if (outdir / "scoring" / "position_lists").is_dir() else None),
        _slot("position_biochem", "Position biochemistry TSV", False, _first_match(outdir, "scoring/position_biochem.tsv")),
        _slot("caas_pos_cycle_caas", "CAAS null cycle detail", False,
              _first_match(outdir, "caas/perm_pos_cycle_caas.tsv.gz", "scoring/caas_pos_cycle_caas.tsv.gz")),
        _slot("position_background", "CAAS tested positions background", False,
              _first_match(outdir, "caas/background.output", "discovery/background.output")),
        _slot("fcs_stats", "FCS stats TSV", False, _first_match(outdir, "**/fcs_stats.tsv")),
        _slot("cleaned_background", "Background/universe gene list", False,
              _first_match(outdir, "**/cleaned_background*.txt")),
    ]


def _posenrich_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = f"14.Position_enrichment_report_{report.traitname}.html"
    render = _render_call(
        "14.Position_enrichment_report.Rmd",
        [
            f"results_file = '{s['results']}'",
            f"leading_edge_file = '{s['leading_edge']}'",
            f"traitname = '{report.traitname}'",
            "padj_thr = 0.05",
            f"position_scores_file = {_r_arg(s['position_scores'])}",
            f"gene_scores_file     = {_r_arg(s['gene_scores'])}",
            f"vep_primateai_file   = {_r_arg(s['vep_primateai'])}",
            f"vep_cosmic_file      = {_r_arg(s['vep_cosmic'])}",
            f"genomic_info_file    = {_r_arg(s['genomic_info'])}",
            f"fade_sites_top_file    = {_r_arg(s['fade_sites_top'])}",
            f"fade_sites_bottom_file = {_r_arg(s['fade_sites_bottom'])}",
            f"position_biochem_file    = {_r_arg(s['position_biochem'])}",
            f"caas_pos_cycle_caas_file = {_r_arg(s['caas_pos_cycle_caas'])}",
            f"position_background_file = {_r_arg(s['position_background'])}",
            f"fcs_stats_file = {_r_arg(s['fcs_stats'])}",
            f"universe_file  = {_r_arg(s['cleaned_background'])}",
            f"position_lists_dir = {_r_arg(s['position_lists'])}",
            "seed = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="ENRICHMENT — 14.Position_enrichment_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/ENRICHMENT/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "posenrich", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── ENRICHMENT: AMI (SCORING_AMI_REPORT, DOMINO active modules) ───────────────


def _ami_find_slots(outdir: Listing) -> list[InputSlot]:
    scoring = outdir / "scoring"
    ami = outdir / "ami"
    return [
        _slot("gene_lists", "Gene percentile lists directory (scoring)", True,
              (scoring / "gene_lists").path if (scoring / "gene_lists").is_dir() else None),
        _slot("background", "Background/universe gene list", True,
              _first_match(outdir, "**/cleaned_background*.txt")),
        _slot("gene_scores", "Gene scores TSV", False, _first_match(scoring, "gene_scores.tsv")),
        _slot("domino_network_sif", "DOMINO network SIF", True, _first_match(ami, "**/network.sif", "**/*network*.sif")),
        _slot("domino_modules_dir", "DOMINO modules directory", True,
              _first_match(ami, "**/domino_modules")),
        _slot("domino_edge_scores", "DOMINO network edge scores TSV", True,
              _first_match(ami, "**/network_edge_scores.tsv")),
    ]


def _ami_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = f"13.AMI_analysis_{report.traitname}.html"
    bg_name = s["background"].name if s["background"] else "background.txt"
    pre_lines = [
        f'for f in "{s["gene_lists"]}"/slice_*.tsv; do',
        '    [ -f "$f" ] || continue',
        '    b=$(basename "$f" .tsv); name="${b#slice_}"',
        '    tail -n +2 "$f" | cut -f1 | grep -v "^[[:space:]]*$" > "${name}.txt" || true',
        "done",
    ]
    render = _render_call(
        "13.AMI_analysis.Rmd",
        [
            f"background_file     = '{s['background']}'",
            f"background_basename = '{bg_name}'",
            f"project_name        = '{report.traitname}'",
            "species             = 9606",
            "domino_network_score_thr = 700",
            f"gene_scores_file    = {_r_arg(s['gene_scores'])}",
            "string_db_dir       = NULL",
            f"domino_network_sif  = '{s['domino_network_sif']}'",
            f"domino_modules_dir  = '{s['domino_modules_dir']}'",
            f"domino_edge_scores_file = '{s['domino_edge_scores']}'",
            "fade_gene_lists_dir      = NULL",
            "fade_background_file     = NULL",
            "fade_domino_network_sif  = NULL",
            "fade_domino_modules_dir  = NULL",
            "fade_domino_edge_scores_file = NULL",
            "rer_gene_lists_dir       = NULL",
            "rer_background_file      = NULL",
            "rer_domino_network_sif   = NULL",
            "rer_domino_modules_dir   = NULL",
            "rer_domino_edge_scores_file = NULL",
            "seed                     = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="ENRICHMENT — 13.AMI_analysis.Rmd (FADE/RER cross-module sections left NULL — "
                        "regenerate their own AMI reports if you need those sections back)",
        repo_dir=repo_dir,
        local_dir="subworkflows/ENRICHMENT/local",
        pre_lines=pre_lines,
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "ami", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── ENRICHMENT: COMPARE (SCORING_COMPARE_REPORT) ──────────────────────────────


def _compare_find_slots(outdir: Listing) -> list[InputSlot]:
    # fcs.nf publishes into fcs/fcs_results/ with pattern 'fcs_results/**'; the
    # emitted paths already start with "fcs_results/", so the files land one level
    # deeper: fcs/fcs_results/fcs_results/*.tsv.
    fcs = outdir / "fcs" / "fcs_results" / "fcs_results"
    return [
        _slot("caas_fcs", "CAAS FCS all-results TSV", False, _first_match(fcs, "fcs_all_results.tsv")),
        _slot("rer_fcs", "RER FCS all-results TSV", False,
              _first_match(outdir, "rerconverge/**/fcs_results/fcs_all_results.tsv")),
        _slot("caas_le", "CAAS leading-edge TSV", False, _first_match(fcs, "fcs_leading_edge.tsv")),
        _slot("rer_le", "RER leading-edge TSV", False,
              _first_match(outdir, "rerconverge/**/fcs_results/fcs_leading_edge.tsv")),
        _slot("caas_le_comp", "CAAS leading-edge composition TSV", False,
              _first_match(fcs, "fcs_leading_edge_composition.tsv")),
        _slot("rer_le_comp", "RER leading-edge composition TSV", False,
              _first_match(outdir, "rerconverge/**/fcs_results/fcs_leading_edge_composition.tsv")),
        # scoring_enrichment.nf publishes into ami/ami_networks/ with pattern
        # 'ami_networks/**': the same double nesting as fcs above.
        _slot("ami_module_desc", "AMI module descriptions TSV", False,
              _first_match(outdir, "ami/ami_networks/ami_networks/ami_module_descriptions_all_tools.tsv")),
        _slot("ami_term_membership", "AMI term membership TSV", False,
              _first_match(outdir, "ami/ami_networks/ami_networks/ami_term_threshold_membership_all_tools.tsv")),
        _slot("posenrich_dotplot", "Posenrich overall dotplot TSV", False,
              _first_match(outdir, "posenrich/posenrich_overall_dotplot.tsv")),
        _slot("posenrich_le", "Posenrich leading-edge TSV", False, _first_match(outdir, "**/posenrich_leading_edge.tsv")),
        _slot("posenrich_le_summary", "Posenrich leading-edge summary TSV", False,
              _first_match(outdir, "posenrich/posenrich_leading_edge_summary.tsv")),
        _slot("fcs_stats", "FCS stats TSV", False, _first_match(outdir, "**/fcs_stats.tsv")),
        _slot("position_scores", "Position scores TSV", False, _first_match(outdir, "scoring/position_scores.tsv")),
        _slot("vep_primateai", "VEP PrimateAI TSV", False, _first_match(outdir, "**/*primateai*.tsv")),
        _slot("vep_cosmic", "VEP COSMIC TSV", False, _first_match(outdir, "**/*cosmic*.tsv")),
        _slot("fade_sites_top", "FADE per-site BF (top)", False,
              _first_match(outdir, "selection/fade/top/fade_site_bf_top.tsv")),
        _slot("fade_sites_bottom", "FADE per-site BF (bottom)", False,
              _first_match(outdir, "selection/fade/bottom/fade_site_bf_bottom.tsv")),
        _slot("gene_lists", "Gene percentile lists directory", False,
              (outdir / "scoring" / "gene_lists").path if (outdir / "scoring" / "gene_lists").is_dir() else None),
        _slot("position_lists", "Position percentile lists directory", False,
              (outdir / "scoring" / "position_lists").path if (outdir / "scoring" / "position_lists").is_dir() else None),
    ]


def _compare_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = f"15.Comparison_report_{report.traitname}.html"
    pre_lines = ["mkdir -p cmp_fcs cmp_fcs_le"]
    if s["caas_fcs"]:
        pre_lines.append(f'cp "{s["caas_fcs"]}" cmp_fcs/caas_fcs_all.tsv')
    if s["rer_fcs"]:
        pre_lines.append(f'cp "{s["rer_fcs"]}" cmp_fcs/rer_fcs_all.tsv')
    if s["caas_le"]:
        pre_lines.append(f'cp "{s["caas_le"]}" cmp_fcs_le/caas_leading_edge.tsv')
    if s["rer_le"]:
        pre_lines.append(f'cp "{s["rer_le"]}" cmp_fcs_le/rer_leading_edge.tsv')
    if s["caas_le_comp"]:
        pre_lines.append(f'cp "{s["caas_le_comp"]}" cmp_fcs_le/caas_leading_edge_composition.tsv')
    if s["rer_le_comp"]:
        pre_lines.append(f'cp "{s["rer_le_comp"]}" cmp_fcs_le/rer_leading_edge_composition.tsv')

    render = _render_call(
        "15.Comparison_report.Rmd",
        [
            "fcs_dir    = 'cmp_fcs'",
            "fcs_le_dir = 'cmp_fcs_le'",
            "fdr_thr    = 0.15",
            "pperm_thr  = 0.025",
            "domino_module_thr = 0.05",
            "top_n      = 100",
            f"traitname  = '{report.traitname}'",
            f"ami_module_desc_file      = {_r_arg(s['ami_module_desc'])}",
            f"ami_term_membership_file  = {_r_arg(s['ami_term_membership'])}",
            f"posenrich_dotplot_file     = {_r_arg(s['posenrich_dotplot'])}",
            f"posenrich_leading_edge_file = {_r_arg(s['posenrich_le'])}",
            f"posenrich_leading_edge_summary_file = {_r_arg(s['posenrich_le_summary'])}",
            f"fcs_stats_file           = {_r_arg(s['fcs_stats'])}",
            f"position_scores_file     = {_r_arg(s['position_scores'])}",
            f"vep_primateai_file       = {_r_arg(s['vep_primateai'])}",
            f"vep_cosmic_file          = {_r_arg(s['vep_cosmic'])}",
            f"fade_sites_top_file      = {_r_arg(s['fade_sites_top'])}",
            f"fade_sites_bottom_file   = {_r_arg(s['fade_sites_bottom'])}",
            f"gene_lists_dir     = {_r_arg(s['gene_lists'])}",
            f"position_lists_dir = {_r_arg(s['position_lists'])}",
            "seed               = '1998'",
        ],
        output_file,
    )
    return _wrap_script(
        header_comment="ENRICHMENT — 15.Comparison_report.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/ENRICHMENT/local",
        pre_lines=pre_lines,
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent.parent / "compare", report.html_path.parent],
        use_singularity=use_singularity,
    )


# ── TRAIT_ANALYSIS (5 reports) ────────────────────────────────────────────────
# trait_file and tree_file are inputs of the pipeline, not outputs, so they are
# found only when a copy (original or pruned) sits under outdir; otherwise the user
# selects them.


def _ta_common_slots(outdir: Listing, extra_dir_glob: Optional[str] = None) -> list[InputSlot]:
    # pruned_trait_file.tsv exists only when pruning ran (--prune_data, ta_data_prune.nf).
    # original_trait_file.tsv is the copy of the raw trait file made by
    # ta_dataset_exploration.nf: the fallback when pruning was skipped.
    slots = [
        _slot("trait_file", "Original trait file (--trait_file)", True,
              _first_match(outdir, "**/pruned_trait_file.tsv", "**/original_trait_file.tsv"),
              note="original pipeline input — browse to it if not auto-found"),
        _slot("tree_file", "Original tree file (--tree_file)", True,
              _first_match(outdir, "**/pruned_tree_file.nwk", "**/original_tree_file.nwk"),
              note="original pipeline input — browse to it if not auto-found"),
    ]
    return slots


def _ta_data_prune_find_slots(outdir: Listing) -> list[InputSlot]:
    return _ta_common_slots(outdir)


def _ta_data_prune_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "2.Phenotype_exploration_pruned.html"
    render = _render_call(
        "0.Data_pruning.Rmd",
        [
            f"trait_file = '{s['trait_file']}'",
            f"tree_file = '{s['tree_file']}'",
            "output_dir = '.'",
            f"traitname = '{report.traitname}'",
            "seed = ''", "clade_name = ''", "taxon_of_interest = ''",
            "n_trait = ''", "c_trait = ''", "tax_id = ''",
            "secondary_trait = ''", "branch_trait = ''",
            "prune_list = ''", "prune_list_secondary = ''",
            "trait_type = ''",
            "pss_top_pct = '0.05'", "perm_strategy = 'best_model'",
            "max_contrasts = '0'",
        ],
        output_file,
        extra_arg="envir = new.env()",
    )
    return _wrap_script(
        header_comment="TRAIT_ANALYSIS — 0.Data_pruning.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/TRAIT_ANALYSIS/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent],
        use_singularity=use_singularity,
    )


def _ta_dataset_exploration_find_slots(outdir: Listing) -> list[InputSlot]:
    slots = _ta_common_slots(outdir)
    slots.append(_slot("prune_results_dir", "Pruned data_exploration/ directory", False,
                        _first_match(outdir, "data_exploration")))
    return slots


def _ta_dataset_exploration_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "1.Dataset_exploration.html"
    pre_lines = ["mkdir -p data_exploration"]
    if s["prune_results_dir"]:
        pre_lines.append(f'cp -R "{s["prune_results_dir"]}"/* data_exploration/ 2>/dev/null || true')
    render = _render_call(
        "1.Dataset_exploration.Rmd",
        [
            f"trait_file = '{s['trait_file']}'",
            f"tree_file = '{s['tree_file']}'",
            "output_dir = 'data_exploration'",
            f"traitname = '{report.traitname}'",
            "seed = ''", "clade_name = ''", "taxon_of_interest = ''",
            "n_trait = ''", "c_trait = ''", "tax_id = ''",
            "secondary_trait = ''", "branch_trait = ''",
            "trait_type = ''",
            "pss_top_pct = '0.05'", "perm_strategy = 'best_model'",
        ],
        output_file,
        extra_arg="envir = new.env()",
    )
    return _wrap_script(
        header_comment="TRAIT_ANALYSIS — 1.Dataset_exploration.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/TRAIT_ANALYSIS/local",
        pre_lines=pre_lines,
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent],
        use_singularity=use_singularity,
    )


def _ta_phenotype_exploration_find_slots(outdir: Listing) -> list[InputSlot]:
    slots = _ta_common_slots(outdir)
    slots.append(_slot("results_dir", "data_exploration/ directory", True,
                        _first_match(outdir, "data_exploration")))
    return slots


def _ta_phenotype_exploration_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "2.Phenotype_exploration_complete.html"
    pre_lines = ["mkdir -p data_exploration", f'cp -R "{s["results_dir"]}"/* data_exploration/']
    render = _render_call(
        "2.Phenotype_exploration.Rmd",
        [
            f"trait_file = '{s['trait_file']}'",
            f"tree_file = '{s['tree_file']}'",
            "output_dir = 'data_exploration'",
            f"traitname = '{report.traitname}'",
            "seed = ''", "clade_name = ''", "taxon_of_interest = ''",
            "n_trait = ''", "c_trait = ''", "tax_id = ''",
            "secondary_trait = ''", "branch_trait = ''",
            "trait_type = ''",
            "pss_top_pct = '0.05'", "perm_strategy = 'best_model'",
            "max_contrasts = '0'",
        ],
        output_file,
        extra_arg="envir = new.env()",
    )
    return _wrap_script(
        header_comment="TRAIT_ANALYSIS — 2.Phenotype_exploration.Rmd (regeneration always applies the "
                        "singularity branch's results_dir merge step, fixing a bare-mode asymmetry in "
                        "ta_phenotype_exploration.nf)",
        repo_dir=repo_dir,
        local_dir="subworkflows/TRAIT_ANALYSIS/local",
        pre_lines=pre_lines,
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent],
        use_singularity=use_singularity,
    )


def _ta_ci_find_slots(outdir: Listing) -> list[InputSlot]:
    slots = _ta_common_slots(outdir)
    slots.append(_slot("results_dir", "data_exploration/ directory", True,
                        _first_match(outdir, "data_exploration")))
    return slots


def _ta_ci_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "3.CI-composition.html"
    render = _render_call(
        "3.CI-composition.Rmd",
        [
            f"trait_file = '{s['trait_file']}'",
            f"tree_file = '{s['tree_file']}'",
            f"output_dir = '{s['results_dir']}'",
            f"traitname = '{report.traitname}'",
            "seed = ''", "clade_name = ''", "taxon_of_interest = ''",
            "n_trait = ''", "c_trait = ''", "tax_id = ''",
            "secondary_trait = ''", "branch_trait = ''",
            "trait_type = ''",
            "pss_top_pct = '0.05'", "perm_strategy = 'best_model'",
        ],
        output_file,
        extra_arg="envir = new.env()",
    )
    return _wrap_script(
        header_comment="TRAIT_ANALYSIS — 3.CI-composition.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/TRAIT_ANALYSIS/local",
        pre_lines=[],
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent],
        use_singularity=use_singularity,
    )


def _ta_contrast_find_slots(outdir: Listing) -> list[InputSlot]:
    slots = _ta_common_slots(outdir)
    slots.append(_slot("results_dir", "data_exploration/ directory", True,
                        _first_match(outdir, "data_exploration")))
    return slots


def _ta_contrast_build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    s = {slot.key: slot.path for slot in report.slots}
    output_file = "4.Independent_contrasts.html"
    render = _render_call(
        "4.Independent_contrasts.Rmd",
        [
            f"trait_file = '{s['trait_file']}'",
            f"tree_file = '{s['tree_file']}'",
            f"output_dir = '{s['results_dir']}'",
            f"traitname = '{report.traitname}'",
            "seed = ''", "clade_name = ''", "taxon_of_interest = ''",
            "n_trait = ''", "c_trait = ''", "tax_id = ''",
            "secondary_trait = ''", "branch_trait = ''",
            "trait_type = ''",
            "pss_top_pct = '0.05'", "perm_strategy = 'best_model'",
            "max_contrasts = '0'",
        ],
        output_file,
        extra_arg="envir = new.env()",
    )
    post_lines = [
        f'mkdir -p "{s["results_dir"]}/2.CT/3.Tree"',
        f'cp "{s["tree_file"]}" "{s["results_dir"]}/2.CT/3.Tree/pruned_tree_file.nwk"',
    ]
    # selection_algorithm.R (sourced by 4.Independent_contrasts.Rmd) locates
    # lean_contrast_selector.R at a list of candidate paths, most of them relative to
    # the full repository tree. The scratch directory holds only the staged
    # subworkflows/TRAIT_ANALYSIS/local/ (see _wrap_script), so the script is copied
    # to the first candidate, ./src/, instead of reproducing the tree.
    pre_lines = [
        f'cp "{repo_dir}/subworkflows/CT/local/scripts/lean_contrast_selector.R" src/',
    ]
    return _wrap_script(
        header_comment="TRAIT_ANALYSIS — 4.Independent_contrasts.Rmd",
        repo_dir=repo_dir,
        local_dir="subworkflows/TRAIT_ANALYSIS/local",
        pre_lines=pre_lines,
        render_block=render,
        output_file=output_file,
        publish_targets=[report.html_path.parent],
        use_singularity=use_singularity,
        post_lines=post_lines,
    )


# ── Registry ──────────────────────────────────────────────────────────────────

REPORTS: list[ReportSpec] = [
    ReportSpec(
        id="scoring",
        display_name="SCORING — Scoring report",
        html_regex=re.compile(r"^11\.Scoring_report_(?P<trait>.+)\.html$"),
        find_slots=_scoring_find_slots,
        build_script=_scoring_build_script,
        guess_traitname=lambda name: _guess_trait_from_suffix(name, "11.Scoring_report_"),
    ),
    ReportSpec(
        id="rer",
        display_name="RERCONVERGE — RER report",
        html_regex=re.compile(r"^5\.RERconverge_report\.html$"),
        find_slots=_rer_find_slots,
        build_script=_rer_build_script,
    ),
    _make_fade_spec("top"),
    _make_fade_spec("bottom"),
    ReportSpec(
        id="accumulation",
        display_name="CT_ACCUMULATION — Accumulation report",
        html_regex=re.compile(r"^10\.Accumulation_report\.html$"),
        find_slots=_accumulation_find_slots,
        build_script=_accumulation_build_script,
    ),
    ReportSpec(
        id="postproc",
        display_name="CT_POSTPROC — Characterization report",
        html_regex=re.compile(r"^8\.Characterization_report\.html$"),
        find_slots=_postproc_find_slots,
        build_script=_postproc_build_script,
    ),
    ReportSpec(
        id="meta_caas",
        display_name="CT_META_CAAS — CAAS Pattern Annotation",
        # Also matches 7.CT_signification.html and 7.CT_Pattern_Annotation.html, the
        # other names under which this report appears in existing outdirs.
        html_regex=re.compile(r"^7\.(?:CT_signification|CAAS_pattern_annotation|CT_Pattern_Annotation)\.html$"),
        find_slots=_meta_caas_find_slots,
        build_script=_meta_caas_build_script,
    ),
    ReportSpec(
        id="signif_significance",
        display_name="CT_META_CAAS — CAAS significance report (post-scoring)",
        html_regex=re.compile(r"^16\.CAAS_significance_report\.html$"),
        find_slots=_signif_significance_find_slots,
        build_script=_signif_significance_build_script,
    ),
    ReportSpec(
        id="asr_diagnostics",
        display_name="ASR_DIAGNOSTICS — Diagnostics report",
        html_regex=re.compile(r"^9\.ASR_diagnostics\.html$"),
        find_slots=_asr_find_slots,
        build_script=_asr_build_script,
    ),
    ReportSpec(
        id="fcs_scoring",
        display_name="ENRICHMENT — FCS report (Scoring/CAAS)",
        html_regex=re.compile(r"^12\.FCS_scoring_(?P<trait>.+)\.html$"),
        find_slots=_fcs_scoring_find_slots,
        build_script=_fcs_scoring_build_script,
        guess_traitname=lambda name: _guess_trait_from_suffix(name, "12.FCS_scoring_"),
    ),
    ReportSpec(
        id="posenrich",
        display_name="ENRICHMENT — Position enrichment report",
        html_regex=re.compile(r"^14\.Position_enrichment_report_(?P<trait>.+)\.html$"),
        find_slots=_posenrich_find_slots,
        build_script=_posenrich_build_script,
        guess_traitname=lambda name: _guess_trait_from_suffix(name, "14.Position_enrichment_report_"),
    ),
    ReportSpec(
        id="ami",
        display_name="ENRICHMENT — AMI (DOMINO modules) report",
        html_regex=re.compile(r"^13\.AMI_analysis_(?P<trait>.+)\.html$"),
        find_slots=_ami_find_slots,
        build_script=_ami_build_script,
        guess_traitname=lambda name: _guess_trait_from_suffix(name, "13.AMI_analysis_"),
    ),
    ReportSpec(
        id="compare",
        display_name="ENRICHMENT — Comparison report",
        html_regex=re.compile(r"^15\.Comparison_report_(?P<trait>.+)\.html$"),
        find_slots=_compare_find_slots,
        build_script=_compare_build_script,
        guess_traitname=lambda name: _guess_trait_from_suffix(name, "15.Comparison_report_"),
    ),
    ReportSpec(
        id="ta_data_prune",
        display_name="TRAIT_ANALYSIS — Data pruning report",
        html_regex=re.compile(r"^2\.Phenotype_exploration_pruned\.html$"),
        find_slots=_ta_data_prune_find_slots,
        build_script=_ta_data_prune_build_script,
    ),
    ReportSpec(
        id="ta_dataset_exploration",
        display_name="TRAIT_ANALYSIS — Dataset exploration report",
        html_regex=re.compile(r"^1\.Dataset_exploration\.html$"),
        find_slots=_ta_dataset_exploration_find_slots,
        build_script=_ta_dataset_exploration_build_script,
    ),
    ReportSpec(
        id="ta_phenotype_exploration",
        display_name="TRAIT_ANALYSIS — Phenotype exploration report",
        html_regex=re.compile(r"^2\.Phenotype_exploration_complete\.html$"),
        find_slots=_ta_phenotype_exploration_find_slots,
        build_script=_ta_phenotype_exploration_build_script,
    ),
    ReportSpec(
        id="ta_ci",
        display_name="TRAIT_ANALYSIS — CI-composition report",
        html_regex=re.compile(r"^3\.CI-composition\.html$"),
        find_slots=_ta_ci_find_slots,
        build_script=_ta_ci_build_script,
    ),
    ReportSpec(
        id="ta_contrast",
        display_name="TRAIT_ANALYSIS — Independent contrasts report",
        html_regex=re.compile(r"^4\.Independent_contrasts\.html$"),
        find_slots=_ta_contrast_find_slots,
        build_script=_ta_contrast_build_script,
    ),
]

# Reports that cannot be detected automatically, shown as a note in the regenerate
# dialog instead of being omitted: RER_FCS_REPORT (subworkflows/ENRICHMENT/fcs.nf)
# takes a caller-supplied report_label, so its HTML filename has no fixed pattern.
UNSUPPORTED_NOTE = (
    "Not auto-detected: RER's own FCS report (RER_FCS_REPORT) uses a "
    "caller-chosen output filename with no fixed pattern to match — "
    "regenerate it manually if needed."
)


def detect_reports(
    outdir: str,
    repo_dir: Path,
    traitname_override: str = "",
    remote_entries: Optional[dict[str, bool]] = None,
) -> list[DetectedReport]:
    """Return one DetectedReport per report HTML found in {outdir}/html_reports/.

    Every report process publishes its HTML there as well as in its own module
    subdirectory. The input slots of each report are resolved on a best-effort basis.
    The local filesystem is scanned by default; pass remote_entries (the {relpath:
    is_dir} map of everything under outdir, from one gui.remote.list_all_remote()
    round trip) to scan a remote outdir instead. The function makes no network call.
    """
    root = Listing(str(outdir), remote_entries)
    html_dir = root / "html_reports"
    if not html_dir.is_dir():
        return []

    found: list[DetectedReport] = []
    for html_file in sorted(html_dir.glob("*.html")):
        for spec in REPORTS:
            m = spec.html_regex.match(html_file.name)
            if not m:
                continue
            traitname = traitname_override.strip() or spec.guess_traitname(html_file.name)
            slots = spec.find_slots(root)
            found.append(
                DetectedReport(
                    id=spec.id,
                    display_name=spec.display_name,
                    html_path=html_file,
                    rmd_filename="",
                    local_dir="",
                    traitname=traitname,
                    slots=slots,
                    spec=spec,
                    script_filename=f"regenerate_{spec.id}.sh",
                )
            )
            break
    return found


def build_script(report: DetectedReport, repo_dir: Path, use_singularity: bool) -> str:
    """Build the re-rendering script of a detected report (dispatches to its spec)."""
    return report.spec.build_script(report, repo_dir, use_singularity)
