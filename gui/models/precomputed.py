#!/usr/bin/env python3
# precomputed.py — Precomputed-input reuse, consolidated from every module tab.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Precomputed: reuse of outputs from an earlier run, as one base_path plus one checkbox
per producing stage.

A single global path per input cannot serve a batch with several phenotypes, because
the TRAIT changes per row. Every file the pipeline writes sits in a stable layout under
an outdir (RESULTS_BASE in run_single.sh.j2), so the per-phenotype path is always
base_path/<TRAIT>/<subpath>. The generated shell script builds it from $TRAIT
(PRECOMP_* variables of sbatch_array.sh.j2 and run_single.sh.j2); nothing is typed per path.

Checking a box supplies the precomputed input and turns off the module that would
recompute it (see the PRECOMP_* wiring in gui/generation/templates/run_single.sh.j2):
feeding a precomputed result while also recomputing it is never intended.

Path templates, relative to base_path/<TRAIT>, matched to the publishDir of each process:
  CT/CAAS       : caastools/{discovery,resample}.tab, caastools/background_genes.output,
                  caastools/background.output, meta_caas/meta_caas/meta_caas.tsv
                  (signification/meta_caas/ is also searched, for outdirs that use that name),
                  caas_permulation/ (dir: null_positions.parquet, blocks/, cycle_roster.txt, null_metadata.json, counts.json,
                  cycle_summary.tsv, cycle_hist.tsv)
  Disambiguation: ct_disambiguation/caas_convergence_master.csv, ct_disambiguation/ (dir)
  Post-processing: postproc/gene_filtering/filtered_discovery.tsv,
                  postproc/cleaned_backgrounds/cleaned_background_main.txt
  RER           : rerconverge/rer_results/<TRAIT>.continuous.{output,perms.rds},
                  rerconverge/rer_results/rerconverge_summary_<TRAIT>.tsv (resolved by glob
                  at generation time, because the filename carries a variable suffix)
  FADE          : selection/fade/{top,bottom} (dir),
                  selection/fade/{top,bottom}/fade_summary_{top,bottom}.tsv,
                  selection/fade/{top,bottom}/fade_site_bf_{top,bottom}.tsv
  VEP           : vep/primateai_mapped.tsv, vep/cosmic_scores.tsv, vep/ensembl_vep_mapped.tsv

Some of these files are published conditionally (the post-processing pair only when
gene_filter_mode != 'none'; the optional emits for RER perms, FADE summary and site
tables and VEP scores). Checking the box
wires the path in any case, but the file exists only if the source run had the matching
settings on.

Imported by: gui/models/project.py, gui/widgets/tabs/precomputed_tab.py,
gui/generation/validate.py (derive_paths)
"""

# ── Standard library ──────────────────────────────────────────────────────────
import glob
import os
from dataclasses import dataclass


@dataclass(kw_only=True)
class PrecomputedConfig:
    """Base path and per-stage reuse toggles. The module docstring gives the path layout."""

    base_path: str = ""  # per-phenotype dir = base_path/<TRAIT> (no toy/postproc-mode tag)

    # --- CT / CAAS (one general toggle and the two step toggles of the CAAS tab:
    # discovery and resample) ---
    use_ct: bool = False  # turns off CAAS; also wires background_input, meta_caas_from,
    # posenrich_background_file, which come from CT's concatenation and
    # permulation outputs rather than from either single step.
    use_discovery: bool = False  # wires discovery_from
    use_resample: bool = False  # wires resample_from

    # --- Disambiguation (ASR and convergence computation) ---
    use_disambiguation: bool = False  # turns off Disambiguation; wires disambiguation_input/_dir

    # --- Post-processing (no enable toggle: it runs with Disambiguation unless
    # this box supplies its outputs; see ct_postproc_enabled in gui/generation/context.py) ---
    use_postproc: bool = False  # wires the filtered-discovery/cleaned-background pair
    # that Accumulation (POSENRICH), VEP and Scoring take as input; it does not turn off
    # Disambiguation (use_disambiguation does)

    # --- RERconverge ---
    use_rer: bool = False  # turns off RER; wires rer_continuous_file/rer_perms_file plus
    # the per-phenotype scoring_rer_input/scoring_rer_perms_input Scoring falls back to

    # --- FADE ---
    use_fade: bool = False  # turns off FADE; wires fade_tables_dir_top/bottom plus the
    # per-phenotype scoring_fade_summary_top/bottom and scoring_fade_site_top/bottom

    # --- VEP ---
    use_vep: bool = False  # turns off VEP; wires scoring_vep_primateai/scoring_vep_cosmic/scoring_vep_ensembl


def derive_paths(config: "PrecomputedConfig", trait: str) -> list[tuple[str, str, str]]:
    """List (label, path, kind) for every path implied by the checked boxes, for one TRAIT.

    kind is "file", "dir", or "path" (either kind). Python mirror of the PRECOMP_OUTDIR construction of
    run_single.sh.j2, used by gui/generation/validate.py for existence checks. The shell
    script builds the same paths at run time, so both must follow the layout in this
    module's docstring.
    """
    if not config.base_path or not trait:
        return []
    outdir = os.path.join(config.base_path, trait)
    if not os.path.isdir(outdir) and os.path.isdir(f"{outdir}_complete"):
        outdir = f"{outdir}_complete"
    entries: list[tuple[str, str, str]] = []

    if config.use_discovery:
        entries.append(("discovery_from", os.path.join(outdir, "caastools", "discovery.tab"), "file"))
        entries.append(("meta_caas_from", os.path.join(outdir, "meta_caas", "meta_caas", "meta_caas.tsv"), "file"))
        entries.append(("background_input", os.path.join(outdir, "caastools", "background_genes.output"), "file"))
        entries.append(("posenrich_background_file", os.path.join(outdir, "caastools", "background.output"), "file"))
        # The position-level permulation null travels as one directory: the positions the null scores, its blocks of
        # scored positions, the roster and the summaries behind the position p-values and the null panels of the
        # report. Mirrors the PRECOMP_USE_DISCOVERY block of run_single.sh.j2.
        entries.append(("caas_null_dir", os.path.join(outdir, "caas_permulation"), "dir"))
        entries.append(("caas_perm_manifest_file", os.path.join(outdir, "caastools", "permulation_manifest.tsv"), "file"))
        entries.append(("caas_perm_harvest_file", os.path.join(outdir, "caastools", "permulation_harvest.tsv"), "file"))
    if config.use_resample:
        entries.append(("resample_from", os.path.join(outdir, "caastools", "resample.tab"), "file"))

    if config.use_disambiguation:
        entries.append(
            ("disambiguation_input", os.path.join(outdir, "ct_disambiguation", "caas_convergence_master.csv"), "file")
        )
        entries.append(("disambiguation_dir", os.path.join(outdir, "ct_disambiguation"), "dir"))

    if config.use_postproc:
        gene_list = os.path.join(outdir, "postproc", "gene_filtering", "filtered_discovery.tsv")
        background = os.path.join(outdir, "postproc", "cleaned_backgrounds", "cleaned_background_main.txt")
        entries.append(("accumulation_caas_input / scoring_postproc_input", gene_list, "file"))
        entries.append(("accumulation_background_input / scoring_background_input", background, "file"))

    if config.use_rer:
        rer_dir = os.path.join(outdir, "rerconverge", "rer_results")
        # Same rule as run_single.sh.j2: binary output only when no continuous one exists.
        binary_only = os.path.isfile(os.path.join(rer_dir, f"{trait}.binary.output")) and not os.path.isfile(
            os.path.join(rer_dir, f"{trait}.continuous.output")
        )
        kind = "binary" if binary_only else "continuous"
        entries.append((f"rer_{kind}_file", os.path.join(rer_dir, f"{trait}.{kind}.output"), "file"))
        entries.append(("rer_perms_file / scoring_rer_perms_input", os.path.join(rer_dir, f"{trait}.{kind}.perms.rds"), "path"))
        matches = sorted(glob.glob(os.path.join(rer_dir, "rerconverge_summary_*.tsv")))
        if matches:
            entries.append(("scoring_rer_input", matches[0], "file"))

    if config.use_fade:
        for direction in ("top", "bottom"):
            fade_dir = os.path.join(outdir, "selection", "fade", direction)
            entries.append((f"fade_tables_dir_{direction}", fade_dir, "dir"))
            entries.append((f"scoring_fade_summary_{direction}", os.path.join(fade_dir, f"fade_summary_{direction}.tsv"), "file"))
            entries.append((f"scoring_fade_site_{direction}", os.path.join(fade_dir, f"fade_site_bf_{direction}.tsv"), "file"))

    if config.use_vep:
        vep_dir = os.path.join(outdir, "vep")
        entries.append(("scoring_vep_primateai", os.path.join(vep_dir, "primateai_mapped.tsv"), "file"))
        entries.append(("scoring_vep_cosmic", os.path.join(vep_dir, "cosmic_scores.tsv"), "file"))
        entries.append(("scoring_vep_ensembl", os.path.join(vep_dir, "ensembl_vep_mapped.tsv"), "file"))

    return entries
