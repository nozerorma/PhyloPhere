#!/usr/bin/env python3
# validate.py — Pre-render validation: required fields per enabled module, CLASS rules.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure function, no PySide6 import. Returns a list of human-readable error strings;
an empty list means the project is renderable. Used both by tests/gui and by the
GUI's "Generate Scripts..." action (which should show these in a dialog rather than
letting render.py fail with a StrictUndefined KeyError deep inside a template).

Scope is deliberately the "essential fields" surface only (see implementation plan
§5) — this checks the fields actually exposed as GUI widgets, not the full
conf/*.config parameter space.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig


def validate(project: ProjectConfig) -> list[str]:
    errors: list[str] = []

    def require(value: str, message: str) -> None:
        if not value.strip():
            errors.append(message)

    # --- General ---
    require(project.general.repo_dir, "General: PhyloPhere repo directory is required.")
    require(
        project.general.nextflow_plugins_dir,
        "General: Nextflow plugins directory is required (symlinked into every run's NXF_HOME).",
    )

    # --- Runtime ---
    rt = project.runtime
    require(rt.alignment_dir, "Runtime: alignment directory is required.")
    require(rt.tree_file, "Runtime: species tree file is required.")
    require(rt.work_dir, "Runtime: work directory is required.")
    require(rt.results_dir, "Runtime: results directory is required.")

    if not rt.phenotype_rows:
        errors.append("Runtime: the phenotype catalogue must have at least one row.")

    has_class1 = any(row.trait_class == 1 for row in rt.phenotype_rows)
    has_class2 = any(row.trait_class == 2 for row in rt.phenotype_rows)
    if has_class1:
        require(rt.trait_file, "Runtime: trait file is required (a CLASS 1 phenotype row exists).")
        require(rt.prune_dir, "Runtime: prune directory is required (a CLASS 1 phenotype row exists).")
    if has_class2:
        require(
            rt.simple_trait_file,
            "Runtime: simple trait file is required (a CLASS 2 phenotype row exists).",
        )

    for i, row in enumerate(rt.phenotype_rows, start=1):
        require(row.trait, f"Phenotype row {i}: trait name is required.")
        if row.trait_class not in (1, 2):
            errors.append(f"Phenotype row {i}: CLASS must be 1 or 2, got {row.trait_class!r}.")
        if row.trait_class == 2:
            require(row.discrete_method, f"Phenotype row {i}: discrete method is required for CLASS 2.")

    # --- CAAS ---
    caas = project.modules.caas
    if caas.enabled:
        require(caas.caas_config_path, "CAAS: config file is required when CAAS is enabled.")
    else:
        if not any([caas.discovery_from, caas.resample_from, caas.bootstrap_from]):
            errors.append(
                "CAAS is disabled but no discovery_from/resample_from/bootstrap_from fallback "
                "was supplied — downstream modules (Disambiguation, Accumulation) have no input."
            )

    # --- Disambiguation (+ Post-processing) ---
    disambig = project.modules.disambiguation
    if disambig.enabled:
        if disambig.ct_disambig_asr_mode == "precomputed":
            require(
                disambig.ct_disambig_asr_cache_dir,
                "Disambiguation: ASR cache directory is required when asr_mode=precomputed.",
            )
        if not caas.enabled and not any(
            [caas.discovery_from, caas.resample_from, caas.bootstrap_from]
        ):
            errors.append(
                "Disambiguation is enabled but CAAS is disabled with no fallback input supplied."
            )
    else:
        if not any(
            [
                disambig.signification_from,
                disambig.disambiguation_input,
                disambig.disambiguation_dir,
                disambig.background_input,
            ]
        ):
            errors.append(
                "Disambiguation is disabled but no fallback input was supplied — Accumulation "
                "and Scoring may have no input."
            )

    # --- Accumulation ---
    accum = project.modules.accumulation
    if accum.enabled:
        require(
            accum.accumulation_entropy_dir,
            "Accumulation: entropy directory is required when Accumulation is enabled.",
        )
        if not disambig.enabled and not accum.accumulation_caas_input:
            errors.append(
                "Accumulation is enabled but Disambiguation is disabled with no "
                "accumulation_caas_input fallback supplied."
            )
    else:
        if not any([accum.accumulation_caas_input, accum.accumulation_background_input]):
            errors.append(
                "Accumulation is disabled but no fallback input was supplied — Enrichment's "
                "accumulation gene lists may have no input."
            )

    # --- RERconverge ---
    rer = project.modules.rer
    if rer.enabled:
        require(rer.gene_trees, "RERconverge: gene trees file is required when RER is enabled.")

    # --- FADE ---
    # fade_mode always has a default; nothing else strictly required.

    # --- VEP ---
    vep = project.modules.vep
    if vep.enabled:
        require(vep.vep_map_dir, "VEP: per-gene MAP directory is required when VEP is enabled.")

    # --- Scoring ---
    scoring = project.modules.scoring
    if scoring.enabled:
        require(
            scoring.gene_ensembl_file,
            "Scoring: gene-Ensembl mapping file is required when Scoring is enabled.",
        )
        if not disambig.ct_postproc and not scoring.scoring_postproc_input:
            errors.append(
                "Scoring is enabled but Post-processing is disabled with no "
                "scoring_postproc_input fallback supplied."
            )
        if not rer.enabled and not all(row.scoring_rer_input for row in rt.phenotype_rows):
            errors.append(
                "Scoring is enabled and RERconverge is disabled — every phenotype row needs a "
                "scoring_rer_input fallback (some rows are missing one)."
            )
        if not project.modules.fade.enabled and not all(
            row.scoring_fade_summary_top and row.scoring_fade_summary_bottom
            for row in rt.phenotype_rows
        ):
            errors.append(
                "Scoring is enabled and FADE is disabled — every phenotype row needs "
                "scoring_fade_summary_top/bottom fallbacks (some rows are missing one)."
            )

    # --- Enrichment (+ POSENRICH) ---
    enrichment = project.modules.enrichment
    if enrichment.enabled:
        require(enrichment.gmt_dir, "Enrichment: GMT directory is required when Enrichment is enabled.")
        if not scoring.enabled and not any(
            [enrichment.enrichment_gene_lists_input, enrichment.enrichment_background_input]
        ):
            errors.append(
                "Enrichment is enabled but Scoring is disabled with no gene-list/background "
                "fallback supplied."
            )
    if enrichment.posenrich_enabled:
        for field_name, label in [
            ("egg_members_file", "eggNOG members file"),
            ("egg_annotations_file", "eggNOG annotations file"),
            ("cosmic_db", "COSMIC database"),
            ("domain_variability_file", "domain variability file"),
            ("ucr_positions_file", "UCR positions file"),
            ("fubar_sites_file", "FUBAR sites file"),
        ]:
            require(getattr(enrichment, field_name), f"POSENRICH: {label} is required.")

    return errors


def path_entries(project: ProjectConfig) -> list[tuple[str, str, str]]:
    """(label, path, kind) for every path-like field in the project. kind is
    "file" or "dir". Empty values are the caller's job to skip — validate() above
    already reports missing-but-required fields; this only checks paths that ARE
    filled in but don't exist on disk.

    Shared with gui/remote.py's SSH-based existence check, so local and remote
    validation always agree on exactly which fields count as paths.
    """
    general = project.general
    runtime = project.runtime
    m = project.modules

    entries: list[tuple[str, str, str]] = [
        ("General: repo directory", general.repo_dir, "dir"),
        ("General: Nextflow plugins directory", general.nextflow_plugins_dir, "dir"),
        ("Runtime: work directory", runtime.work_dir, "dir"),
        ("Runtime: results directory", runtime.results_dir, "dir"),
        ("Runtime: alignment directory", runtime.alignment_dir, "dir"),
        ("Runtime: species tree", runtime.tree_file, "file"),
        ("Runtime: CLASS 1 trait file", runtime.trait_file, "file"),
        ("Runtime: CLASS 2 trait file", runtime.simple_trait_file, "file"),
        ("Runtime: prune directory", runtime.prune_dir, "dir"),
        ("Runtime: alignment species names", runtime.ali_sp_names, "file"),
        ("Runtime: taxonomy ID mapping", runtime.tax_id_file, "file"),
        ("CAAS: config file", m.caas.caas_config_path, "file"),
        ("CAAS: discovery_from", m.caas.discovery_from, "dir"),
        ("CAAS: resample_from", m.caas.resample_from, "dir"),
        ("CAAS: bootstrap_from", m.caas.bootstrap_from, "dir"),
        ("Disambiguation: ASR cache directory", m.disambiguation.ct_disambig_asr_cache_dir, "dir"),
        ("Disambiguation: signification_from", m.disambiguation.signification_from, "dir"),
        ("Disambiguation: disambiguation_input", m.disambiguation.disambiguation_input, "file"),
        ("Disambiguation: disambiguation_dir", m.disambiguation.disambiguation_dir, "dir"),
        ("Disambiguation: background_input", m.disambiguation.background_input, "file"),
        ("Accumulation: entropy directory", m.accumulation.accumulation_entropy_dir, "dir"),
        ("Accumulation: accumulation_caas_input", m.accumulation.accumulation_caas_input, "file"),
        (
            "Accumulation: accumulation_background_input",
            m.accumulation.accumulation_background_input,
            "file",
        ),
        ("RERconverge: gene trees", m.rer.gene_trees, "file"),
        ("VEP: PrimateAI-3D database", m.vep.vep_primateai_db, "file"),
        ("VEP: MAP directory", m.vep.vep_map_dir, "dir"),
        ("VEP: COSMIC database", m.vep.scoring_vep_cosmic, "file"),
        ("VEP: vep_caas_input", m.vep.vep_caas_input, "file"),
        ("Scoring: gene-Ensembl file", m.scoring.gene_ensembl_file, "file"),
        ("Scoring: FADE site fallback (top)", m.scoring.scoring_fade_site_top, "file"),
        ("Scoring: FADE site fallback (bottom)", m.scoring.scoring_fade_site_bottom, "file"),
        ("Scoring: scoring_postproc_input", m.scoring.scoring_postproc_input, "file"),
        ("Scoring: scoring_accum_dir", m.scoring.scoring_accum_dir, "dir"),
        ("Scoring: scoring_vep_primateai", m.scoring.scoring_vep_primateai, "file"),
        ("Scoring: scoring_background_input", m.scoring.scoring_background_input, "file"),
        ("Scoring: caas_perms_file", m.scoring.caas_perms_file, "file"),
        ("Enrichment: GMT directory", m.enrichment.gmt_dir, "dir"),
        ("Enrichment: STRING database directory", m.enrichment.string_db_dir, "dir"),
        ("Enrichment: eggNOG members file", m.enrichment.egg_members_file, "file"),
        ("Enrichment: eggNOG annotations file", m.enrichment.egg_annotations_file, "file"),
        ("Enrichment: COSMIC database", m.enrichment.cosmic_db, "file"),
        ("Enrichment: domain variability file", m.enrichment.domain_variability_file, "file"),
        ("Enrichment: UCR positions file", m.enrichment.ucr_positions_file, "file"),
        ("Enrichment: FUBAR sites file", m.enrichment.fubar_sites_file, "file"),
        (
            "Enrichment: enrichment_background_input",
            m.enrichment.enrichment_background_input,
            "file",
        ),
        (
            "Enrichment: enrichment_gene_lists_input",
            m.enrichment.enrichment_gene_lists_input,
            "file",
        ),
        (
            "Enrichment: accumulation_enrichment_gene_lists_input",
            m.enrichment.accumulation_enrichment_gene_lists_input,
            "file",
        ),
        ("Enrichment: posenrich_background_file", m.enrichment.posenrich_background_file, "file"),
    ]

    for i, row in enumerate(runtime.phenotype_rows, start=1):
        if row.prune:
            entries.append((f"Phenotype row {i}: prune", os.path.join(runtime.prune_dir, row.prune), "file"))
        if row.prune_secondary:
            entries.append(
                (
                    f"Phenotype row {i}: prune_secondary",
                    os.path.join(runtime.prune_dir, row.prune_secondary),
                    "file",
                )
            )
        if row.scoring_rer_input:
            entries.append((f"Phenotype row {i}: scoring_rer_input", row.scoring_rer_input, "file"))
        if row.scoring_fade_summary_top:
            entries.append(
                (f"Phenotype row {i}: scoring_fade_summary_top", row.scoring_fade_summary_top, "file")
            )
        if row.scoring_fade_summary_bottom:
            entries.append(
                (
                    f"Phenotype row {i}: scoring_fade_summary_bottom",
                    row.scoring_fade_summary_bottom,
                    "file",
                )
            )

    return entries


def validate_paths(project: ProjectConfig) -> list[str]:
    """Check that every filled-in path field actually exists on disk.

    Complements validate() (which only checks required-ness) — a field can be
    non-empty but point at a typo'd or not-yet-mounted path, which validate()
    can't catch since it never touches the filesystem.
    """
    problems: list[str] = []
    for label, path, kind in path_entries(project):
        if not path.strip():
            continue  # required-ness is validate()'s job, not this function's
        exists = os.path.isfile(path) if kind == "file" else os.path.isdir(path)
        if not exists:
            noun = "file" if kind == "file" else "directory"
            problems.append(f"{label}: {noun} not found — {path}")
    return problems
