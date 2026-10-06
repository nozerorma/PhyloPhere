#!/usr/bin/env python3
# validate.py — Pre-render validation: required fields per enabled module, and path existence.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure functions, no PySide6 import. validate() returns a list of human-readable
error strings; an empty list means the project is renderable. The GUI's "Generate
Scripts..." action shows them in a dialog, so a missing field is reported there
instead of surfacing as a StrictUndefined error inside a template.

The scope is the fields exposed as GUI widgets, not the full conf/*.config
parameter space. validate_paths() and path_entries() cover the complementary
question: whether the filled-in paths exist.

Imported by: gui/widgets/main_window.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.precomputed import derive_paths
from gui.models.project import ProjectConfig


def validate(project: ProjectConfig) -> list[str]:
    """Return the error messages of a project that cannot be rendered (empty if valid).

    Each enabled module is checked only for inputs that have no fallback: a blank
    field that the pipeline can auto-generate is not an error.
    """
    errors: list[str] = []

    def require(value: str, message: str) -> None:
        if not value.strip():
            errors.append(message)

    # ── General ───────────────────────────────────────────────────────────────
    require(project.general.repo_dir, "General: PhyloPhere repo directory is required.")
    require(
        project.general.nextflow_plugins_dir,
        "General: Nextflow plugins directory is required (symlinked into every run's NXF_HOME).",
    )

    # ── Runtime ───────────────────────────────────────────────────────────────
    rt = project.runtime
    pc = project.precomputed
    require(rt.alignment_dir, "Runtime: alignment directory is required.")
    require(rt.tree_file, "Runtime: species tree file is required.")
    require(rt.work_dir, "Runtime: work directory is required.")
    require(rt.results_dir, "Runtime: results directory is required.")

    if rt.toy_mode and not (rt.toy_perms.strip().isascii() and rt.toy_perms.strip().isdigit() and int(rt.toy_perms) >= 0):
        errors.append(f"Runtime: toy permulation cycles must be a non-negative integer (0 replays only the real labeling) (got {rt.toy_perms!r}).")

    if rt.toy_mode and not (rt.toy_n.strip().isascii() and rt.toy_n.strip().isdigit() and int(rt.toy_n) > 0):
        errors.append(f"Runtime: the toy sample size must be a positive integer; a blank or zero value makes the pipeline sample 50 alignments (got {rt.toy_n!r}).")

    if not rt.phenotype_rows:
        errors.append("Runtime: the phenotype catalogue must have at least one row.")

    require(rt.trait_file, "Runtime: trait file is required.")
    any_pruned = any(row.prune.strip() or row.prune_secondary.strip() for row in rt.phenotype_rows)
    if any_pruned:
        require(rt.prune_dir, "Runtime: prune directory is required (a phenotype row has PRUNE/PRUNE_SEC set).")

    for i, row in enumerate(rt.phenotype_rows, start=1):
        require(row.trait, f"Phenotype row {i}: trait name is required.")
        if str(row.trait_type).strip().lower() not in ("", "auto", "ordinal", "continuous"):
            errors.append(
                f"Phenotype row {i}: TRAIT_TYPE must be blank, auto, ordinal, or continuous "
                f"(got {row.trait_type!r})."
            )

    # ── Precomputed Run: base_path is required whenever a reuse box is checked 
    any_precomp_checked = any(
        [
            pc.use_discovery, pc.use_resample, pc.use_disambiguation,
            pc.use_postproc, pc.use_accumulation, pc.use_rer, pc.use_fade, pc.use_vep,
        ]
    )
    if any_precomp_checked:
        require(
            pc.base_path,
            "Precomputed Run: base path is required — at least one reuse checkbox is checked.",
        )

    # ── CAAS ──────────────────────────────────────────────────────────────────
    # caas_config_path is never required. CONTRAST_SELECTION supplies the trait file
    # and run_single.sh.j2 turns it on whenever CAAS or Disambiguation is enabled, so
    # a trait file exists for whichever consumer needs one:
    #   * CAAS on: the trait and tree of CONTRAST_SELECTION feed CT(...) directly.
    #   * CAAS off, Disambiguation on (reuse via the Precomputed Run tab): CT() does
    #     not run (ct_tool is empty) but CONTRAST_SELECTION does, and main.nf hands
    #     its trait and tree files to the observed scoring (CT_OBSERVED).
    # The --caas_config fallbacks are therefore reached only by standalone, non-GUI
    # invocations, which is why the GUI can drive phenotypes through --my_traits alone.
    #
    # The output of CAAS (discovery, resample) is consumed downstream only by
    # Disambiguation, so a disabled CAAS is a problem only when Disambiguation is
    # enabled and nothing is reused in its place.
    caas = project.modules.caas
    disambig = project.modules.disambiguation
    if not caas.enabled and disambig.enabled:
        if not any([pc.use_discovery, pc.use_resample]):
            errors.append(
                "CAAS is disabled but no Discovery/Resample reuse box is checked on "
                "the Precomputed Run tab — downstream modules (Disambiguation, Accumulation) "
                "have no input."
            )

    # ── Disambiguation (+ Post-processing) ────────────────────────────────────
    if disambig.enabled:
        require(
            disambig.ct_disambig_asr_cache_dir,
            "Disambiguation: the ASR cache directory is required.",
        )
        if not caas.enabled and not any([pc.use_discovery, pc.use_resample]):
            errors.append(
                "Disambiguation is enabled but CAAS is disabled with no reuse box checked on "
                "the Precomputed Run tab."
            )
        # Post-processing always runs with Disambiguation (it has no toggle of its
        # own, see DisambiguationConfig in gui/models/modules.py). Its gene filtering
        # and characterization reports need gene_ensembl_file, which is not required
        # here: left blank, it is generated from the alignment gene list by an Ensembl
        # BioMart query (bin/generate_ensembl_mapping.py, called through
        # bin/resolve_core_inputs.py). The alignment directory is required above, so
        # the source of that gene list is always present.
    # With Disambiguation disabled, Post-processing is disabled too. The Accumulation
    # and Scoring checks below cover the case that still needs its output through
    # pc.use_postproc.

    # ── Accumulation ──────────────────────────────────────────────────────────
    accum = project.modules.accumulation
    scoring = project.modules.scoring
    if accum.enabled:
        # accumulation_entropy_dir is not required: left blank, Valdar variability
        # files are generated from the alignment (bin/compute_alignment_entropy.py).
        # That needs --tax_id, itself optional because it is generated from the
        # tree; without a tax_id, Accumulation falls back to a coarser conservation
        # measure computed from the alignment (raw majority-residue conservation).
        if not disambig.enabled and not pc.use_postproc:
            errors.append(
                "Accumulation is enabled but Post-processing is disabled with no "
                "'Use precomputed Post-processing output' box checked on the Precomputed Run tab."
            )
    # Scoring does not require Accumulation: with no accumulation channel, main.nf
    # passes null, scoring.nf resolves a NO_ACCUM sentinel and scoring_compute.R
    # skips every accumulation code path (has_accum_dir). The only hard upstream of
    # Scoring is CT post-processing, checked below.

    # ── RERconverge ───────────────────────────────────────────────────────────
    rer = project.modules.rer
    if rer.enabled:
        require(rer.gene_trees, "RERconverge: gene trees file is required when RER is enabled.")

    # ── FADE ──────────────────────────────────────────────────────────────────
    # Every FADE parameter has a default; nothing is required.

    # ── VEP ───────────────────────────────────────────────────────────────────
    vep = project.modules.vep
    if vep.enabled:
        require(project.modules.disambiguation.caas_map_dir,
                "VEP: per-gene MAP directory (Disambiguation, post-processing) is required when VEP is enabled.")
        # vep_cache_dir is not required even with vep_ensembl checked: left blank,
        # ENSEMBL_VEP_ANNOTATE resolves a persistent default location and fills it
        # with vep_install on first use (subworkflows/VEP/ensembl_vep.nf).

    # ── Scoring ───────────────────────────────────────────────────────────────
    if scoring.enabled:
        # gene_ensembl_file is not required, as in the Disambiguation section: left
        # blank, it is generated from the alignment gene list by an Ensembl BioMart
        # query (bin/generate_ensembl_mapping.py), whether or not Disambiguation is
        # enabled.
        if not disambig.enabled and not pc.use_postproc:
            errors.append(
                "Scoring is enabled but Post-processing is disabled with no 'Use precomputed "
                "Post-processing output' box checked on the Precomputed Run tab."
            )
        # RER and FADE are optional inputs to Scoring. When either is absent, main.nf
        # passes null, the Scoring process resolves a NO_* sentinel and the
        # file_exists() guard of scoring_compute.R (it rejects any basename starting
        # with NO_) sets has_rer / has_fade, which skip every RER/FADE code path.
        # A CAAS-only Scoring run is valid and is not blocked here.

    # ── Evidence of the N best positions (CAAS_EVIDENCE, after SCORING) ───────
    n_evidence = scoring.caas_evidence_top_n.strip()
    if not (n_evidence.isascii() and n_evidence.isdigit()):
        errors.append(f"Scoring: the number of positions to explain (evidence) must be a non-negative integer (got {scoring.caas_evidence_top_n!r}).")
    elif int(n_evidence) > 0:
        if not scoring.enabled:
            errors.append("Scoring: evidence of the best positions explains position_scores.tsv, so Scoring must be enabled.")
        # main.nf re-scores the rows of those positions from the observed
        # discovery.tab: the run's own (CAAS feeding Disambiguation) or the reused
        # one (Precomputed Run, discovery_from).
        if not disambig.enabled or not (caas.enabled or pc.use_discovery):
            errors.append("Scoring: evidence of the best positions needs the observed discovery.tab: enable CAAS and "
                          "Disambiguation, or reuse a discovery on the Precomputed Run tab with Disambiguation enabled.")

    # ── Enrichment (+ POSENRICH) ──────────────────────────────────────────────
    enrichment = project.modules.enrichment
    # gmt_dir is not required: left blank, it defaults to the curated gene sets of
    # subworkflows/ENRICHMENT/dat/ (plus the downloaded sets when auto_fetch_gmt is
    # on, through bin/resolve_gmts.py).
    if enrichment.posenrich_enabled:
        # Nothing is required for POSENRICH. cosmic_db and fubar_sites_file fall back
        # to NO_FILE sentinels (workflows/enrichment.nf, posenrich.nf) and their
        # layers are skipped when absent. egg_members_file and egg_annotations_file,
        # both blank, use the pair versioned in subworkflows/ENRICHMENT/dat/, or a
        # download when auto_fetch_eggnog is set (eggnog_resolution.nf).
        # ucr_positions_file and domain_variability_file, when blank, are generated
        # from the alignment (ucr_generation.nf, domain_variability_generation.nf);
        # ucr_positions_file also needs runtime.tax_id, itself optional.
        pass

    return errors


def path_entries(project: ProjectConfig) -> list[tuple[str, str, str]]:
    """List (label, path, kind) for every path-like field of the project.

    kind is "file" or "dir". Empty values are included and left to the caller to
    skip: validate() reports the required ones that are missing. Shared with the
    SSH-based existence check of gui/remote.py, so local and remote validation
    cover the same fields.
    """
    general = project.general
    runtime = project.runtime
    m = project.modules
    pc = project.precomputed

    entries: list[tuple[str, str, str]] = [
        ("General: repo directory", general.repo_dir, "dir"),
        ("General: Nextflow plugins directory", general.nextflow_plugins_dir, "dir"),
        ("Runtime: work directory", runtime.work_dir, "dir"),
        ("Runtime: results directory", runtime.results_dir, "dir"),
        ("Runtime: alignment directory", runtime.alignment_dir, "dir"),
        ("Runtime: species tree", runtime.tree_file, "file"),
        ("Runtime: trait file", runtime.trait_file, "file"),
        ("Runtime: prune directory", runtime.prune_dir, "dir"),
        ("Runtime: alignment species names", runtime.ali_sp_names, "file"),
        ("Runtime: taxonomy ID mapping", runtime.tax_id_file, "file"),
        ("CAAS: config file", m.caas.caas_config_path, "file"),
        ("Disambiguation: ASR cache directory", m.disambiguation.ct_disambig_asr_cache_dir, "dir"),
        ("Accumulation: entropy directory", m.accumulation.accumulation_entropy_dir, "dir"),
        ("RERconverge: gene trees", m.rer.gene_trees, "file"),
        ("RERconverge: tested-gene universe file", m.rer.rer_universe_file, "file"),
        ("RERconverge: cross-module gene scores", m.rer.rer_gene_scores, "file"),
        ("FADE: custom fg/bg species file", m.fade.fade_species_file, "file"),
        ("FADE: Substitution matrix", m.fade.lg_dat_path, "file"),
        ("FADE: tested-gene universe file", m.fade.fade_universe_file, "file"),
        ("VEP: PrimateAI-3D database", m.vep.vep_primateai_db, "file"),
        ("Post-processing: MAP directory", m.disambiguation.caas_map_dir, "dir"),
        ("VEP: COSMIC database", m.vep.cosmic_db, "file"),
        ("VEP: Ensembl VEP cache directory", m.vep.vep_cache_dir, "dir"),
        ("Scoring: gene-Ensembl file", m.scoring.gene_ensembl_file, "file"),
        ("Scoring: contrast hypotheses pairs file", m.scoring.scoring_hypotheses_pairs, "file"),
        ("Disambiguation: contrast hypotheses pairs file", m.disambiguation.ct_disambig_hypotheses_pairs, "file"),
        ("Enrichment: GMT directory", m.enrichment.gmt_dir, "dir"),
        ("Enrichment: Pfam cache directory", m.enrichment.pfam_cache_dir, "dir"),
        ("Enrichment: STRING database directory", m.enrichment.string_db_dir, "dir"),
        ("Enrichment: eggNOG members file", m.enrichment.egg_members_file, "file"),
        ("Enrichment: eggNOG annotations file", m.enrichment.egg_annotations_file, "file"),
        ("Enrichment: domain variability file", m.enrichment.domain_variability_file, "file"),
        ("Enrichment: UCR positions file", m.enrichment.ucr_positions_file, "file"),
        ("Enrichment: FUBAR sites file", m.enrichment.fubar_sites_file, "file"),
    ]

    # Precomputed Run tab: the per-phenotype paths derived from base_path, not stored strings.
    for i, row in enumerate(runtime.phenotype_rows, start=1):
        if not row.trait:
            continue
        for name, path, kind in derive_paths(pc, row.trait):
            entries.append((f"Precomputed ({row.trait}): {name}", path, kind))

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

    return entries


def validate_paths(project: ProjectConfig) -> list[str]:
    """Return one message per filled-in path field that does not exist on disk.

    Complements validate(), which checks only that required fields are set and never
    touches the filesystem: a field can be non-empty and still point at a mistyped
    or unmounted path.
    """
    problems: list[str] = []
    for label, path, kind in path_entries(project):
        if not path.strip():
            continue  # a blank path is validate()'s concern
        exists = os.path.isfile(path) if kind == "file" else os.path.isdir(path)
        if not exists:
            noun = "file" if kind == "file" else "directory"
            problems.append(f"{label}: {noun} not found — {path}")
    return problems
