#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: Functional Enrichment Workflow
 *
 * Runs gene-set and position-level functional enrichment analyses (FCS, AMI).
 */

include { SCORING_AMI_REPORT; SCORING_COMPARE_REPORT } from '../subworkflows/ENRICHMENT/scoring_enrichment.nf'
include { SCORING_FCS_REPORT }                            from '../subworkflows/ENRICHMENT/fcs.nf'
include { RER_FCS_REPORT  as MODULE_FCS_RER }             from '../subworkflows/ENRICHMENT/fcs.nf'
// CAAS calibrates per-gene rankings; both sources count pathway AUC by null block.
include { FCS_COMPUTE_CALIBRATED as FCS_COMPUTE_SCORING } from '../subworkflows/ENRICHMENT/fcs_blocks.nf'
include { FCS_RER as FCS_COMPUTE_RER } from '../subworkflows/ENRICHMENT/fcs_inputs.nf'
include { POSITION_ANNOTATIONS } from '../subworkflows/ENRICHMENT/position_annotations.nf'
include { POSENRICH_MODULAR } from '../subworkflows/ENRICHMENT/posenrich_modular.nf'
include { EGGNOG_RESOLUTION }                                from '../subworkflows/ENRICHMENT/eggnog_resolution.nf'
// Nextflow forbids invoking the same process/subworkflow more than once in one workflow
// scope without a distinct alias per call site (DuplicateProcessInvocation) -- DOMINO_MODULES
// is called up to 3 times below (main/CAAS, FADE, RER), each against a different
// gene-list/universe/background, so each needs its own imported name.
include { DOMINO_MODULES as DOMINO_MODULES_MAIN }          from '../subworkflows/ENRICHMENT/domino.nf'
include { DOMINO_MODULES as DOMINO_MODULES_FADE }          from '../subworkflows/ENRICHMENT/domino.nf'
include { DOMINO_MODULES as DOMINO_MODULES_RER }           from '../subworkflows/ENRICHMENT/domino.nf'

// ── FADE Global background/gene-list union ──────────────────────────────────
// FADE's AMI network is built once per trait (not once per direction) against
// the union of both directions' own universes -- they differ by only ~12
// genes out of ~16k on real data, so the union is a safe superset either way.
// The Global gene list (for the unified 13.AMI_analysis.Rmd's "Global" FADE
// tab) is likewise the union of both directions' significant-gene lists.
process FADE_UNION_BACKGROUND {
    label 'process_low'

    input:
    path top_bg,    stageAs: 'fade_top_background.txt'
    path bottom_bg, stageAs: 'fade_bottom_background.txt'

    output:
    path 'fade_global_background.txt', emit: background

    script:
    """
    cat fade_top_background.txt fade_bottom_background.txt | sort -u > fade_global_background.txt
    """
}

process FADE_UNION_SIGNIFICANT {
    label 'process_low'

    input:
    path top_sig,    stageAs: 'fade_top_significant.txt'
    path bottom_sig, stageAs: 'fade_bottom_significant.txt'

    output:
    path 'fade_global_significant.txt', emit: significant

    script:
    """
    cat fade_top_significant.txt fade_bottom_significant.txt | sort -u > fade_global_significant.txt
    """
}

// ── Shared "universes" directory ─────────────────────────────────────────────
// Each tool (CAAS/FADE/RER) has its own correct gene universe -- they legitimately
// differ (verified on real data: RER ~16.1k, FADE ~16.1k-per-direction, CAAS's
// own cleaned_background a different count again) and none should be derived
// from another. Published together, trait-suffixed, purely for easy browsing/
// comparison -- plain copies, no new logic. NO_-prefixed sentinel inputs (tool
// didn't run) are simply skipped.
process PUBLISH_UNIVERSES {
    label 'process_low'
    publishDir path: "${params.outdir}/universes", mode: 'copy', overwrite: true

    input:
    // stageAs gives each universe a distinct staged name -- RER, FADE (single-
    // direction runs pass their raw background.txt through) and CAAS-discovery
    // upstream processes all emit plain 'background.txt'/'cleaned_background_main.txt',
    // so unnamed staging collides whenever two tools ran in the same invocation.
    path rer_bg,  stageAs: 'rer_background.universe'
    path fade_bg, stageAs: 'fade_background.universe'
    path caas_bg, stageAs: 'caas_background.universe'
    val  traitname

    output:
    // Explicit dynamic filenames (not a "*.txt" glob) -- a glob here would also
    // match the Nextflow-staged INPUT files themselves (the *.universe files
    // staged above all sit in the task work dir), publishing confusing
    // un-renamed duplicates alongside the trait-suffixed copies below.
    path "rer_background_${traitname}.txt",             emit: rer_universe,  optional: true
    path "fade_background_${traitname}.txt",            emit: fade_universe, optional: true
    path "cleaned_background_main_${traitname}.txt",     emit: caas_universe, optional: true

    script:
    // stageAs renames every staged input, so the original NO_* sentinel name is
    // lost by the time the script block sees it -- a Groovy name check on the
    // staged variable can never match.  Sentinel inputs arrive as dangling
    // symlinks (file('NO_FADE_BACKGROUND') points to a file that does not
    // exist), so a shell -f test (which follows symlinks) is the reliable way
    // to skip them.
    """
    if [ -f rer_background.universe ];  then cp rer_background.universe  "rer_background_${traitname}.txt"; fi
    if [ -f fade_background.universe ]; then cp fade_background.universe "fade_background_${traitname}.txt"; fi
    if [ -f caas_background.universe ];  then cp caas_background.universe  "cleaned_background_main_${traitname}.txt"; fi
    """
}

// Channel of a `.collect()`-wrapped single-element List<Path> -> channel of the bare Path.
// Precomputed fallback of a VEP table: the file of the param when it exists, else the sentinel.
// A param that is set but points to no file is reported, since the sentinel hides it from the process.
def precomputedOrSentinel(String param_name, String param_value, String sentinel) {
    if (!param_value) return file(sentinel)
    if (file(param_value).exists()) return file(param_value)
    log.warn "${param_name} is set to '${param_value}', which does not exist: the table is treated as absent"
    return file(sentinel)
}

def unwrap1(ch) {
    return ch.map { it instanceof List ? it[0] : it }
}

// Validate supplied references before any expensive workflow tasks are submitted.
def validatePositionReferences() {
    def checkedAlignments = [] as Set
    def checkedMaps = [] as Set
    def hashFile = { path ->
        def digest = java.security.MessageDigest.getInstance('SHA-256')
        java.nio.file.Files.newInputStream(path).withCloseable { stream ->
            byte[] buffer = new byte[1024 * 1024]
            int count
            while ((count = stream.read(buffer)) > 0) digest.update(buffer, 0, count)
        }
        digest.digest().encodeHex().toString()
    }
    def sources = [
        [key: 'domain_variability_file', enabled: params.posenrich_domains, required: ['gene','pfam_id','target_name','description','clan_acc','clan_name','ali_start','ali_end','coordinate_space']],
        [key: 'ucr_positions_file', enabled: params.posenrich_ucr, required: ['gene','position','region_type','method','variability']],
        [key: 'fubar_sites_file', enabled: true, required: ['gene','is_pos_hit','is_neg_hit','hg38_aa_pos']],
        [key: 'posenrich_custom_file', enabled: true, required: []]
    ]
    sources.each { spec ->
        def source = params[spec.key]
        if (spec.enabled && source && !source.toString().startsWith('NO_')) {
            def input = file(source)
            if (!input.exists() && !source.toString().endsWith('.gz')) input = file("${source}.gz")
            if (!input.exists()) error "Supplied ${spec.key} does not exist: ${source}"
            if (spec.required) {
                def stream = java.nio.file.Files.newInputStream(input)
                if (input.toString().endsWith('.gz')) stream = new java.util.zip.GZIPInputStream(stream)
                def header = new BufferedReader(new InputStreamReader(stream, 'UTF-8')).withCloseable { it.readLine() }
                def columns = header ? header.split('\t').toList() : []
                def missing = spec.required.findAll { !columns.contains(it) }
                if (missing) {
                    if (spec.key == 'domain_variability_file' && missing.contains('coordinate_space'))
                        error "Pfam table '${input}' has no coordinate_space. Regenerate on the productive CAAS alignments; do not relabel legacy endpoints without verification."
                    error "${spec.key} '${input}' missing columns: ${missing}"
                }
            }
            def manifest = file("${input}.manifest.json")
            if (manifest.exists()) {
                def meta = new groovy.json.JsonSlurper().parseText(manifest.text)
                if (meta.table_sha256 && hashFile(input) != meta.table_sha256)
                    error "Reference checksum mismatch: ${input}"
                if (meta.alignment_dir && params.alignment) {
                    def expected = file(meta.alignment_dir).toRealPath()
                    def productive = file(params.alignment).toRealPath()
                    if (expected != productive) error "${spec.key} belongs to alignment '${expected}', but this run uses '${productive}'"
                    if (meta.alignments && !checkedAlignments.contains(productive.toString())) {
                        def actualNames = java.nio.file.Files.list(productive).withCloseable { paths ->
                            paths.filter { java.nio.file.Files.isRegularFile(it) && it.fileName.toString().toLowerCase() ==~ /.*\.(fa|fasta|faa|fas)/ }
                                .map { it.fileName.toString() }.collect(java.util.stream.Collectors.toSet())
                        }
                        if (actualNames != meta.alignments.keySet()) error "Alignment inventory changed since reference generation: ${productive}"
                        meta.alignments.each { name, expectedHash ->
                            def alignment = productive.resolve(name)
                            if (!java.nio.file.Files.exists(alignment) || hashFile(alignment) != expectedHash)
                                error "Productive alignment changed since reference generation: ${alignment}"
                        }
                        checkedAlignments.add(productive.toString())
                    }
                }
                if (meta.map_dir && params.caas_map_dir && file(meta.map_dir).toRealPath() != file(params.caas_map_dir).toRealPath())
                    error "${spec.key} was validated with a different MAP directory: ${meta.map_dir}"
                if (meta.map_dir && meta.maps && params.caas_map_dir && !checkedMaps.contains(meta.map_dir)) {
                    def mapRoot = file(params.caas_map_dir)
                    def actualMaps = java.nio.file.Files.list(mapRoot).withCloseable { paths ->
                        paths.filter { java.nio.file.Files.isRegularFile(it) && it.fileName.toString().endsWith('.map.tsv') }
                            .map { it.fileName.toString() }.collect(java.util.stream.Collectors.toSet())
                    }
                    if (actualMaps != meta.maps.keySet()) error "MAP inventory changed since reference generation: ${mapRoot}"
                    meta.maps.each { name, expectedHash ->
                        def mapPath = file(params.caas_map_dir).resolve(name)
                        if (!java.nio.file.Files.exists(mapPath) || hashFile(mapPath) != expectedHash)
                            error "Productive MAP changed since reference generation: ${mapPath}"
                    }
                    checkedMaps.add(meta.map_dir)
                }
            }
        }
    }
}

workflow ENRICHMENT {
    take:
        fcs_stats
        fcs_stats_rer
        fcs_stats_fade
        fcs_stats_accum
        gene_lists
        gene_scores
        cleaned_background_ch
        rer_perms_ch
        caas_perms_ch
        caas_blocks_ch
        gene_calibration
        caas_pos_sample_ch
        position_scores
        scoring_postproc_ch        // filtered_discovery.tsv (evaluated positions; POSENRICH modular registry)
        scoring_null_universe_ch   // null_universe_summary.tsv (positional null universe; POSENRICH modular registry)
        background_output_ch
        vep_primateai_ch     // optional: PrimateAI-3D score TSV (null when --vep not run)
        vep_cosmic_ch        // optional: COSMIC Mutant Census score TSV (null when --vep not run)
        vep_ensembl_ch       // optional: Ensembl VEP consequence TSV (null when --vep not run)
        fade_sites_top_ch    // optional: fade_sites_top.csv (null when --fade not run)
        fade_sites_bottom_ch // optional: fade_sites_bottom.csv (null when --fade not run)
        rer_gene_lists_bg_ch        // RER's own background.txt (empty when RER didn't run)
        rer_gene_lists_interest_ch  // rer_significant/accelerating/decelerating.txt
        fade_gene_lists_bg_top_ch     // FADE top direction's own background.txt
        fade_gene_lists_bg_bottom_ch  // FADE bottom direction's own background.txt
        fade_gene_lists_sig_top_ch    // fade_top_significant.txt
        fade_gene_lists_sig_bottom_ch // fade_bottom_significant.txt
        taxid_map_ch                  // tax_id map (curated by NAME_CURATION) for the clade variability of UCR, or NO_FILE

    main:
        validatePositionReferences()
        def fcs_universe_ch = (cleaned_background_ch ?: Channel.empty())
            .ifEmpty { file('NO_BACKGROUND') }
            .collect()

        // rerconverge.nf/fade.nf emit these already `.collect()`-wrapped (each a
        // single-element List<Path>) since they're built via a `.filter(...).collect()`
        // chain there; normalize to a single Path here so every downstream use
        // (stageAs, DOMINO_MODULES' background arg, PUBLISH_UNIVERSES) sees one
        // consistent shape regardless of whether the tool ran this invocation.
        def rer_bg_r         = unwrap1(rer_gene_lists_bg_ch)
        def fade_bg_top_r    = unwrap1(fade_gene_lists_bg_top_ch)
        def fade_bg_bottom_r = unwrap1(fade_gene_lists_bg_bottom_ch)
        def fade_sig_top_r    = unwrap1(fade_gene_lists_sig_top_ch)
        def fade_sig_bottom_r = unwrap1(fade_gene_lists_sig_bottom_ch)

        // RER has its own natural gene universe (every gene RERConverge actually
        // attempted), independent of the CAAS-discovery background. Priority:
        // explicit --rer_universe_file override > RER's own live-run background
        // (rer_bg_r, from RER_MAIN.out.gene_lists_bg -- already computed above for
        // PUBLISH_UNIVERSES/DOMINO but previously never reached the FCS report) >
        // CAAS's background as a last-resort fallback for runs where RER didn't
        // run this invocation and no override was given (preserves prior behaviour).
        def rer_ran = (params.rer_tool || params.rer_continuous_file || params.rer_binary_file) as boolean
        def rer_universe_ch  = params.rer_universe_file
            ? Channel.value(file(params.rer_universe_file))
            : (rer_ran ? rer_bg_r : fcs_universe_ch)

        def caas_perms_resolved = (caas_perms_ch ?: Channel.empty())
            .ifEmpty { file('NO_FILE') }
            .collect()
            .map { it[0] }

        def rer_perms_resolved = (rer_perms_ch ?: Channel.empty())
            .ifEmpty { file(params.rer_perms_file ?: 'NO_FILE') }
            .collect()
            .map { it[0] }

        // SCORING's own published nominal active gene_lists/ (scoring_evidence.R),
        // consumed by the DOMINO CAAS seed lists below.
        def gene_lists_ch = gene_lists.ifEmpty { file('NO_GENE_LISTS') }

        // FCS ranked-Wilcoxon enrichment (gated by params.fcs_enabled)
        def run_fcs = (params.fcs_enabled != null) ? params.fcs_enabled : true
        def fcs_out_report = Channel.empty()
        def rer_fcs_report = Channel.empty()
        def caas_all = Channel.value(file('NO_FCS_ALL'))
        def rer_all  = Channel.value(file('NO_FCS_ALL'))
        def caas_le  = Channel.value(file('NO_LEADING_EDGE'))
        def rer_le   = Channel.value(file('NO_LEADING_EDGE'))
        def caas_le_comp = Channel.value(file('NO_LEADING_EDGE_COMPOSITION'))
        def rer_le_comp  = Channel.value(file('NO_LEADING_EDGE_COMPOSITION'))

        def annot_ch = fcs_stats
        def trait_lbl = params.traitname ?: 'trait'

        if (run_fcs) {
            def scoring_fcs_compute = FCS_COMPUTE_SCORING(fcs_universe_ch, gene_calibration)
            def fcs_out = SCORING_FCS_REPORT(fcs_stats, fcs_universe_ch,
                                              scoring_fcs_compute.enrich_file)
            fcs_out_report = fcs_out.report
            caas_all = fcs_out.fcs_all_results.ifEmpty { file('NO_FCS_ALL') }
            caas_le  = fcs_out.fcs_leading_edge.ifEmpty { file('NO_LEADING_EDGE') }
            caas_le_comp = fcs_out.fcs_leading_edge_composition.ifEmpty { file('NO_LEADING_EDGE_COMPOSITION') }

            if (rer_ran) {
                def rer_fcs_compute = FCS_COMPUTE_RER(fcs_stats_rer, rer_universe_ch, rer_perms_resolved)
                def rer_fcs = MODULE_FCS_RER(
                    Channel.value('scoring/rer'),   fcs_stats_rer,
                    rer_universe_ch, Channel.value("12.FCS_rer_${trait_lbl}"),
                    rer_perms_resolved, annot_ch, rer_fcs_compute.enrich_file)
                rer_fcs_report = rer_fcs.report
                rer_all  = rer_fcs.fcs_all_results.ifEmpty  { file('NO_FCS_ALL') }
                rer_le   = rer_fcs.fcs_leading_edge.ifEmpty { file('NO_LEADING_EDGE') }
                rer_le_comp  = rer_fcs.fcs_leading_edge_composition.ifEmpty  { file('NO_LEADING_EDGE_COMPOSITION') }
            }
        }

        // FADE and RER each have their own gene universe (verified to differ
        // meaningfully from CAAS's cleaned_background and from each other on
        // real data) -- computed once here regardless of whether AMI itself
        // runs, since PUBLISH_UNIVERSES below wants it too. NO_-prefixed
        // sentinels flow through when a tool didn't run this invocation.
        // FADE "ran" for AMI purposes whenever its gene lists/background are
        // actually available -- either a live --fade this invocation, or a
        // precomputed summary TSV supplied via --fade_tables_dir_top/_bottom
        // (main.nf derives fade_gene_lists_bg/sig_* from that TSV the same way
        // it derives scoring_fade_top_ch/scoring_fade_bot_ch), mirroring how
        // rer_ran already treats rer_continuous_file as equivalent to rer_tool.
        def fade_ran = (params.fade || params.fade_tables_dir_top || params.fade_tables_dir_bottom) as boolean
        def fade_has_top    = params.fade ? ((params.fade_direction ?: 'both') in ['top', 'both']) : (params.fade_tables_dir_top as boolean)
        def fade_has_bottom = params.fade ? ((params.fade_direction ?: 'both') in ['bottom', 'both']) : (params.fade_tables_dir_bottom as boolean)
        // rer_ran already hoisted above (needed earlier for rer_universe_ch).

        def fade_union_bg_ch
        def fade_union_sig_ch
        if (fade_has_top && fade_has_bottom) {
            def fade_union_bg  = FADE_UNION_BACKGROUND(fade_bg_top_r, fade_bg_bottom_r)
            def fade_union_sig = FADE_UNION_SIGNIFICANT(fade_sig_top_r, fade_sig_bottom_r)
            fade_union_bg_ch  = fade_union_bg.background
            fade_union_sig_ch = fade_union_sig.significant
        } else if (fade_has_top) {
            fade_union_bg_ch  = fade_bg_top_r
            fade_union_sig_ch = fade_sig_top_r
        } else if (fade_has_bottom) {
            fade_union_bg_ch  = fade_bg_bottom_r
            fade_union_sig_ch = fade_sig_bottom_r
        } else {
            fade_union_bg_ch  = Channel.value(file('NO_FADE_BACKGROUND'))
            fade_union_sig_ch = Channel.value(file('NO_FADE_SIGNIFICANT'))
        }
        // Manual override, mirroring rer_universe_ch above -- unlike RER (which
        // falls back to CAAS's background when unset), FADE's own union background
        // computed above is already correct on its own; this only matters when you
        // want to substitute a hand-picked universe file instead.
        if (params.fade_universe_file) {
            fade_union_bg_ch = Channel.value(file(params.fade_universe_file))
        }

        def rer_universe_for_publish_ch = rer_ran ? rer_bg_r : Channel.value(file('NO_RER_BACKGROUND'))
        PUBLISH_UNIVERSES(rer_universe_for_publish_ch, fade_union_bg_ch, fcs_universe_ch, Channel.value(trait_lbl))

        // AMI active module identification (DOMINO + STRING labelling)
        def run_ami = params.scoring_ami ?: params.scoring_string ?: false
        def ami_out = null
        if (run_ami) {
            def domino_out = DOMINO_MODULES_MAIN(
                fcs_universe_ch,
                gene_lists_ch,
                params.domino_network_score_thr ?: 700,
                params.domino_slice_thr ?: 0.3,
                params.domino_module_thr ?: 0.05,
                params.string_db_dir ?: '',
                'caas'
            )

            // FADE and RER each get their own DOMINO network, built against
            // the universe resolved above -- never reuse CAAS's network for
            // these. Each is gated on whether that tool actually ran this
            // invocation; when it didn't, a NO_-prefixed sentinel flows
            // through so 13.AMI_analysis.Rmd's NULL-arg checks (mirroring the
            // existing gene_scores/gs_arg pattern) skip that tool's section
            // gracefully.
            def fade_gene_list_files_arg
            def fade_background_arg
            def fade_domino_sif_arg
            def fade_domino_modules_arg
            def fade_domino_edges_arg
            if (fade_ran) {
                def fade_sig_mix = Channel.empty()
                if (fade_has_top)    fade_sig_mix = fade_sig_mix.mix(fade_sig_top_r)
                if (fade_has_bottom) fade_sig_mix = fade_sig_mix.mix(fade_sig_bottom_r)
                if (fade_has_top && fade_has_bottom) fade_sig_mix = fade_sig_mix.mix(fade_union_sig_ch)
                def fade_gene_lists_final_ch = fade_sig_mix.collect()
                def fade_domino_out = DOMINO_MODULES_FADE(
                    fade_union_bg_ch,
                    fade_gene_lists_final_ch,
                    params.domino_network_score_thr ?: 700,
                    params.domino_slice_thr ?: 0.3,
                    params.domino_module_thr ?: 0.05,
                    params.string_db_dir ?: '',
                    'fade'
                )
                fade_gene_list_files_arg = fade_gene_lists_final_ch
                fade_background_arg      = fade_union_bg_ch
                fade_domino_sif_arg      = fade_domino_out.network_sif
                fade_domino_modules_arg  = fade_domino_out.modules_dir
                fade_domino_edges_arg    = fade_domino_out.edge_scores
            } else {
                fade_gene_list_files_arg = Channel.value(file('NO_FADE_GENE_LISTS'))
                fade_background_arg      = Channel.value(file('NO_FADE_BACKGROUND'))
                fade_domino_sif_arg      = Channel.value(file('NO_FADE_NETWORK_SIF'))
                fade_domino_modules_arg  = Channel.value(file('NO_FADE_MODULES_DIR'))
                fade_domino_edges_arg    = Channel.value(file('NO_FADE_EDGE_SCORES'))
            }

            def rer_gene_list_files_arg
            def rer_background_arg
            def rer_domino_sif_arg
            def rer_domino_modules_arg
            def rer_domino_edges_arg
            if (rer_ran) {
                def rer_domino_out = DOMINO_MODULES_RER(
                    rer_bg_r,
                    rer_gene_lists_interest_ch,
                    params.domino_network_score_thr ?: 700,
                    params.domino_slice_thr ?: 0.3,
                    params.domino_module_thr ?: 0.05,
                    params.string_db_dir ?: '',
                    'rer'
                )
                rer_gene_list_files_arg = rer_gene_lists_interest_ch
                rer_background_arg      = rer_bg_r
                rer_domino_sif_arg     = rer_domino_out.network_sif
                rer_domino_modules_arg = rer_domino_out.modules_dir
                rer_domino_edges_arg   = rer_domino_out.edge_scores
            } else {
                rer_gene_list_files_arg = Channel.value(file('NO_RER_GENE_LISTS'))
                rer_background_arg      = Channel.value(file('NO_RER_BACKGROUND'))
                rer_domino_sif_arg     = Channel.value(file('NO_RER_NETWORK_SIF'))
                rer_domino_modules_arg = Channel.value(file('NO_RER_MODULES_DIR'))
                rer_domino_edges_arg   = Channel.value(file('NO_RER_EDGE_SCORES'))
            }

            ami_out = SCORING_AMI_REPORT(
                gene_lists_ch, fcs_universe_ch, gene_scores,
                domino_out.network_sif, domino_out.modules_dir, domino_out.edge_scores,
                fade_gene_list_files_arg, fade_background_arg,
                fade_domino_sif_arg, fade_domino_modules_arg, fade_domino_edges_arg,
                rer_gene_list_files_arg, rer_background_arg,
                rer_domino_sif_arg, rer_domino_modules_arg, rer_domino_edges_arg,
                Channel.value(fade_ran), Channel.value(rer_ran)
            )
        }

        def final_reports = fcs_out_report.mix(rer_fcs_report)
        if (run_ami) {
            final_reports = final_reports.mix(ami_out.report)
        }

        // ── Position-level FCS ──
        // POSENRICH took its p.perm from the position-level null as a table of every (position, cycle). The null is kept as
        // counters per block of cycles, with no such table, so the module is not run (docs/POSENRICH_PORTE_EVALUACION.md
        // describes how to port it to the counters).
        def posenrich_out = null
        if (params.posenrich) {
            log.warn "[POSENRICH] skipped: the position-level permulation null is no longer written as a table of every (position, cycle), which its p.perm needs. The position-level enrichment (report 14) is not produced."
        }

        // ── Comparison report — pulls in CAAS/RER FCS, plus (when available) AMI's
        // module descriptions and posenrich's position-level results, all called
        // above so their outputs exist as real objects by this point in the DAG. ──

        def ami_module_desc_ch = (run_ami ? ami_out.module_descriptions : Channel.empty())
            .ifEmpty { file('NO_AMI_MODULE_DESC') }
        def ami_term_membership_ch = (run_ami ? ami_out.term_membership : Channel.empty())
            .ifEmpty { file('NO_AMI_TERM_MEMBERSHIP') }
        def posenrich_dotplot_ch = (posenrich_out ? posenrich_out.overall_dotplot : Channel.empty())
            .ifEmpty { file('NO_POSENRICH_DOTPLOT') }
        def posenrich_le_ch = (posenrich_out ? posenrich_out.leading_edge : Channel.empty())
            .ifEmpty { file('NO_POSENRICH_LE') }
        def posenrich_le_summary_ch = (posenrich_out ? posenrich_out.leading_edge_summary : Channel.empty())
            .ifEmpty { file('NO_POSENRICH_LE_SUMMARY') }
        // Interesting Genes/Positions tables reuse the SAME optional gene/position-
        // level inputs POSENRICH's own report already consumes (fcs_stats,
        // position_scores, VEP, FADE sites) -- all workflow-level take: params, so
        // available here whether or not POSENRICH itself ran this invocation.
        // Reuses annot_ch (fcs_stats.first(), already established above as the
        // broadcastable form of fcs_stats -- see its own definition) rather than
        // touching the raw fcs_stats channel a further time.
        def cmp_fcs_stats_ch  = annot_ch ?: file('NO_FCS_STATS')
        def cmp_pos_scores_ch = position_scores ? position_scores : file('NO_POS_SCORES')
        // Same precomputed-VEP fallback as pos_vep_primateai_ch/pos_vep_cosmic_ch above.
        def cmp_vep_pai_ch     = (vep_primateai_ch ?: Channel.empty())
            .ifEmpty { (params.scoring_vep_primateai && file(params.scoring_vep_primateai).exists()) ? file(params.scoring_vep_primateai) : file('NO_VEP_PAI') }
        def cmp_vep_cosmic_ch  = (vep_cosmic_ch ?: Channel.empty())
            .ifEmpty { precomputedOrSentinel('scoring_vep_cosmic', params.scoring_vep_cosmic, 'NO_VEP_COSMIC') }
        def cmp_vep_ensembl_ch = (vep_ensembl_ch ?: Channel.empty())
            .ifEmpty { precomputedOrSentinel('scoring_vep_ensembl', params.scoring_vep_ensembl, 'NO_VEP_ENSEMBL') }
        def cmp_fade_top_ch   = (fade_sites_top_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_TOP') }
        def cmp_fade_bot_ch   = (fade_sites_bottom_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_BOTTOM') }
        // Annotations do not require a positional permulation null. Use supplied
        // Pfam/UCR tables only: this path never starts downloads or heavy fitting.
        def posenrich_position_char_ch = Channel.empty()
        if (!params.posenrich_modular) {
            POSITION_ANNOTATIONS(
                Channel.value(precomputedOrSentinel('gene_ensembl_file', params.gene_ensembl_file, 'NO_ANNOT_ENSEMBL')),
                Channel.value(precomputedOrSentinel('domain_variability_file', params.posenrich_domains ? params.domain_variability_file : '', 'NO_ANNOT_DOMAIN')),
                Channel.value(precomputedOrSentinel('ucr_positions_file', params.posenrich_ucr ? params.ucr_positions_file : '', 'NO_ANNOT_UCR')),
                Channel.value(precomputedOrSentinel('fubar_sites_file', params.fubar_sites_file, 'NO_ANNOT_FUBAR')),
                Channel.value(precomputedOrSentinel('caas_map_dir', params.caas_map_dir, 'NO_ANNOT_MAP')),
                fcs_universe_ch.map { it[0] },
                cmp_fade_top_ch, cmp_fade_bot_ch,
                Channel.value(precomputedOrSentinel('posenrich_custom_file', params.posenrich_custom_file, 'NO_ANNOT_CUSTOM'))
            )
            posenrich_position_char_ch = POSITION_ANNOTATIONS.out.position_char
        }

        // ── POSENRICH modular (WP10) ──
        // Common registry + characterization + within-protein localization + global
        // association + proximity, gated by params.posenrich_modular. The historical
        // POSENRICH inference (Lachenbruch) stays omitted; no legacy null is reused.
        if (params.posenrich_modular) {
            def modular_out = POSENRICH_MODULAR(
                cmp_pos_scores_ch,
                (scoring_postproc_ch ?: Channel.empty()).ifEmpty { file('NO_POSTPROC') },
                (scoring_null_universe_ch ?: Channel.empty()).ifEmpty { file('NO_NULL_UNIVERSE') },
                Channel.value(precomputedOrSentinel('caas_map_dir', params.caas_map_dir, 'NO_ANNOT_MAP')),
                fcs_universe_ch.map { it[0] },
                Channel.value(precomputedOrSentinel('gene_ensembl_file', params.gene_ensembl_file, 'NO_ANNOT_ENSEMBL')),
                Channel.value(precomputedOrSentinel('domain_variability_file', params.posenrich_domains ? params.domain_variability_file : '', 'NO_ANNOT_DOMAIN')),
                Channel.value(precomputedOrSentinel('ucr_positions_file', params.posenrich_ucr ? params.ucr_positions_file : '', 'NO_ANNOT_UCR')),
                Channel.value(precomputedOrSentinel('fubar_sites_file', params.fubar_sites_file, 'NO_ANNOT_FUBAR')),
                cmp_fade_top_ch, cmp_fade_bot_ch,
                Channel.value(precomputedOrSentinel('posenrich_custom_file', params.posenrich_custom_file, 'NO_ANNOT_CUSTOM'))
            )
            posenrich_position_char_ch = modular_out.position_char
            final_reports = final_reports.mix(modular_out.report)
        }

        def cmp_out = SCORING_COMPARE_REPORT(caas_all, rer_all, caas_le, rer_le,
                                              caas_le_comp, rer_le_comp,
                                              ami_module_desc_ch, ami_term_membership_ch,
                                              posenrich_dotplot_ch, posenrich_le_ch,
                                              posenrich_le_summary_ch,
                                              cmp_fcs_stats_ch, cmp_pos_scores_ch,
                                              cmp_vep_pai_ch, cmp_vep_cosmic_ch,
                                              cmp_vep_ensembl_ch,
                                              cmp_fade_top_ch, cmp_fade_bot_ch,
                                              posenrich_position_char_ch)
        final_reports = final_reports.mix(cmp_out.report)

    emit:
        report = final_reports
}
