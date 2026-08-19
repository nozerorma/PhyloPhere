#!/usr/bin/env nextflow

/*
#
#  CT_ACCUMULATION Workflow: Tests for gene-level accumulation of CAAS using
#  permutation-based randomization.
#
#  Inputs (dual-mode: integrated pipeline OR standalone params):
#    - meta_caas_channel   : global_meta_caas.tsv from CT_SIGNIFICATION
#    - background_channel  : cleaned_background_main.txt from CT_POSTPROC
#    - trait_file_channel  : traitfile emitted by CT (pruned when contrast_selection is on)
#
#  Static inputs (always from params):
#    - params.alignment        : alignment directory (same as CT discovery)
#    - params.gene_ensembl_file: genomic info TSV
#    - params.caas_config      : traitfile fallback (3-col, no header: species trait pair)
*/

include { CT_ACCUMULATION_AGGREGATE; CT_ACCUMULATION_RANDOMIZE } from "${baseDir}/subworkflows/CT_ACCUMULATION/ctacc_run"
include { ACCUMULATION_REPORT } from "${baseDir}/subworkflows/CT_ACCUMULATION/accum_report"
include { ACCUMULATION_GENE_LISTS } from "${baseDir}/subworkflows/CT_ACCUMULATION/accum_gene_lists.nf"

workflow CT_ACCUMULATION {
    take:
        caas_channel          // filtered_discovery.tsv from CT_POSTPROC gene_filtering (or Channel.empty())
        background_channel    // cleaned_background_main.txt from CT_POSTPROC  (or Channel.empty())
        trait_file_channel    // traitfile from CT (pruned when contrast_selection is on; or Channel.empty())
        tested_positions_channel // caastools background.output — TESTED positions; the
                                 // randomization null's eligible pool (or Channel.empty())

    main:
        // ── Resolve CAAS source (filtered_discovery.tsv from postproc) ────────
        def meta_caas_ch = (caas_channel ?: Channel.empty())
        if (params.accumulation_caas_input) {
            meta_caas_ch = meta_caas_ch.ifEmpty {
                def meta_path = file(params.accumulation_caas_input)
                assert meta_path.exists() \
                    : "Error: accumulation_caas_input not found: ${params.accumulation_caas_input}"
                assert meta_path.isFile() \
                    : "Error: accumulation_caas_input must be a file: ${params.accumulation_caas_input}"

                log.info "📄 CT_ACCUMULATION: using standalone CAAS file: ${meta_path}"
                meta_path
            }
        }

        // Make it a value channel so it can be used in both AGGREGATE and RANDOMIZE
        // without using `.first()` (which warns when applied to a value channel).
        def meta_caas_val = meta_caas_ch
            .collect()
            .filter { files -> files && files.size() > 0 }
            .map { files ->
                files[0]
            }

        // ── Resolve background source ─────────────────────────────────────────
        def background_ch = (background_channel ?: Channel.empty())
        if (params.accumulation_background_input) {
            background_ch = background_ch.ifEmpty {
                def bg_path = file(params.accumulation_background_input)
                assert bg_path.exists() \
                    : "Error: accumulation_background_input not found: ${params.accumulation_background_input}"

                if (bg_path.isDirectory()) {
                    def preferred = file("${params.accumulation_background_input}/cleaned_background_main.txt")
                    if (preferred.exists()) {
                        log.info "📂 CT_ACCUMULATION: using cleaned_background_main.txt from directory"
                        preferred
                    } else {
                        def txts = bg_path.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
                        assert txts.size() > 0 \
                            : "Error: No .txt background files found in ${params.accumulation_background_input}"
                        log.info "📂 CT_ACCUMULATION: using background file: ${txts[0]}"
                        txts[0]
                    }
                } else {
                    log.info "📄 CT_ACCUMULATION: using background file: ${bg_path}"
                    bg_path
                }
            }
        }

        // Value channel so the cleaned background can feed BOTH AGGREGATE and
        // RANDOMIZE (as the null pool's gene universe). A queue channel can only be
        // consumed once — same reason meta_caas_val exists above.
        def background_val = background_ch
            .collect()
            .filter { files -> files && files.size() > 0 }
            .map { files -> files[0] }

        // ── Resolve species-list / traitfile ─────────────────────────────────
        def species_list_ch = (trait_file_channel ?: Channel.empty())
        if (!params.contrast_selection && params.caas_config) {
            species_list_ch = species_list_ch.ifEmpty {
                def tf_path = file(params.caas_config)
                assert tf_path.exists() \
                    : "Error: caas_config not found: ${params.caas_config}"

                log.info "📄 CT_ACCUMULATION: using standalone caas_config for species list: ${tf_path}"
                tf_path
            }
        }
        if (params.my_traits) {
            // Last-resort fallback. trait_file_channel only carries a value when CT
            // actually ran live (main.nf only populates ct_results when --ct_tool is
            // set); the moment CT/CONTRAST_SELECTION is fully precomputed-reused
            // instead — whether via the Precomputed Run tab's Discovery/Resample/
            // Bootstrap boxes, or a _complete pass reusing an _exploratory sibling's
            // CAAS output — trait_file_channel is empty. --caas_config above doesn't
            // help either: the GUI always drives CLASS 1/2 phenotypes through
            // --my_traits, never the legacy --caas_config path, so that condition is
            // never true for a GUI-generated run. Without this, species_list_ch never
            // emits and CT_ACCUMULATION_AGGREGATE (a process with a plain, non-optional
            // `path species_list` input) silently never executes — no error, just a
            // permanently pending process. params.my_traits is always required for
            // both CLASS 1 and CLASS 2 rows, so it's always safe to fall back to here.
            species_list_ch = species_list_ch.ifEmpty {
                def tf_path = file(params.my_traits)
                assert tf_path.exists() \
                    : "Error: my_traits not found: ${params.my_traits}"

                log.info "📄 CT_ACCUMULATION: using --my_traits for species list (CT precomputed-reused): ${tf_path}"
                tf_path
            }
        }

        def species_list_val = species_list_ch
            .collect()
            .filter { files -> files && files.size() > 0 }
            .map { files ->
                files[0]
            }

        // ── Static inputs (always from params) ───────────────────────────────
        assert params.alignment        : "CT_ACCUMULATION requires --alignment (alignment directory)"
        assert params.gene_ensembl_file: "CT_ACCUMULATION requires --gene_ensembl_file (genomic info TSV)"

        def alignment_dir_val = Channel.value(params.alignment.toString())
        def genomic_info_ch   = Channel.value(file(params.gene_ensembl_file))

        log.info "🔬 CT_ACCUMULATION"
        log.info "   Alignment dir  : ${params.alignment}"
        log.info "   Genomic info   : ${params.gene_ensembl_file}"
        log.info "   Rand type      : ${params.accumulation_randomization_type ?: 'naive'}"
        log.info "   Randomizations : ${params.accumulation_n_randomizations   ?: 10000}"

        // ── Phase 1: Aggregate ────────────────────────────────────────────────
        aggregate_out = CT_ACCUMULATION_AGGREGATE(
            alignment_dir_val,
            genomic_info_ch,
            species_list_val,
            meta_caas_val,
            background_val,
            params.accumulation_entropy_dir ?: ""
        )

        // ── Phase 2: Randomize — run once per phenotype direction + once for all ──
        // AGGREGATE is direction-agnostic (alignment/conservation data).
        // RANDOMIZE filters the CAAS pool by change_side and names outputs
        // accumulation_{direction}_{cat}_aggregated_results.csv.
        // "all" uses --change-side both (retains all non-none positions) for the
        // global accumulation score used in gene-level significance characterisation.
        // Eligible null pool inputs, broadcast to every direction. Sentinels (not
        // Channel.empty()) so a run without them still executes — randomize.py warns
        // and falls back to the ungapped-column pool rather than failing.
        def tested_pos_val = (tested_positions_channel ?: Channel.empty())
            .ifEmpty { file('NO_TESTED_POSITIONS') }
            .collect()
            .map { it[0] }
        def bg_universe_val = background_val.ifEmpty { file('NO_BG_UNIVERSE') }

        def rand_in = Channel.of("top", "bottom", "all")
            .combine(aggregate_out.global_csv)
            .combine(meta_caas_val)
            .combine(tested_pos_val)
            .combine(bg_universe_val)
            .multiMap { dir, global_csv, caas_csv, tested_pos, bg_universe ->
                direction: dir
                global:    global_csv
                caas:      caas_csv
                tested:    tested_pos
                universe:  bg_universe
            }

        randomize_out = CT_ACCUMULATION_RANDOMIZE(
            rand_in.direction,
            rand_in.global,
            rand_in.caas,
            rand_in.tested,
            rand_in.universe
        )

        // ── Accumulation report (rendered after all directions complete) ──────
        def all_rand_csvs = randomize_out.results.collect()
        def agg_csvs      = aggregate_out.global_csv.collect()
        ACCUMULATION_REPORT(all_rand_csvs, agg_csvs)

        // ── Gene list extraction (always automatic) ───────────────────────
        // Accumulation no longer runs its own FCS ranking — its FCS "significance"
        // was inflated by the missing permulation null (no compensating check against
        // the RER/CAAS-tuned FDR threshold). It now contributes as a cross-module
        // corroboration flag on CAAS's own leading edge instead. Its standalone AMI
        // report (DOMINO active modules) was removed too — never produced usable
        // output and Accumulation isn't part of the unified 13.AMI_analysis.Rmd
        // (only CAAS/FADE/RER are). There's no separate --ami flag anymore; gene
        // lists are simply always computed and published whenever Accumulation runs.
        ACCUMULATION_GENE_LISTS(randomize_out.direction, randomize_out.results)

    emit:
        results      = randomize_out.results
        accum_report = ACCUMULATION_REPORT.out.report
}
