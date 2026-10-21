#!/usr/bin/env nextflow

/*
#
#  CT_ACCUMULATION Workflow: tests for gene-level accumulation of CAAS with an exact
#  hypergeometric test on the universe U of positions (the same U as SCORING and POSENRICH).
#
#  Inputs (dual-mode: integrated pipeline OR standalone params):
#    - caas_channel            : filtered_discovery.tsv from CT_POSTPROC (gene_filtering)
#    - background_channel      : cleaned_background_main.txt from CT_POSTPROC
#    - tested_positions_channel: caastools background.output (required)
#    - null_cycle_caas_channel : perm_pos_cycle_caas.tsv.gz from CAAS_CORE_MERGE (required, or the
#                                --caas_pos_cycle_caas_file fallback)
*/

include { CT_ACCUMULATION_TEST } from '../subworkflows/CT_ACCUMULATION/ctacc_run'
include { ACCUMULATION_REPORT } from '../subworkflows/CT_ACCUMULATION/accum_report'
include { ACCUMULATION_GENE_LISTS } from '../subworkflows/CT_ACCUMULATION/accum_gene_lists.nf'

workflow CT_ACCUMULATION {
    take:
        caas_channel             // filtered_discovery.tsv from CT_POSTPROC gene_filtering (or Channel.empty())
        background_channel       // cleaned_background_main.txt from CT_POSTPROC  (or Channel.empty())
        tested_positions_channel // caastools background.output — TESTED positions, the starting point of the
                                 // universe U (required: Channel.empty() makes the run stop with an error)
        null_cycle_caas_channel  // perm_pos_cycle_caas.tsv.gz from CAAS_CORE_MERGE, the filtered null SCORING also
                                 // uses; U keeps the tested positions it scores (required, or the
                                 // --caas_pos_cycle_caas_file fallback)

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

        // Make it a value channel so it feeds the three directions without `.first()`
        // (which warns when applied to a value channel).
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

        // Value channel so the cleaned background (the gene universe of U) feeds the three
        // directions. A queue channel can only be consumed once — same reason meta_caas_val exists above.
        def background_val = background_ch
            .collect()
            .filter { files -> files && files.size() > 0 }
            .map { files -> files[0] }

        log.info "🔬 CT_ACCUMULATION"

        // ── Exact test — run once per phenotype direction + once for all ──
        // The test filters the observed CAAS by direction and names outputs
        // accumulation_{direction}_results.csv.
        // "all" uses --change-side both (a position detected on both sides counts once) for the
        // global accumulation score used in gene-level significance characterisation.
        // Universe U inputs, broadcast to every direction. All are required: a missing one stops the run.
        def tested_pos_val = (tested_positions_channel ?: Channel.empty())
            .ifEmpty { error "CT_ACCUMULATION requires the tested positions (caastools background.output): run CAAS discovery or pass --posenrich_background_file." }
            .collect()
            .map { it[0] }
        def bg_universe_val = background_val
            .ifEmpty { error "CT_ACCUMULATION requires the cleaned background (cleaned_background_main.txt): run CT post-processing or pass --accumulation_background_input." }
        def null_cycle_caas_val = (null_cycle_caas_channel ?: Channel.empty())
            .ifEmpty { params.caas_pos_cycle_caas_file ? file(params.caas_pos_cycle_caas_file) : error("CT_ACCUMULATION requires the filtered permulation null (perm_pos_cycle_caas.tsv.gz): run CAAS permulation or pass --caas_pos_cycle_caas_file.") }
            .collect()
            .map { it[0] }

        def test_in = Channel.of("top", "bottom", "all")
            .combine(meta_caas_val)
            .combine(tested_pos_val)
            .combine(bg_universe_val)
            .combine(null_cycle_caas_val)
            .multiMap { dir, caas_csv, tested_pos, bg_universe, null_cycle_caas ->
                direction:         dir
                caas:              caas_csv
                tested:            tested_pos
                universe:          bg_universe
                null_cycle_caas:   null_cycle_caas
            }

        test_out = CT_ACCUMULATION_TEST(
            test_in.direction,
            test_in.caas,
            test_in.tested,
            test_in.universe,
            test_in.null_cycle_caas
        )

        // ── Accumulation report (rendered after all directions complete) ──────
        ACCUMULATION_REPORT(test_out.results.collect())

        // ── Gene list extraction (always automatic) ───────────────────────
        // Accumulation no longer runs its own FCS ranking — its FCS "significance"
        // was inflated by the missing permulation null (no compensating check against
        // the RER/CAAS-tuned FDR threshold). It now contributes as a cross-module
        // corroboration flag on CAAS's own leading edge instead. Its standalone AMI
        // report (DOMINO active modules) was removed too — never produced usable
        // output and Accumulation isn't part of the unified 13.AMI_analysis.Rmd
        // (only CAAS/FADE/RER are). There's no separate --ami flag anymore; gene
        // lists are simply always computed and published whenever Accumulation runs.
        ACCUMULATION_GENE_LISTS(test_out.direction, test_out.results)

    emit:
        results      = test_out.results
        accum_report = ACCUMULATION_REPORT.out.report
}
