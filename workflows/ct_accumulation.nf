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

workflow CT_ACCUMULATION {
    take:
        caas_channel          // filtered_discovery.tsv from CT_POSTPROC gene_filtering (or Channel.empty())
        background_channel    // cleaned_background_main.txt from CT_POSTPROC  (or Channel.empty())
        trait_file_channel    // traitfile from CT (pruned when contrast_selection is on; or Channel.empty())

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
            background_ch
        )

        // ── Phase 2: Randomize ────────────────────────────────────────────────
        // global_csv already contains masked + iscaas (single-group accumulation)
        randomize_out = CT_ACCUMULATION_RANDOMIZE(
            aggregate_out.global_csv,
            meta_caas_val
        )

    emit:
        results    = randomize_out.results
}
