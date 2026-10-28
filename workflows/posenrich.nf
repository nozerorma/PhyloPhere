#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: POSENRICH Workflow
 *
 * Position-level enrichment of CAAS detections on the universe U of positions
 * (the same U as SCORING):
 *   - the universe table (universe.tsv.gz), the one universe of levels 1 and 2;
 *   - level 1: Accumulation, an exact hypergeometric test per gene and direction, on that table;
 *   - the position registry and its characterization, proximity and level 2 (concentration of the signal inside
 *     genes, calibrated with the null cycles) components (they need SCORING);
 *   - the Pfam and UCR tables the registry consumes, generated when none is supplied;
 *   - the report (14.POSENRICH_report).
 *
 * Data order: SCORING, POSENRICH, ENRICHMENT. ENRICHMENT reads position_characterization.tsv
 * as an optional input.
 */

include { ACCUMULATION_TEST } from '../subworkflows/POSENRICH/accumulation_test.nf'
include { DOMAIN_VARIABILITY_GENERATION } from '../subworkflows/POSENRICH/domain_variability_generation.nf'
include { UCR_GENERATION } from '../subworkflows/POSENRICH/ucr_generation.nf'
include { POSITION_REGISTRY; POSENRICH_UNIVERSE; POSENRICH_CHARACTERIZE; POSENRICH_WITHIN_GENE; POSENRICH_PROXIMITY; POSENRICH_REPORT } from '../subworkflows/POSENRICH/position_modules.nf'

// Precomputed fallback of an optional table: the file of the param when it exists, else the sentinel.
// A param that is set but points to no file is reported, since the sentinel hides it from the process.
def precomputedOrSentinel(String param_name, String param_value, String sentinel) {
    if (!param_value) return file(sentinel)
    if (file(param_value).exists()) return file(param_value)
    if (!param_value.endsWith('.gz') && file("${param_value}.gz").exists()) return file("${param_value}.gz")
    log.warn "${param_name} is set to '${param_value}', which does not exist: the table is treated as absent"
    return file(sentinel)
}

// The standalone cleaned background (--accumulation_background_input: a file or a directory holding
// cleaned_background_main.txt), or null.
def standaloneBackground() {
    if (!params.accumulation_background_input) return null
    def input = file(params.accumulation_background_input)
    return input.isDirectory() ? file("${input}/cleaned_background_main.txt") : input
}

// Validate supplied references before any expensive workflow tasks are submitted.
def validatePositionReferences() {
    if ((params.posenrich_domains && !params.domain_variability_file) ||
        (params.posenrich_ucr && !params.ucr_positions_file)) {
        if (!params.alignment || !java.nio.file.Files.isDirectory(file(params.alignment)))
            error "Automatic Pfam/UCR annotations require --alignment pointing to the productive protein FASTA directory"
    }
    def checkedAlignments = [] as Set
    def checkedMaps = [] as Set
    def hashFile = { path ->
        def digest = java.security.MessageDigest.getInstance('SHA-256')
        new java.security.DigestInputStream(java.nio.file.Files.newInputStream(path), digest).withCloseable { stream ->
            stream.transferTo(java.io.OutputStream.nullOutputStream())
        }
        digest.digest().encodeHex().toString()
    }
    def sources = [
        [key: 'domain_variability_file', enabled: params.posenrich_domains, required: ['gene','pfam_id','target_name','description','clan_acc','clan_name','ali_start','ali_end']],
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
                    error "${spec.key} '${input}' missing columns: ${missing}"
                }
                if (spec.key in ['domain_variability_file', 'ucr_positions_file'] && !columns.contains('coordinate_space')) {
                    error "${spec.key} '${input}' has no coordinate_space: convert the verified source with bin/build_reference_tables.py project-ppga or regenerate it before running"
                }
            }
            def manifest = file("${input}.manifest.json")
            if (manifest.exists()) {
                def meta = new groovy.json.JsonSlurper().parseText(manifest.text)
                if (meta.table_sha256 && hashFile.call(input) != meta.table_sha256)
                    error "Reference checksum mismatch: ${input}"
                if (meta.alignment_dir && params.alignment) {
                    def expected = file(meta.alignment_dir).toRealPath()
                    def productive = file(params.alignment).toRealPath()
                    if (expected != productive) error "${spec.key} belongs to alignment '${expected}', but this run uses '${productive}'"
                    if (meta.alignments && !checkedAlignments.contains(productive.toString())) {
                        log.info "Validating ${meta.alignments.size()} productive alignments for positional references: ${productive}"
                        def actualNames = java.nio.file.Files.list(productive).withCloseable { paths ->
                            paths.filter { java.nio.file.Files.isRegularFile(it) && it.fileName.toString().toLowerCase() ==~ /.*\.(fa|fasta|faa|fas)/ }
                                .map { it.fileName.toString() }.collect(java.util.stream.Collectors.toSet())
                        }
                        if (actualNames != meta.alignments.keySet()) error "Alignment inventory changed since reference generation: ${productive}"
                        meta.alignments.each { name, expectedHash ->
                            def alignment = productive.resolve(name)
                            if (!java.nio.file.Files.exists(alignment) || hashFile.call(alignment) != expectedHash)
                                error "Productive alignment changed since reference generation: ${alignment}"
                        }
                        checkedAlignments.add(productive.toString())
                    }
                }
                if (meta.map_dir && params.caas_map_dir && file(meta.map_dir).toRealPath() != file(params.caas_map_dir).toRealPath())
                    error "${spec.key} was validated with a different MAP directory: ${meta.map_dir}"
                if (meta.map_dir && meta.maps && params.caas_map_dir && !checkedMaps.contains(meta.map_dir)) {
                    log.info "Validating ${meta.maps.size()} productive MAP files for positional references"
                    def mapRoot = file(params.caas_map_dir)
                    def actualMaps = java.nio.file.Files.list(mapRoot).withCloseable { paths ->
                        paths.filter { java.nio.file.Files.isRegularFile(it) && it.fileName.toString().endsWith('.map.tsv') }
                            .map { it.fileName.toString() }.collect(java.util.stream.Collectors.toSet())
                    }
                    if (actualMaps != meta.maps.keySet()) error "MAP inventory changed since reference generation: ${mapRoot}"
                    meta.maps.each { name, expectedHash ->
                        def mapPath = file(params.caas_map_dir).resolve(name)
                        if (!java.nio.file.Files.exists(mapPath) || hashFile.call(mapPath) != expectedHash)
                            error "Productive MAP changed since reference generation: ${mapPath}"
                    }
                    checkedMaps.add(meta.map_dir)
                }
            }
        }
    }
}


workflow POSENRICH {
    take:
        run_accumulation         // Boolean: level 1 runs (CT post-processing or the standalone accumulation inputs exist)
        caas_channel             // filtered_discovery.tsv from CT_POSTPROC gene_filtering (or Channel.empty())
        background_channel       // cleaned_background_main.txt from CT_POSTPROC (or Channel.empty())
        tested_positions_channel // caastools background.output: TESTED positions, the starting point of the
                                 // universe U (required: Channel.empty() makes the run stop with an error)
        null_positions_channel   // null_positions.parquet from CAAS_CORE_MERGE_REDUCE, the filtered null SCORING also
                                 // uses (required, or the --caas_null_dir fallback)
        null_blocks_channel      // value channel: the list of block directories of the null (positions.parquet and
                                 // cycles.txt each), or Channel.empty() to use those of --caas_null_dir
        position_scores          // SCORING position_scores.tsv (Channel.empty() without --scoring)
        postproc                 // SCORING filtered_discovery.tsv: positions evaluated by CT
        null_universe            // SCORING null_universe_summary.tsv
        fade_sites_top_ch        // fade_sites_top.csv (or Channel.empty())
        fade_sites_bottom_ch     // fade_sites_bottom.csv (or Channel.empty())

    main:
        validatePositionReferences()

        // The universe U: one table for level 1 and level 2. Its observed positives come from filtered_discovery.tsv
        // (the one SCORING reads when it runs), so level 1 does not depend on SCORING.
        def accum_results = Channel.empty()
        def universe_out = null
        if (run_accumulation || params.scoring) {
            def observed_v = (params.scoring ? (postproc ?: Channel.empty()) : (caas_channel ?: Channel.empty()))
                .ifEmpty { params.accumulation_caas_input ? file(params.accumulation_caas_input, checkIfExists: true) : error("POSENRICH requires the observed CAAS (filtered_discovery.tsv): run CT post-processing or pass --accumulation_caas_input.") }
                .collect().map { it[0] }
            universe_out = POSENRICH_UNIVERSE(
                (background_channel ?: Channel.empty())
                    .ifEmpty { standaloneBackground() ?: error("POSENRICH requires the cleaned background (cleaned_background_main.txt): run CT post-processing or pass --accumulation_background_input.") }
                    .collect().map { it[0] },
                (tested_positions_channel ?: Channel.empty())
                    .ifEmpty { error "POSENRICH requires the tested positions (caastools background.output): run CAAS discovery or pass --posenrich_background_file." }
                    .collect().map { it[0] },
                (null_positions_channel ?: Channel.empty())
                    .ifEmpty { params.caas_null_dir ? file("${params.caas_null_dir}/null_positions.parquet") : error("POSENRICH requires null_positions.parquet: run CAAS permulation or pass --caas_null_dir.") }
                    .collect().map { it[0] },
                observed_v
            )
        }

        // Level 1: the exact test per gene and direction on the universe table.
        if (run_accumulation) {
            log.info "POSENRICH level 1: Accumulation"
            accum_results = ACCUMULATION_TEST(Channel.of('top', 'bottom', 'all'), universe_out.universe).results
        }

        // The registry and its components read SCORING's tables.
        def position_char_ch = Channel.empty()
        def report_ch = Channel.empty()
        if (params.scoring) {
            // Supplied tables take priority. Enabled layers with blank files are
            // generated once on the productive alignment columns.
            def domain_annotations_ch = Channel.value(file('NO_ANNOT_DOMAIN'))
            if (params.posenrich_domains) {
                if (params.domain_variability_file) {
                    domain_annotations_ch = Channel.value(precomputedOrSentinel('domain_variability_file', params.domain_variability_file, 'NO_ANNOT_DOMAIN'))
                } else {
                    DOMAIN_VARIABILITY_GENERATION(Channel.value(file(params.alignment)))
                    domain_annotations_ch = DOMAIN_VARIABILITY_GENERATION.out.domain_variability_file
                }
            }
            def ucr_annotations_ch = Channel.value(file('NO_ANNOT_UCR'))
            if (params.posenrich_ucr) {
                if (params.ucr_positions_file) {
                    ucr_annotations_ch = Channel.value(precomputedOrSentinel('ucr_positions_file', params.ucr_positions_file, 'NO_ANNOT_UCR'))
                } else {
                    UCR_GENERATION(Channel.value(file(params.alignment)))
                    ucr_annotations_ch = UCR_GENERATION.out.ucr_positions_file
                }
            }
            def fade_top_ch = (fade_sites_top_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_TOP') }
            def fade_bot_ch = (fade_sites_bottom_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_BOTTOM') }

            def registry_out = POSITION_REGISTRY(
                position_scores,
                (postproc ?: Channel.empty()).ifEmpty { file('NO_POSTPROC') },
                (null_universe ?: Channel.empty()).ifEmpty { file('NO_NULL_UNIVERSE') },
                Channel.value(precomputedOrSentinel('caas_map_dir', params.caas_map_dir, 'NO_ANNOT_MAP')),
                (background_channel ?: Channel.empty()).ifEmpty { file('NO_BACKGROUND') }.collect().map { it[0] },
                Channel.value(precomputedOrSentinel('gene_ensembl_file', params.gene_ensembl_file, 'NO_ANNOT_ENSEMBL')),
                domain_annotations_ch,
                ucr_annotations_ch,
                Channel.value(precomputedOrSentinel('fubar_sites_file', params.fubar_sites_file, 'NO_ANNOT_FUBAR')),
                fade_top_ch, fade_bot_ch,
                Channel.value(precomputedOrSentinel('posenrich_custom_file', params.posenrich_custom_file, 'NO_ANNOT_CUSTOM'))
            )
            def characterization_out = POSENRICH_CHARACTERIZE(registry_out.registry)
            def proximity_out        = POSENRICH_PROXIMITY(registry_out.registry)

            // Level 2 needs the cycles of the null.
            // The block directories of a prior null need type: 'dir' (files() matches regular files by default).
            def prior_blocks = params.caas_null_dir ? files("${params.caas_null_dir}/blocks/*", type: 'dir').sort() : []
            def blocks_v = (null_blocks_channel ?: Channel.empty())
                .ifEmpty { prior_blocks ?: error("POSENRICH level 2 requires the block directories of the null: run CAAS permulation or pass --caas_null_dir with a blocks/ directory.") }
            def within_gene_out = POSENRICH_WITHIN_GENE(universe_out.universe, universe_out.manifest,
                                                        position_scores, registry_out.registry, blocks_v)

            // The report closes the module. Level 1 is a sentinel when it does not run (the report states it).
            def accumulation_v = accum_results.collect().ifEmpty { [file('NO_ACCUMULATION')] }
            def report_out = POSENRICH_REPORT(universe_out.universe_dir, registry_out.registry,
                                              characterization_out.characterization, accumulation_v,
                                              within_gene_out.within_gene, proximity_out.proximity)
            position_char_ch = characterization_out.position_char
            report_ch        = report_out.report
        }

    emit:
        accumulation  = accum_results
        position_char = position_char_ch
        report        = report_ch
}
