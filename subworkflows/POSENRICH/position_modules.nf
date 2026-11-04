#!/usr/bin/env nextflow
// POSENRICH position modules: universe U, common position registry, characterization,
// within-gene concentration (level 2), proximity and the report.
// The workflow that wires them is workflows/posenrich.nf.

def modularSourceHash(List names) {
    names = (names + ['src/posenrich_contract.py', '../../SCORING/local/src/null_position_stats.py', '../../SCORING/local/src/scoring_input.py', '../../CT_DISAMBIGUATION/local/src/core/scores.py', '../../VEP/local/src/vep_common.py']).unique()
    def sources = names.collect { file("${baseDir}/subworkflows/POSENRICH/local/${it}").text }.join('\n') +
        '\nintegrity=' + (System.getenv('PHYLOPHERE_VERIFY_CHECKSUMS') ?: '0')
    return java.security.MessageDigest.getInstance('SHA-256')
        .digest(sources.getBytes('UTF-8')).encodeHex().toString()
}

// Command-line fragment for an optional source; NO_* sentinels and dangling paths mean absent.
// A function rather than a closure so the linter resolves it.
def optionalFlag(flag, source) {
    return source.name.startsWith('NO_') ? '' : "--${flag} '${source}'"
}

process POSITION_REGISTRY {
    cache true
    label 'process_registry'
    memory { params.position_registry_memory ?: (background.size() < 20000 ? '8 GB' : '64 GB') }
    publishDir path: { universe.name.startsWith('NO_') ? "${params.outdir}/position_characterization" : "${params.outdir}/posenrich" }, mode: 'copy', overwrite: true, pattern: 'registry'

    input:
    path position_scores
    path postproc
    path null_universe
    path map_dir
    path background
    path gene_ensembl
    path domain_file
    path ucr_file
    path fubar_file
    path fade_top
    path fade_bottom
    path custom_file
    path universe
    path universe_manifest
    path annotation_manifests
    path position_calibration

    output:
    path 'registry', emit: registry

    script:
    def py3 = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def source_hash = modularSourceHash(['src/posenrich_registry.py', 'src/position_sources.py', 'src/posenrich_universe.py', 'src/posenrich_members.py'])
    // The registry tolerates every absent source (see optionalFlag).
    """
    # POSENRICH sources SHA256: ${source_hash}
    ${py3} ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_registry.py \\
        ${optionalFlag('position_scores', position_scores)} \\
        ${optionalFlag('postproc', postproc)} \\
        ${optionalFlag('null_universe', null_universe)} \\
        ${optionalFlag('map_dir', map_dir)} \\
        ${optionalFlag('cleaned_background', background)} \\
        ${optionalFlag('gene_ensembl_file', gene_ensembl)} \\
        ${optionalFlag('domain_variability_file', domain_file)} \\
        ${optionalFlag('ucr_positions_file', ucr_file)} \\
        ${optionalFlag('fubar_sites_file', fubar_file)} \\
        ${optionalFlag('fade_sites_top_file', fade_top)} \\
        ${optionalFlag('fade_sites_bottom_file', fade_bottom)} \\
        ${optionalFlag('custom_marker_file', custom_file)} \\
        --window_size_bp ${params.scoring_window_size_bp ?: 1000000} \\
        ${universe.name.startsWith('NO_') ? '' : "--universe '${universe}' --universe-manifest '${universe_manifest}' --position-calibration '${position_calibration}'"} \\
        --output_dir registry
    """
}

process POSENRICH_CHARACTERIZE {
    cache true
    label 'process_low'
    memory { params.position_module_memory ?: (registry.toRealPath().resolve('positions.tsv').size() < 200000000 ? '4 GB' : '32 GB') }
    publishDir path: { task.process.contains('POSITION_CHARACTERIZATION') ? "${params.outdir}/position_characterization" : "${params.outdir}/posenrich" }, mode: 'copy', overwrite: true, pattern: 'characterization'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'characterization', emit: characterization
    path 'characterization/position_characterization.tsv', emit: position_char

    script:
    def py3 = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def source_hash = modularSourceHash(['src/posenrich_members.py', 'src/posenrich_characterize.py'])
    """
    # POSENRICH sources SHA256: ${source_hash}
    ${py3} ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_characterize.py \\
        --registry_dir registry --output_dir characterization --fdr_levels ${params.evidence_fdr_levels ?: '0.05,0.1,0.15'}
    """
}

process POSENRICH_UNIVERSE {
    cache true
    tag "posenrich_universe"
    label 'process_medium'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'universe'

    input:
    path background_genes       // cleaned_background_main.txt (the gene list of U)
    path background_positions   // caastools background.output (tested positions) or a NO_* sentinel
    path null_positions         // null_positions.parquet (the filtered null SCORING uses)
    path observed               // filtered_discovery.tsv: observed positives
    path observed_contract
    path roster
    path counts
    path blocks, stageAs: 'blk_*'

    output:
    path 'universe',                            emit: universe_dir
    path 'universe/universe.tsv.gz',            emit: universe
    path 'universe/universe_manifest.json',     emit: manifest
    path 'universe/position_calibration.tsv', emit: calibration

    script:
    def py3 = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def rscript = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = modularSourceHash(['src/posenrich_universe.py', '../../SCORING/local/src/validate_position_matrices.R'])
    def tested_arg = background_positions.name.startsWith('NO_') ? '' : "--background-positions '${background_positions}'"
    """
    # POSENRICH sources SHA256: ${source_hash}
    ${py3} ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_universe.py \\
        --background-genes '${background_genes}' ${tested_arg} \\
        --null-positions '${null_positions}' \\
        --observed '${observed}' \\
        --observed-contract '${observed_contract}' --roster '${roster}' --counts '${counts}' --blocks blk_* \\
        --output-dir universe
    ${rscript} ${baseDir}/subworkflows/SCORING/local/src/validate_position_matrices.R universe blk_*
    """
}

// Level 2: concentration of the signal inside genes, calibrated with the cycles of the null. The null blocks are read
// one at a time, so the memory follows the size of a block (1.5 GB for 500 cycles and 474,450 term-gene pairs).
process POSENRICH_WITHIN_GENE {
    cache true
    tag "posenrich_within_gene"
    label 'process_medium'
    memory { params.posenrich_within_gene_memory ?: '16 GB' }
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'within_gene'

    input:
    path universe, stageAs: 'universe.tsv.gz'
    path universe_manifest, stageAs: 'universe_manifest.json'
    path observed                                   // position_scores.tsv
    path registry, stageAs: 'registry'
    path blocks, stageAs: 'blk_*'                   // block directories of the null (positions.parquet, cycles.txt)

    output:
    path 'within_gene', emit: within_gene

    script:
    def py3 = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def source_hash = modularSourceHash(['src/posenrich_within_gene.py', 'src/posenrich_members.py', 'src/posenrich_universe.py'])
    def collections_arg = params.posenrich_collections ? "--collections '${params.posenrich_collections}'" : ''
    """
    # POSENRICH sources SHA256: ${source_hash}
    ${py3} ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_within_gene.py \\
        --universe universe.tsv.gz --universe-manifest universe_manifest.json --observed '${observed}' --registry-dir registry \\
        --blocks blk_* --output-dir within_gene \\
        --min-genes ${params.posenrich_min_genes ?: 2} \\
        --holdout ${params.posenrich_holdout ?: 0} ${collections_arg}
    """
}

process POSENRICH_PROXIMITY {
    cache true
    label 'process_low'
    memory { params.position_module_memory ?: (registry.toRealPath().resolve('positions.tsv').size() < 200000000 ? '4 GB' : '32 GB') }
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'proximity'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'proximity', emit: proximity

    script:
    def py3 = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def source_hash = modularSourceHash(['src/posenrich_members.py', 'src/posenrich_proximity.py'])
    """
    # POSENRICH sources SHA256: ${source_hash}
    ${py3} ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_proximity.py \\
        --registry_dir registry --output_dir proximity \\
        --radios '${params.posenrich_radios ?: '0,5,10'}'
    """
}

// The report closes the module. Absent components arrive as NO_* sentinels or empty directories and are stated as unavailable.
process POSENRICH_REPORT {
    cache true
    tag "posenrich_report|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: '{*.html,report_settings.json}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path universe,         stageAs: 'universe'
    path registry,         stageAs: 'registry'
    path characterization, stageAs: 'characterization'
    path accumulation,     stageAs: 'accumulation/*'    // accumulation_<direction>_results.csv, or the NO_ACCUMULATION sentinel
    path within_gene,      stageAs: 'within_gene'
    path proximity,        stageAs: 'proximity'

    output:
    path "13.POSENRICH_report_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "report_settings.json", emit: settings

    script:
    def rscript = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = modularSourceHash(['13.POSENRICH_report.Rmd']) + modularSourceHash(['../../ENRICHMENT/local/src/report_presentation.R', '../../TRAIT_ANALYSIS/local/src/report_style.R'])
    def local_dir = "${baseDir}/subworkflows/POSENRICH/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def fdr_levels = params.evidence_fdr_levels ?: '0.05,0.1,0.15'
    """
    # POSENRICH sources SHA256: ${source_hash}
    cp -R ${local_dir}/* .
    mkdir -p src
    cp ${baseDir}/subworkflows/ENRICHMENT/local/src/report_presentation.R ${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R src/
    ${rscript} -e "
        rmarkdown::render(
            '13.POSENRICH_report.Rmd',
            params = list(
                universe_dir         = 'universe',
                registry_dir         = 'registry',
                characterization_dir = 'characterization',
                accumulation_dir     = 'accumulation',
                within_gene_dir      = 'within_gene',
                proximity_dir        = 'proximity',
                traitname            = '${traitname}',
                fdr_levels           = '${fdr_levels}'
            ),
            output_file = '13.POSENRICH_report_${traitname}.html'
        )
        jsonlite::write_json(list(schema='phylophere.posenrich_report.v1',
                                 traitname='${traitname}',fdr_levels=as.numeric(strsplit('${fdr_levels}',',')[[1]])),
                             'report_settings.json',auto_unbox=TRUE,pretty=TRUE,digits=NA)
    "
    """
}
