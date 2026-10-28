#!/usr/bin/env nextflow
// POSENRICH position modules: universe U, common position registry, characterization,
// within-gene concentration (level 2), proximity and the report.
// The workflow that wires them is workflows/posenrich.nf.

def modularSourceHash(List names) {
    def sources = names.collect { file("${baseDir}/subworkflows/POSENRICH/local/${it}").text }.join('\n')
    return java.security.MessageDigest.getInstance('SHA-256')
        .digest(sources.getBytes('UTF-8')).encodeHex().toString()
}

process POSITION_REGISTRY {
    label 'process_registry'
    memory { params.position_registry_memory ?: (background.size() < 20000 ? '8 GB' : '64 GB') }
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'registry'

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

    output:
    path 'registry', emit: registry

    script:
    def source_hash = modularSourceHash(['src/posenrich_registry.py', 'src/position_sources.py', 'src/posenrich_universe.py', 'src/posenrich_members.py'])
    // NO_* sentinels and dangling paths mean absent; the registry tolerates every absent source.
    def optional = { flag, source -> source.name.startsWith('NO_') ? '' : "--${flag} '${source}'" }
    """
    # POSENRICH sources SHA256: ${source_hash}
    python3 ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_registry.py \\
        ${optional('position_scores', position_scores)} \\
        ${optional('postproc', postproc)} \\
        ${optional('null_universe', null_universe)} \\
        ${optional('map_dir', map_dir)} \\
        ${optional('cleaned_background', background)} \\
        ${optional('gene_ensembl_file', gene_ensembl)} \\
        ${optional('domain_variability_file', domain_file)} \\
        ${optional('ucr_positions_file', ucr_file)} \\
        ${optional('fubar_sites_file', fubar_file)} \\
        ${optional('fade_sites_top_file', fade_top)} \\
        ${optional('fade_sites_bottom_file', fade_bottom)} \\
        ${optional('custom_marker_file', custom_file)} \\
        --window_size_bp ${params.scoring_window_size_bp ?: 1000000} \\
        --output_dir registry
    """
}

process POSENRICH_CHARACTERIZE {
    label 'process_low'
    memory { params.position_module_memory ?: (registry.toRealPath().resolve('positions.tsv').size() < 200000000 ? '4 GB' : '32 GB') }
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'characterization'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'characterization', emit: characterization
    path 'characterization/position_characterization.tsv', emit: position_char

    script:
    def source_hash = modularSourceHash(['src/posenrich_members.py', 'src/posenrich_characterize.py'])
    """
    # POSENRICH sources SHA256: ${source_hash}
    python3 ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_characterize.py \\
        --registry_dir registry --output_dir characterization
    """
}

process POSENRICH_UNIVERSE {
    tag "posenrich_universe"
    label 'process_medium'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'universe'

    input:
    path background_genes       // cleaned_background_main.txt (the gene list of U)
    path background_positions   // caastools background.output (tested positions) or a NO_* sentinel
    path null_positions         // null_positions.parquet (the filtered null SCORING uses)
    path observed               // filtered_discovery.tsv (or position_scores.tsv): the observed positives

    output:
    path 'universe',                            emit: universe_dir
    path 'universe/universe.tsv.gz',            emit: universe
    path 'universe/universe_manifest.json',     emit: manifest

    script:
    def source_hash = modularSourceHash(['src/posenrich_universe.py'])
    def tested_arg = background_positions.name.startsWith('NO_') ? '' : "--background-positions '${background_positions}'"
    """
    # POSENRICH sources SHA256: ${source_hash}
    python3 ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_universe.py \\
        --background-genes '${background_genes}' ${tested_arg} \\
        --null-positions '${null_positions}' \\
        --observed '${observed}' \\
        --output-dir universe
    """
}

// Level 2: concentration of the signal inside genes, calibrated with the cycles of the null. The null blocks are read
// one at a time, so the memory follows the size of a block (1.5 GB for 500 cycles and 474,450 term-gene pairs).
process POSENRICH_WITHIN_GENE {
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
    def source_hash = modularSourceHash(['src/posenrich_within_gene.py', 'src/posenrich_members.py', 'src/posenrich_universe.py'])
    def collections_arg = params.posenrich_collections ? "--collections '${params.posenrich_collections}'" : ''
    """
    # POSENRICH sources SHA256: ${source_hash}
    python3 ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_within_gene.py \\
        --universe universe.tsv.gz --observed '${observed}' --registry-dir registry \\
        --blocks blk_* --output-dir within_gene \\
        --min-nonzero-cycles ${params.posenrich_min_nonzero_cycles} \\
        --holdout ${params.posenrich_holdout} ${collections_arg}
    """
}

process POSENRICH_PROXIMITY {
    label 'process_low'
    memory { params.position_module_memory ?: (registry.toRealPath().resolve('positions.tsv').size() < 200000000 ? '4 GB' : '32 GB') }
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'proximity'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'proximity', emit: proximity

    script:
    def source_hash = modularSourceHash(['src/posenrich_members.py', 'src/posenrich_proximity.py'])
    """
    # POSENRICH sources SHA256: ${source_hash}
    python3 ${baseDir}/subworkflows/POSENRICH/local/src/posenrich_proximity.py \\
        --registry_dir registry --output_dir proximity \\
        --radios '${params.posenrich_radios ?: '0,5,10'}'
    """
}

// The report closes the module. Absent components arrive as NO_* sentinels or empty directories and are stated as unavailable.
process POSENRICH_REPORT {
    tag "posenrich_report|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path universe,         stageAs: 'universe'
    path registry,         stageAs: 'registry'
    path characterization, stageAs: 'characterization'
    path accumulation,     stageAs: 'accumulation/*'    // accumulation_<direction>_results.csv, or the NO_ACCUMULATION sentinel
    path within_gene,      stageAs: 'within_gene'
    path proximity,        stageAs: 'proximity'

    output:
    path "14.POSENRICH_report_${params.traitname ?: 'unknown_trait'}.html", emit: report

    script:
    def source_hash = modularSourceHash(['14.POSENRICH_report.Rmd']) + modularSourceHash(['../../ENRICHMENT/local/src/report_presentation.R'])
    def local_dir = "${baseDir}/subworkflows/POSENRICH/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def fdr = params.accumulation_fdr ?: 0.1
    """
    # POSENRICH sources SHA256: ${source_hash}
    cp -R ${local_dir}/* .
    mkdir -p src
    cp ${baseDir}/subworkflows/ENRICHMENT/local/src/report_presentation.R src/
    Rscript -e "
        rmarkdown::render(
            '14.POSENRICH_report.Rmd',
            params = list(
                universe_dir         = 'universe',
                registry_dir         = 'registry',
                characterization_dir = 'characterization',
                accumulation_dir     = 'accumulation',
                within_gene_dir      = 'within_gene',
                proximity_dir        = 'proximity',
                traitname            = '${traitname}',
                fdr_threshold        = ${fdr}
            ),
            output_file = '14.POSENRICH_report_${traitname}.html'
        )
    "
    """
}
