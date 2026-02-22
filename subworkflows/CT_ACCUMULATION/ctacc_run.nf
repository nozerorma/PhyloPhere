#!/usr/bin/env nextflow

/*
#  CT Accumulation: Aggregate and Randomize processes
*/

process CT_ACCUMULATION_AGGREGATE {
    tag "ct_accumulation_aggregate"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/accumulation/aggregation", mode: 'copy', overwrite: true,
               pattern: '*.csv'

    input:
    val  alignment_dir
    path genomic_info
    path species_list
    path metadata_caas
    path bg_caas

    output:
    path "*_global.csv",    emit: global_csv
    path "*_deciles.csv",   emit: deciles, optional: true

    script:
    def local_dir    = "${baseDir}/subworkflows/CT_ACCUMULATION/local"
    def ali_fmt      = params.ali_format ?: 'phylip-relaxed'
    def out_pfx      = 'accumulation'
    def log_level    = params.accumulation_log_level ?: 'INFO'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh python main.py \\
            --tool aggregate \\
            --alignment-dir '${alignment_dir}' \\
            --alignment-format '${ali_fmt}' \\
            --genomic-info '${genomic_info}' \\
            --species-list '${species_list}' \\
            --metadata-caas '${metadata_caas}' \\
            --bg-caas '${bg_caas}' \\
            --output-prefix '${out_pfx}' \\
            --log-level '${log_level}'
        """
    } else {
        """
        cp -R ${local_dir}/* .

        python main.py \\
            --tool aggregate \\
            --alignment-dir '${alignment_dir}' \\
            --alignment-format '${ali_fmt}' \\
            --genomic-info '${genomic_info}' \\
            --species-list '${species_list}' \\
            --metadata-caas '${metadata_caas}' \\
            --bg-caas '${bg_caas}' \\
            --output-prefix '${out_pfx}' \\
            --log-level '${log_level}'
        """
    }
}

process CT_ACCUMULATION_RANDOMIZE {
    tag "ct_accumulation_randomize"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/accumulation/randomization", mode: 'copy', overwrite: true,
               pattern: '{*_aggregated_results.csv,gene_lists/**}'

    input:
    path global_csv
    path caas_csv

    output:
    path "*_aggregated_results.csv", emit: results
    path "gene_lists/**",            emit: gene_lists, optional: true

    script:
    def local_dir    = "${baseDir}/subworkflows/CT_ACCUMULATION/local"
    def out_pfx      = 'accumulation'
    def rand_type    = params.accumulation_randomization_type ?: 'naive'
    def n_rands      = params.accumulation_n_randomizations   ?: 10000
    def fdr_thr      = params.accumulation_fdr                ?: 0.05
    def log_level    = params.accumulation_log_level          ?: 'INFO'
    def workers_flag = params.accumulation_workers ? "--workers ${params.accumulation_workers}" : ''
    def seed_flag    = params.accumulation_seed    ? "--global-seed ${params.accumulation_seed}"   : ''

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh python main.py \\
            --tool randomize \\
            --global-csv '${global_csv}' \\
            --caas-csv '${caas_csv}' \\
            --output-prefix '${out_pfx}' \\
            --randomization-type '${rand_type}' \\
            --n-randomizations ${n_rands} \\
            --fdr-threshold ${fdr_thr} \\
            ${workers_flag} ${seed_flag} \\
            --log-level '${log_level}'
        """
    } else {
        """
        cp -R ${local_dir}/* .

        python main.py \\
            --tool randomize \\
            --global-csv '${global_csv}' \\
            --caas-csv '${caas_csv}' \\
            --output-prefix '${out_pfx}' \\
            --randomization-type '${rand_type}' \\
            --n-randomizations ${n_rands} \\
            --fdr-threshold ${fdr_thr} \\
            ${workers_flag} ${seed_flag} \\
            --log-level '${log_level}'
        """
    }
}
