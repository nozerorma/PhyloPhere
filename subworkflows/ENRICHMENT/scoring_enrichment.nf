#!/usr/bin/env nextflow

/*
 * SCORING_STRING / SCORING_COMPARE
 * ────────────────────────────────────────────
 * STRING enrichment runs on the gated directional 9-slice gene lists produced
 * by SCORING_COMPUTE. Ranked enrichment is handled by SCORING_FCS_REPORT
 * (subworkflows/ENRICHMENT/fcs.nf).
 *
 *   SCORING_STRING_REPORT  → string_results/**, HTML (network/set question)
 *   SCORING_COMPARE_REPORT → compare_results/**, HTML (top vs bottom: FCS + STRING)
 */


// ─────────────────────────────────────────────────────────────────────────────
// SCORING_STRING_REPORT
// Runs STRING_scoring.Rmd (STRINGdb walktrap clustering + enrichment) on the
// 9 ranked slices.
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_STRING_REPORT {
    tag "scoring_string|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries    3

    publishDir path: "${params.outdir}/string",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/string/string_results",
               mode: 'copy', overwrite: true,
               pattern: 'string_results/**'
    publishDir path: "${params.outdir}/string/string_summary",
               mode: 'copy', overwrite: true,
               pattern: 'string_summary/**'
    publishDir path: "${params.outdir}/string/string_plots",
               mode: 'copy', overwrite: true,
               pattern: 'string_plots/**'
    publishDir path: "${params.outdir}/string/string_networks",
               mode: 'copy', overwrite: true,
               pattern: 'string_networks/**'

    input:
    path gene_lists
    path background
    path gene_scores

    output:
    path "STRING_general_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "string_results/**",                  emit: string_results,          optional: true
    path "string_summary/**",                  emit: string_summary,          optional: true
    path "string_plots/**",                    emit: string_plots,            optional: true
    path "string_networks/**",                 emit: string_networks,         optional: true
    path "string_results/*_enrichment.tsv",    emit: string_enrichment_tsvs,  optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def species   = params.string_species           ?: 9606
    def req_score = params.string_required_score    ?: 400
    def net_score = params.string_network_score_thr ?: 700
    def fdr_thr   = params.string_fdr              ?: 0.1
    def top_thr   = params.string_top_thr          ?: 15
    def bg_name   = background.getName().replace("'", "\\'")
    def gs_arg    = (gene_scores.name =~ /^NO_GENE_SCORES/) ? 'NULL' : "'${gene_scores}'"

    def stage_cmd = """
        cp -R ${local_dir}/* .

        # Convert slice_*.tsv to *.txt in the current directory
        for f in ${gene_lists}/slice_*.tsv; do
            if [ -f "\$f" ]; then
                basename=\$(basename "\$f" .tsv)
                name=\${basename#slice_}
                # Extract first column (Gene) except header (1st line)
                tail -n +2 "\$f" | cut -f1 | grep -v "^[[:space:]]*\$" > "\${name}.txt"
            fi
        done
    """

    def render_cmd = """
        Rscript -e "
            rmarkdown::render(
                'STRING_general.Rmd',
                params = list(
                    background_file     = '${background}',
                    background_basename = '${bg_name}',
                    output_dir          = 'string_results',
                    project_name        = '${traitname}',
                    species             = ${species},
                    required_score      = ${req_score},
                    network_score_thr   = ${net_score},
                    fdr_thr             = ${fdr_thr},
                    top_thr             = ${top_thr},
                    gene_scores_file    = ${gs_arg}
                ),
                output_file = 'STRING_general_${traitname}.html'
            )
        "
    """

    def compat_cmd = """
        mkdir -p string_results
        if [ -d string_summary ]; then
            for f in string_summary/*_results.tsv; do
                if [ -f "\$f" ]; then
                    basename=\$(basename "\$f" _results.tsv)
                    cp "\$f" "string_results/\${basename}_enrichment.tsv"
                fi
            done
        fi
    """

    if (params.use_singularity || params.use_apptainer) {
        """
        ${stage_cmd}
        /usr/local/bin/_entrypoint.sh ${render_cmd}
        ${compat_cmd}
        """
    } else {
        """
        ${stage_cmd}
        ${render_cmd}
        ${compat_cmd}
        """
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SCORING_COMPARE_REPORT
// Comparative top-vs-bottom analysis across FCS and STRING outputs.
//
// Inputs (all staged flat in work dir by Nextflow):
//   fcs_all_results : fcs_results/fcs_all_results.tsv  (single file)
//   string_tsvs     : string_results/*_enrichment.tsv  (collection, optional)
//
// The script sorts staged files into cmp_fcs/ and cmp_string/ before
// invoking COMPARE_scoring.Rmd.
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_COMPARE_REPORT {
    tag "scoring_compare|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/compare",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/compare/compare_results",
               mode: 'copy', overwrite: true,
               pattern: 'compare_results/**'

    input:
    path fcs_all_results   // fcs_results/fcs_all_results.tsv or NO_FCS_ALL placeholder
    path string_tsvs       // string_results/*_enrichment.tsv files (may be placeholder)

    output:
    path "COMPARE_scoring_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "compare_results/**", emit: compare_results, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def fdr_thr   = params.scoring_compare_fdr   ?: params.fcs_fdr ?: 0.1
    def top_n     = params.scoring_compare_top_n  ?: 20

    def render_cmd = """
        Rscript -e "
            rmarkdown::render(
                'COMPARE_scoring.Rmd',
                params = list(
                    fcs_dir    = 'cmp_fcs',
                    string_dir = 'cmp_string',
                    fdr_thr    = ${fdr_thr},
                    top_n      = ${top_n},
                    traitname  = '${traitname}'
                ),
                output_file = 'COMPARE_scoring_${traitname}.html'
            )
        "
    """

    def stage_cmd = """
        cp -R ${local_dir}/* .

        mkdir -p cmp_fcs cmp_string

        # FCS all-results table
        [ -f "fcs_all_results.tsv" ] && cp "fcs_all_results.tsv" cmp_fcs/ || true

        # STRING enrichment TSVs
        for f in *_enrichment.tsv; do
            [ -f "\$f" ] && cp "\$f" cmp_string/ || true
        done
    """

    if (params.use_singularity || params.use_apptainer) {
        """
        ${stage_cmd}
        /usr/local/bin/_entrypoint.sh ${render_cmd}
        """
    } else {
        """
        ${stage_cmd}
        ${render_cmd}
        """
    }
}
