#!/usr/bin/env nextflow

/*
 * POSENRICH — Position-wise enrichment (Fisher-exact, fixed cutoffs)
 * ────────────────────────────────────────────────────────────────────
 * NOT the gene-level FCS. Builds position-level GMTs (Pfam, Bins, Orthogroups,
 * COSMIC, UCR core/flank, positive/purifying selection) and Fisher-tests all of
 * them, plus the broad functional characterization layers, at fixed top-10/5/1%
 * cutoffs, per direction (global/top/bottom). The background is large
 * (~1.47M tested positions), large enough that Fisher's exact test has power
 * to flag biologically negligible deviations as "significant" on p-value
 * alone, so significance is dual-gated on p_adj AND a minimum fold-enrichment
 * magnitude (posenrich_padj_thr / posenrich_fold_thr).
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

process POSENRICH_BUILD_GMT {
    label 'process_low'
    publishDir "${params.outdir}/posenrich/gmts", mode: 'copy', overwrite: true

    input:
    path gene_ensembl_file
    path domain_variability_file
    path ucr_positions_file
    path fubar_sites_file
    path egg_members_file
    path egg_annotations_file
    path map_dir
    path cosmic_db
    path cleaned_background

    output:
    path "*.gmt", emit: gmts
    path "characterization_layers.tsv", emit: charset

    script:
    def cosmic_arg = cosmic_db.name != 'NO_FILE' ? "--cosmic_db ${cosmic_db}" : ""
    def bg_arg     = cleaned_background.name != 'NO_FILE' ? "--cleaned_background ${cleaned_background}" : ""
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/build_position_gmt.py \
        --gene_ensembl_file ${gene_ensembl_file} \
        --domain_variability_file ${domain_variability_file} \
        --ucr_positions_file ${ucr_positions_file} \
        --fubar_sites_file ${fubar_sites_file} \
        --egg_members_file ${egg_members_file} \
        --egg_annotations_file ${egg_annotations_file} \
        --map_dir ${map_dir} \
        ${cosmic_arg} \
        ${bg_arg} \
        --output_dir .
    """
}

process POSENRICH_RUN {
    label 'process_medium'
    publishDir "${params.outdir}/posenrich", mode: 'copy', overwrite: true

    input:
    path caas_file
    path "gmts/*"
    path characterization_layers
    path universe
    path background_output
    path annot_file
    val min_size
    val max_size

    output:
    path "posenrich_characterization.tsv", emit: results

    script:
    // No permulation here: the exact Fisher p already corrects for the fixed
    // cutoff (no cutoff search to correct for), and the CAAS permulation null
    // (naive random relabeling, random cross-tree pairs) does not isolate the
    // phenotype at position level. Significance is instead dual-gated on
    // p_adj AND fold-enrichment (posenrich_padj_thr / posenrich_fold_thr) —
    // see posenrich_enrich.py's module docstring for why the fold gate exists.
    def annot_arg = annot_file.name != 'NO_FILE' ? "--annot-file ${annot_file}" : ""
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_enrich.py \
        --obs-scores ${caas_file} \
        --gmt-dir gmts \
        --characterization ${characterization_layers} \
        --universe ${universe} \
        --background ${background_output} \
        ${annot_arg} \
        --min-size ${min_size} \
        --max-size ${max_size} \
        --padj-thr ${params.posenrich_padj_thr} \
        --fold-thr ${params.posenrich_fold_thr} \
        --output-dir .
    """
}

process POSENRICH_REPORT {
    tag "posenrich_report|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/posenrich",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    path results

    output:
    path "16.Position_enrichment_report_${params.traitname ?: 'unknown_trait'}.html", emit: report

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '16.Position_enrichment_report.Rmd',
                params = list(
                    results_file = '${results}',
                    traitname = '${traitname}',
                    padj_thr = ${params.posenrich_padj_thr},
                    fold_thr = ${params.posenrich_fold_thr}
                ),
                output_file = '16.Position_enrichment_report_${traitname}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "
            rmarkdown::render(
                '16.Position_enrichment_report.Rmd',
                params = list(
                    results_file = '${results}',
                    traitname = '${traitname}',
                    padj_thr = ${params.posenrich_padj_thr},
                    fold_thr = ${params.posenrich_fold_thr}
                ),
                output_file = '16.Position_enrichment_report_${traitname}.html'
            )
        "
        """
    }
}

workflow POSENRICH {
    take:
    gene_ensembl_file
    domain_variability_file
    ucr_positions_file
    fubar_sites_file
    egg_members_file
    egg_annotations_file
    map_dir
    cosmic_db
    cleaned_background
    caas_file
    background_output
    annot_file
    min_size
    max_size

    main:
    POSENRICH_BUILD_GMT(
        gene_ensembl_file,
        domain_variability_file,
        ucr_positions_file,
        fubar_sites_file,
        egg_members_file,
        egg_annotations_file,
        map_dir,
        cosmic_db,
        cleaned_background
    )

    POSENRICH_RUN(
        caas_file,
        POSENRICH_BUILD_GMT.out.gmts,
        POSENRICH_BUILD_GMT.out.charset,
        cleaned_background,
        background_output,
        annot_file,
        min_size,
        max_size
    )

    POSENRICH_REPORT(
        POSENRICH_RUN.out.results
    )

    emit:
    results = POSENRICH_RUN.out.results
    report  = POSENRICH_REPORT.out.report
}
