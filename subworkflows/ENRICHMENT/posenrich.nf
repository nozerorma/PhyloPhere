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
    path pai3d_db
    path cleaned_background
    path fade_sites_top_file
    path fade_sites_bottom_file

    output:
    path "*.gmt", emit: gmts
    path "characterization_layers.tsv", emit: charset
    path "cosmic_coverage_genes.txt", optional: true, emit: cosmic_coverage
    path "pai3d_coverage_genes.txt", optional: true, emit: pai3d_coverage

    script:
    // Each optional input's absent-value sentinel is a uniquely-named 'NO_FILE_*'
    // path (see workflows/enrichment.nf) rather than a shared literal 'NO_FILE' —
    // staging two path inputs under the identical filename in one task directory
    // is a Nextflow "input file name collision", which is exactly what happens
    // if two or more of these optional inputs are absent in the same run.
    def cosmic_arg = !(cosmic_db.name =~ /^NO_FILE/) ? "--cosmic_db ${cosmic_db}" : ""
    def pai3d_arg  = !(pai3d_db.name =~ /^NO_FILE/) ? "--pai3d_db ${pai3d_db}" : ""
    def bg_arg     = !(cleaned_background.name =~ /^NO_FILE/) ? "--cleaned_background ${cleaned_background}" : ""
    // FADE_top_sig/FADE_bottom_sig position group (§ FADE_JSON_TO_CSV) — the
    // classic BF>=100 sites, joined directly on Gene:Position (same
    // coordinate space CAAS's own Position column uses, no map_cache lookup).
    def fade_top_arg    = !(fade_sites_top_file.name    =~ /^NO_FILE/) ? "--fade_sites_top_file ${fade_sites_top_file}"       : ""
    def fade_bottom_arg = !(fade_sites_bottom_file.name =~ /^NO_FILE/) ? "--fade_sites_bottom_file ${fade_sites_bottom_file}" : ""
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
        ${pai3d_arg} \
        ${bg_arg} \
        ${fade_top_arg} \
        ${fade_bottom_arg} \
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
    path cosmic_coverage
    path pai3d_coverage
    val min_size
    val max_size
    path position_lists_dir

    output:
    path "posenrich_characterization.tsv", emit: results
    path "posenrich_leading_edge.tsv", emit: leading_edge

    script:
    // No permulation here: the exact Fisher p already corrects for the fixed
    // cutoff (no cutoff search to correct for), and the CAAS permulation null
    // (naive random relabeling, random cross-tree pairs) does not isolate the
    // phenotype at position level. Significance is instead dual-gated on
    // p_adj AND fold-enrichment (posenrich_padj_thr / posenrich_fold_thr) —
    // see posenrich_enrich.py's module docstring for why the fold gate exists.
    def annot_arg = annot_file.name != 'NO_FILE' ? "--annot-file ${annot_file}" : ""
    // cosmic_orthogroups/pai3d_orthogroups are GMTs derived from external,
    // incompletely-covered databases; restricting their background to genes
    // the database itself could ever annotate avoids diluting the test with
    // structurally-uncoverable genes (see build_position_gmt.py's coverage
    // file comment). Every other GMT keeps the full honest background.
    def cosmic_cov_arg = !(cosmic_coverage.name =~ /^NO_FILE/) ? "--cosmic-coverage ${cosmic_coverage}" : ""
    def pai3d_cov_arg  = !(pai3d_coverage.name =~ /^NO_FILE/) ? "--pai3d-coverage ${pai3d_coverage}" : ""
    // SCORING's own published position_lists/slice_{top,bottom,global}{25,10,5,1}.tsv
    // (scoring_compute.R) is posenrich's SOLE foreground source -- SCORING is a
    // mandatory upstream dependency, never optional, so this is passed
    // unconditionally (no NO_FILE-sentinel guard): if it's absent,
    // posenrich_enrich.py hard-fails with a clear message rather than silently
    // re-deriving its own ranking.
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_enrich.py \
        --obs-scores ${caas_file} \
        --gmt-dir gmts \
        --characterization ${characterization_layers} \
        --universe ${universe} \
        --background ${background_output} \
        ${annot_arg} \
        ${cosmic_cov_arg} \
        ${pai3d_cov_arg} \
        --position-lists-dir ${position_lists_dir} \
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
    publishDir path: "${params.outdir}/posenrich",
               mode: 'copy', overwrite: true,
               pattern: 'posenrich_overall_dotplot.tsv'

    input:
    path results
    path leading_edge
    // Position Characterisation (PrimateAI-3D + COSMIC + FADE validation, moved
    // here from the Scoring report): all optional, NO_FILE-sentinel-tolerant.
    // Section is skipped entirely (has_pos_char = FALSE) when position_scores
    // is absent.
    path position_scores
    path gene_scores
    path vep_primateai
    path vep_cosmic
    path genomic_info
    path fade_sites_top
    path fade_sites_bottom
    // SCORING's published position_lists/ dir -- same channel POSENRICH_RUN
    // already consumes as its foreground source; reused here so Position
    // Characterisation's own 10/5/1% membership (test_glob_*/test_top_*/
    // test_bot_* in 14.Position_enrichment_report.Rmd) reads the identical
    // published cutoff instead of re-deriving its own quantile().
    path position_lists
    // SCORING's fcs_stats.tsv (gene + flag_* columns) -- same file POSENRICH_RUN
    // already consumes as --annot-file, reused here for the Overall dotplot's
    // orthogonal-support composite (background flag rates via orthogonal_score.R).
    path fcs_stats
    // cleaned_background_main.txt -- the SAME universe_file 12.FCS_general_report.Rmd
    // uses (n_universe = length(universe)), so the Overall dotplot's gene-level
    // composite background is pinned to the identical universe everywhere it's
    // computed, not inferred from nrow(fcs_stats.tsv) (a different file that
    // usually, but isn't guaranteed to, agree in size).
    path cleaned_background

    output:
    path "14.Position_enrichment_report_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "posenrich_overall_dotplot.tsv",       emit: overall_dotplot,       optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def pos_scores_arg = (position_scores.name =~ /^NO_FILE/) ? 'NULL' : "'${position_scores}'"
    def gene_scores_arg = (gene_scores.name =~ /^NO_FILE/) ? 'NULL' : "'${gene_scores}'"
    def vep_pai_arg = (vep_primateai.name =~ /^NO_FILE/) ? 'NULL' : "'${vep_primateai}'"
    def vep_cosmic_arg = (vep_cosmic.name =~ /^NO_FILE/) ? 'NULL' : "'${vep_cosmic}'"
    def genomic_info_arg = (genomic_info.name =~ /^NO_FILE/) ? 'NULL' : "'${genomic_info}'"
    def fade_sites_top_arg = (fade_sites_top.name =~ /^NO_FILE/) ? 'NULL' : "'${fade_sites_top}'"
    def fade_sites_bottom_arg = (fade_sites_bottom.name =~ /^NO_FILE/) ? 'NULL' : "'${fade_sites_bottom}'"
    def fcs_stats_arg = (fcs_stats.name =~ /^NO_FILE/) ? 'NULL' : "'${fcs_stats}'"
    def universe_arg  = (cleaned_background.name =~ /^NO_FILE/) ? 'NULL' : "'${cleaned_background}'"
    def position_lists_arg = (position_lists.name =~ /^NO_/) ? 'NULL' : "'${position_lists}'"
    def render = """
        rmarkdown::render(
            '14.Position_enrichment_report.Rmd',
            params = list(
                results_file = '${results}',
                leading_edge_file = '${leading_edge}',
                traitname = '${traitname}',
                padj_thr = ${params.posenrich_padj_thr},
                fold_thr = ${params.posenrich_fold_thr},
                position_scores_file = ${pos_scores_arg},
                gene_scores_file     = ${gene_scores_arg},
                vep_primateai_file   = ${vep_pai_arg},
                vep_cosmic_file      = ${vep_cosmic_arg},
                genomic_info_file    = ${genomic_info_arg},
                fade_sites_top_file    = ${fade_sites_top_arg},
                fade_sites_bottom_file = ${fade_sites_bottom_arg},
                fcs_stats_file = ${fcs_stats_arg},
                universe_file  = ${universe_arg},
                position_lists_dir = ${position_lists_arg}
            ),
            output_file = '14.Position_enrichment_report_${traitname}.html'
        )
    """
    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "${render}"
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "${render}"
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
    pai3d_db
    cleaned_background
    caas_file
    position_lists_file     // SCORING's published position_lists/ dir -- mandatory, posenrich's sole foreground source
    background_output
    annot_file
    min_size
    max_size
    gene_scores_file        // optional: gene_scores.tsv (Position Characterisation)
    vep_primateai_file      // optional: PrimateAI-3D score TSV (Position Characterisation)
    vep_cosmic_file         // optional: COSMIC Mutant Census score TSV (Position Characterisation)
    genomic_info_file       // optional: gene genomic coords TSV (Position Characterisation)
    fade_sites_top_file     // optional: fade_sites_top.csv (FADE_top_sig position group)
    fade_sites_bottom_file  // optional: fade_sites_bottom.csv (FADE_bottom_sig position group)

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
        pai3d_db,
        cleaned_background,
        fade_sites_top_file,
        fade_sites_bottom_file
    )

    def cosmic_coverage_ch = POSENRICH_BUILD_GMT.out.cosmic_coverage.ifEmpty { file('NO_FILE_COSMIC_COV') }
    def pai3d_coverage_ch  = POSENRICH_BUILD_GMT.out.pai3d_coverage.ifEmpty { file('NO_FILE_PAI3D_COV') }

    POSENRICH_RUN(
        caas_file,
        POSENRICH_BUILD_GMT.out.gmts,
        POSENRICH_BUILD_GMT.out.charset,
        cleaned_background,
        background_output,
        annot_file,
        cosmic_coverage_ch,
        pai3d_coverage_ch,
        min_size,
        max_size,
        position_lists_file
    )

    POSENRICH_REPORT(
        POSENRICH_RUN.out.results,
        POSENRICH_RUN.out.leading_edge,
        caas_file,
        gene_scores_file,
        vep_primateai_file,
        vep_cosmic_file,
        genomic_info_file,
        fade_sites_top_file,
        fade_sites_bottom_file,
        position_lists_file,
        annot_file,
        cleaned_background
    )

    emit:
    results               = POSENRICH_RUN.out.results
    leading_edge          = POSENRICH_RUN.out.leading_edge
    report                = POSENRICH_REPORT.out.report
    overall_dotplot       = POSENRICH_REPORT.out.overall_dotplot
}
