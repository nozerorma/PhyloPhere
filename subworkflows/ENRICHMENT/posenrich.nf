#!/usr/bin/env nextflow
// posenrich.nf — Position-level enrichment of CAAS score magnitudes (Position-Level Path Sum Permulation).
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  POSENRICH: position-wise enrichment, not the gene-level FCS.
 *
 *  Builds position-level GMTs (Pfam, bins, orthogroups, COSMIC, PrimateAI-3D, UCR
 *  core/flank, FUBAR positive/purifying selection, FADE sites) plus the broad
 *  functional characterization layers, and tests them with Position-Level Path Sum
 *  Permulation (posenrich_enrich.py). Raw CAAS score magnitudes are summed per term
 *  and compared with a null over the whole tested-position background, per
 *  direction (global, top, bottom). The null is the CAAS permulation cycles
 *  (perm_pos_cycle_caas.tsv.gz), the same null fcs_enrich.R gives FCS's own Permsum
 *  test. Without a null, p_value, p_adj and perm_nes are NA and no term is
 *  significant. A term is significant when p_adj < posenrich_padj_thr with NES > 0.
 *
 *  Consumes:  position scores, SCORING's position_lists/, tested-position background,
 *             cleaned background, functional annotation inputs, CAAS permulation null
 *  Produces:  posenrich/ (posenrich_characterization.tsv, posenrich_leading_edge.tsv,
 *             HTML report), position_characterization.tsv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── GMT construction ───────────────────────────────────────────────────────────

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
    path "position_characterization.tsv", emit: position_char
    path "cosmic_coverage_genes.txt", optional: true, emit: cosmic_coverage
    path "pai3d_coverage_genes.txt", optional: true, emit: pai3d_coverage

    script:
    // Each optional input is absent as a uniquely named 'NO_FILE_*' sentinel (set in
    // workflows/enrichment.nf), not a shared 'NO_FILE': two path inputs staged under
    // the same name in one task directory collide, which happens whenever two or
    // more optional inputs are missing in the same run.
    def ensembl_arg = !(gene_ensembl_file.name =~ /^NO_FILE/) ? "--gene_ensembl_file ${gene_ensembl_file}" : ""
    def domain_arg  = !(domain_variability_file.name =~ /^NO_FILE/) ? "--domain_variability_file ${domain_variability_file}" : ""
    def ucr_arg     = !(ucr_positions_file.name =~ /^NO_FILE/) ? "--ucr_positions_file ${ucr_positions_file}" : ""
    def fubar_arg   = !(fubar_sites_file.name =~ /^NO_FILE/) ? "--fubar_sites_file ${fubar_sites_file}" : ""
    def egg_mem_arg = !(egg_members_file.name =~ /^NO_FILE/) ? "--egg_members_file ${egg_members_file}" : ""
    def egg_ann_arg = !(egg_annotations_file.name =~ /^NO_FILE/) ? "--egg_annotations_file ${egg_annotations_file}" : ""
    def cosmic_arg = !(cosmic_db.name =~ /^NO_FILE/) ? "--cosmic_db ${cosmic_db}" : ""
    def pai3d_arg  = !(pai3d_db.name =~ /^NO_FILE/) ? "--pai3d_db ${pai3d_db}" : ""
    def bg_arg     = !(cleaned_background.name =~ /^NO_FILE/) ? "--cleaned_background ${cleaned_background}" : ""
    // FADE_top_sig/FADE_bottom_sig position groups (from FADE_JSON_TO_CSV): the sites
    // with BF >= fade_bf_threshold, joined directly on Gene:Position, the coordinate
    // space of the CAAS Position column (no map_cache lookup).
    def fade_top_arg    = !(fade_sites_top_file.name    =~ /^NO_FILE/) ? "--fade_sites_top_file ${fade_sites_top_file}"       : ""
    def fade_bottom_arg = !(fade_sites_bottom_file.name =~ /^NO_FILE/) ? "--fade_sites_bottom_file ${fade_sites_bottom_file}" : ""
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/build_position_gmt.py \
        ${ensembl_arg} \
        ${domain_arg} \
        ${ucr_arg} \
        ${fubar_arg} \
        ${egg_mem_arg} \
        ${egg_ann_arg} \
        --map_dir ${map_dir} \
        ${cosmic_arg} \
        ${pai3d_arg} \
        ${bg_arg} \
        ${fade_top_arg} \
        ${fade_bottom_arg} \
        --output_dir .
    """
}


// ── Enrichment test ────────────────────────────────────────────────────────────

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
    path caas_cycle_null    // optional: perm_pos_cycle_caas.tsv.gz, the null (NO_FILE* sentinel for none)

    output:
    path "posenrich_characterization.tsv", emit: results
    path "posenrich_leading_edge.tsv", emit: leading_edge

    script:
    // Raw CAAS score magnitudes are summed per term and compared with a null. When
    // caas_cycle_null is supplied, its permulation cycles are the null; without one,
    // p_value, p_adj and perm_nes are NA and no term is significant. Significance is
    // p_adj < posenrich_padj_thr with NES > 0.
    def annot_arg = annot_file.name != 'NO_FILE' ? "--annot-file ${annot_file}" : ""
    // The cosmic_orthogroups and pai3d_orthogroups GMTs come from external databases
    // that do not cover every gene. Their background is restricted to the genes the
    // database can annotate, so that genes it cannot cover do not dilute the test
    // (coverage files written by build_position_gmt.py). Every other GMT keeps the
    // full tested-position background.
    def cosmic_cov_arg = !(cosmic_coverage.name =~ /^NO_FILE/) ? "--cosmic-coverage ${cosmic_coverage}" : ""
    def pai3d_cov_arg  = !(pai3d_coverage.name =~ /^NO_FILE/) ? "--pai3d-coverage ${pai3d_coverage}" : ""
    def caas_null_arg  = !(caas_cycle_null.name =~ /^NO_FILE/) ? "--caas-cycle-null ${caas_cycle_null}" : ""
    // SCORING's published position_lists/ (slice_{top,bottom,global}{25,10,5,1}.tsv,
    // from scoring_compute.R) is the only foreground source, so it is passed without
    // a NO_FILE guard. --background (caastools background.output, the tested
    // positions) is equally mandatory: without it the null background would collapse
    // to the scored positions only. When the CT stage is skipped, the caller supplies
    // it through params.posenrich_background_file.
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
        --n-perms ${params.posenrich_n_perms ?: 10000} \
        --perm-chunk-size ${params.posenrich_perm_chunk_size ?: 1000} \
        ${caas_null_arg} \
        --seed ${params.seed ?: 1998} \
        --padj-thr ${params.posenrich_padj_thr} \
        --output-dir .
    """
}


// ── Batched enrichment ─────────────────────────────────────────────────────────

process POSENRICH_PREP_NULL {
    label 'process_posenrich_prep_null'
    publishDir path: "${params.outdir}/posenrich",
               mode: 'copy', overwrite: true,
               enabled: params.publish_intermediates

    input:
    path caas_cycle_null    // optional: perm_pos_cycle_caas.tsv.gz, the null (NO_FILE* sentinel for none)

    output:
    path "caas_null_prepped.pkl", emit: prepped

    script:
    // Every POSENRICH_RUN_BATCHED task uses the same null (only the GMT term sets
    // differ), so the null is parsed once here instead of once per batch
    // (posenrich_prep_caas_null.py).
    def caas_null_arg = !(caas_cycle_null.name =~ /^NO_FILE/) ? "--caas-cycle-null ${caas_cycle_null}" : "--caas-cycle-null NO_FILE"
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_prep_caas_null.py \
        ${caas_null_arg} \
        --output caas_null_prepped.pkl
    """
}

process POSENRICH_RUN_BATCHED {
    tag "$batchID (${batchSize} GMTs)"
    label 'process_posenrich_batched'

    publishDir path: "${params.outdir}/posenrich/batches",
               mode: 'copy', overwrite: true,
               enabled: params.publish_intermediates

    input:
    tuple val(batchID), val(batchSize), val(includeCharacterization), path(gmtFiles, stageAs: 'gmts/*')
    path characterization_layers
    path caas_file
    path universe
    path background_output
    path annot_file
    path cosmic_coverage
    path pai3d_coverage
    val min_size
    val max_size
    path position_lists_dir
    path caas_null_prepped    // caas_null_prepped.pkl from POSENRICH_PREP_NULL (shared by all batches)

    output:
    path "posenrich_characterization.tsv", emit: results
    path "posenrich_leading_edge.tsv", emit: leading_edge

    script:
    def annot_arg = annot_file.name != 'NO_FILE' ? "--annot-file ${annot_file}" : ""
    def cosmic_cov_arg = !(cosmic_coverage.name =~ /^NO_FILE/) ? "--cosmic-coverage ${cosmic_coverage}" : ""
    def pai3d_cov_arg  = !(pai3d_coverage.name =~ /^NO_FILE/) ? "--pai3d-coverage ${pai3d_coverage}" : ""
    def caas_null_arg  = "--caas-null-prepped ${caas_null_prepped}"
    // The characterization layers are shared by the whole run, not split per GMT file:
    // giving them to every batch would append them once per batch in POSENRICH_CONCAT.
    // Only the batch with includeCharacterization = true receives them; an omitted
    // --characterization means "no characterization layer" for the others.
    def char_arg = includeCharacterization ? "--characterization ${characterization_layers}" : ""
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_enrich.py \
        --obs-scores ${caas_file} \
        --gmt-dir gmts \
        ${char_arg} \
        --universe ${universe} \
        --background ${background_output} \
        ${annot_arg} \
        ${cosmic_cov_arg} \
        ${pai3d_cov_arg} \
        --position-lists-dir ${position_lists_dir} \
        --min-size ${min_size} \
        --max-size ${max_size} \
        --n-perms ${params.posenrich_n_perms ?: 10000} \
        --perm-chunk-size ${params.posenrich_perm_chunk_size ?: 1000} \
        ${caas_null_arg} \
        --seed ${params.seed ?: 1998} \
        --padj-thr ${params.posenrich_padj_thr} \
        --output-dir .
    """
}


// ── Merge and report ───────────────────────────────────────────────────────────

process POSENRICH_CONCAT {
    tag "Concatenating POSENRICH batch outputs"
    publishDir "${params.outdir}/posenrich", mode: 'copy', overwrite: true

    input:
    path(characterization_files, stageAs: "characterization_*")
    path(leading_edge_files, stageAs: "leading_edge_*")

    output:
    path "posenrich_characterization.tsv", emit: results
    path "posenrich_leading_edge.tsv", emit: leading_edge

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mapfile -t char_files < <(find . -maxdepth 1 -name "characterization_*" ! -name ".*" | sort)
    cat "\${char_files[0]}" > posenrich_characterization.tsv
    for ((i=1; i<\${#char_files[@]}; i++)); do
        tail -n +2 "\${char_files[\$i]}" >> posenrich_characterization.tsv
    done

    mapfile -t le_files < <(find . -maxdepth 1 -name "leading_edge_*" ! -name ".*" | sort)
    cat "\${le_files[0]}" > posenrich_leading_edge.tsv
    for ((i=1; i<\${#le_files[@]}; i++)); do
        tail -n +2 "\${le_files[\$i]}" >> posenrich_leading_edge.tsv
    done
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
    publishDir path: "${params.outdir}/posenrich",
               mode: 'copy', overwrite: true,
               pattern: 'posenrich_leading_edge_summary.tsv'

    input:
    path results
    path leading_edge
    // Inputs of the Position Characterisation section (PrimateAI-3D, COSMIC and FADE
    // validation): all optional and NO_FILE-tolerant. The section is skipped when
    // position_scores is absent.
    path position_scores
    path gene_scores
    path vep_primateai
    path vep_cosmic
    path genomic_info
    path fade_sites_top
    path fade_sites_bottom
    // SCORING's position_lists/, the foreground source of the enrichment test. The
    // report reads the same 10/5/1% cutoffs instead of recomputing quantiles.
    path position_lists
    // SCORING's fcs_stats.tsv (gene + flag_* columns), the file POSENRICH_RUN takes as
    // --annot-file. It is passed to the report as fcs_stats_file.
    path fcs_stats
    // cleaned_background_main.txt, passed as the universe_file parameter: the tested-gene
    // background that restricts the PrimateAI-3D and COSMIC join sections.
    path cleaned_background

    output:
    path "14.Position_enrichment_report_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "posenrich_overall_dotplot.tsv",       emit: overall_dotplot,       optional: true
    path "posenrich_leading_edge_summary.tsv",  emit: leading_edge_summary,  optional: true

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
                position_scores_file = ${pos_scores_arg},
                gene_scores_file     = ${gene_scores_arg},
                vep_primateai_file   = ${vep_pai_arg},
                vep_cosmic_file      = ${vep_cosmic_arg},
                genomic_info_file    = ${genomic_info_arg},
                fade_sites_top_file    = ${fade_sites_top_arg},
                fade_sites_bottom_file = ${fade_sites_bottom_arg},
                fcs_stats_file = ${fcs_stats_arg},
                universe_file  = ${universe_arg},
                position_lists_dir = ${position_lists_arg},
                seed = '${params.seed ?: 1998}'
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


// ── Workflow ───────────────────────────────────────────────────────────────────

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
    position_lists_file     // SCORING's position_lists/ dir, mandatory (sole foreground source)
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
    caas_cycle_null_file    // optional: perm_pos_cycle_caas.tsv.gz, the null

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

    // Batch by GMT file: posenrich_enrich.py applies the BH correction within each
    // (direction, database) group, so spreading the GMT files over independent tasks
    // (each with its own memory limit and retries) does not change the statistics.
    // posenrich_batch_size = 1 runs everything in the single POSENRICH_RUN task.
    def posenrichBatchSize = (params.posenrich_batch_size ?: 1) as int
    def posenrich_results_ch
    def posenrich_leading_edge_ch
    if (posenrichBatchSize > 1) {
        def posenrichBatchCounter = 0
        def posenrich_batches = POSENRICH_BUILD_GMT.out.gmts
            .flatten()
            .collate(posenrichBatchSize)
            .map { batch ->
                def idx = ++posenrichBatchCounter
                def batchID = String.format('posenrich_batch_%05d', idx)
                tuple(batchID, batch.size(), idx == 1, batch)
            }

        // The single-item inputs (the take: channels and the coverage channels rebuilt
        // with .ifEmpty above) are queue channels, not value channels. Paired with the
        // many-item posenrich_batches channel, each would end POSENRICH_RUN_BATCHED
        // after its first batch. .collect().map { it[0] } turns each into a value
        // channel that every batch can read (same mechanism as CAAS_CORE in
        // caas_permulation.nf).
        def caas_file_bc           = caas_file.collect().map { it[0] }
        def cleaned_background_bc  = cleaned_background.collect().map { it[0] }
        def background_output_bc   = background_output.collect().map { it[0] }
        def annot_file_bc          = annot_file.collect().map { it[0] }
        def cosmic_coverage_bc     = cosmic_coverage_ch.collect().map { it[0] }
        def pai3d_coverage_bc      = pai3d_coverage_ch.collect().map { it[0] }
        def position_lists_file_bc = position_lists_file.collect().map { it[0] }

        // The null is parsed once for the whole run (POSENRICH_PREP_NULL).
        POSENRICH_PREP_NULL(caas_cycle_null_file)
        def caas_null_prepped_bc = POSENRICH_PREP_NULL.out.prepped.collect().map { it[0] }

        POSENRICH_RUN_BATCHED(
            posenrich_batches,
            POSENRICH_BUILD_GMT.out.charset,
            caas_file_bc,
            cleaned_background_bc,
            background_output_bc,
            annot_file_bc,
            cosmic_coverage_bc,
            pai3d_coverage_bc,
            min_size,
            max_size,
            position_lists_file_bc,
            caas_null_prepped_bc
        )

        POSENRICH_CONCAT(
            POSENRICH_RUN_BATCHED.out.results.collect(),
            POSENRICH_RUN_BATCHED.out.leading_edge.collect()
        )

        posenrich_results_ch = POSENRICH_CONCAT.out.results
        posenrich_leading_edge_ch = POSENRICH_CONCAT.out.leading_edge
    } else {
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
            position_lists_file,
            caas_cycle_null_file
        )

        posenrich_results_ch = POSENRICH_RUN.out.results
        posenrich_leading_edge_ch = POSENRICH_RUN.out.leading_edge
    }

    POSENRICH_REPORT(
        posenrich_results_ch,
        posenrich_leading_edge_ch,
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
    results                  = posenrich_results_ch
    leading_edge             = posenrich_leading_edge_ch
    report                   = POSENRICH_REPORT.out.report
    overall_dotplot          = POSENRICH_REPORT.out.overall_dotplot
    leading_edge_summary     = POSENRICH_REPORT.out.leading_edge_summary
    // Per-position Pfam domain/clan, UCR core/flank region and variability, and FUBAR
    // selection call, flattened for a Gene/Position join (read by 15.Comparison_report.Rmd
    // and by the SCORING report).
    position_characterization = POSENRICH_BUILD_GMT.out.position_char
}
