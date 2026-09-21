#!/usr/bin/env nextflow

/*
 * POSENRICH - Position-wise enrichment (Position-Level Path Sum Permulation)
 * ────────────────────────────────────────────────────────────────────
 * NOT the gene-level FCS. Builds position-level GMTs (Pfam, Bins, Orthogroups,
 * COSMIC, UCR core/flank, positive/purifying selection) plus the broad
 * functional characterization layers, and tests them with Position-Level Path
 * Sum Permulation (posenrich_enrich.py): raw CAAS score magnitudes are summed
 * per term and compared against a label-permuted null across the full honest
 * ~1.47M position background, per direction (global/top/bottom). Replaces an
 * earlier fixed-cutoff Fisher-exact design, which had power to flag
 * biologically negligible deviations as significant at that background size.
 * Significance is p_adj < posenrich_padj_thr with NES > 0.
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
    path "position_characterization.tsv", emit: position_char
    path "cosmic_coverage_genes.txt", optional: true, emit: cosmic_coverage
    path "pai3d_coverage_genes.txt", optional: true, emit: pai3d_coverage

    script:
    // Each optional input's absent-value sentinel is a uniquely-named 'NO_FILE_*'
    // path (see workflows/enrichment.nf) rather than a shared literal 'NO_FILE':
    // staging two path inputs under the identical filename in one task directory
    // is a Nextflow "input file name collision", which is exactly what happens
    // if two or more of these optional inputs are absent in the same run.
    def cosmic_arg = !(cosmic_db.name =~ /^NO_FILE/) ? "--cosmic_db ${cosmic_db}" : ""
    def pai3d_arg  = !(pai3d_db.name =~ /^NO_FILE/) ? "--pai3d_db ${pai3d_db}" : ""
    def bg_arg     = !(cleaned_background.name =~ /^NO_FILE/) ? "--cleaned_background ${cleaned_background}" : ""
    // FADE_top_sig/FADE_bottom_sig position group (§ FADE_JSON_TO_CSV) - the
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
    path caas_cycle_null    // optional: perm_pos_cycle_caas.tsv.gz -> p.perm

    output:
    path "posenrich_characterization.tsv", emit: results
    path "posenrich_leading_edge.tsv", emit: leading_edge

    script:
    // Position-Level Path Sum Permulation (posenrich_enrich.py): raw CAAS score
    // magnitudes are summed per term and compared against a label-permuted null
    // built from posenrich_n_perms permutations, AND (when caas_cycle_null is
    // supplied) against the CAAS permulation null's real cycles -> p.perm.
    // Significance is p_adj < posenrich_padj_thr with NES > 0, additionally
    // gated on p.perm < posenrich_p_perm_thr whenever p.perm is available.
    def annot_arg = annot_file.name != 'NO_FILE' ? "--annot-file ${annot_file}" : ""
    // cosmic_orthogroups/pai3d_orthogroups are GMTs derived from external,
    // incompletely-covered databases; restricting their background to genes
    // the database itself could ever annotate avoids diluting the test with
    // structurally-uncoverable genes (see build_position_gmt.py's coverage
    // file comment). Every other GMT keeps the full honest background.
    def cosmic_cov_arg = !(cosmic_coverage.name =~ /^NO_FILE/) ? "--cosmic-coverage ${cosmic_coverage}" : ""
    def pai3d_cov_arg  = !(pai3d_coverage.name =~ /^NO_FILE/) ? "--pai3d-coverage ${pai3d_coverage}" : ""
    def caas_null_arg  = !(caas_cycle_null.name =~ /^NO_FILE/) ? "--caas-cycle-null ${caas_cycle_null}" : ""
    // SCORING's own published position_lists/slice_{top,bottom,global}{25,10,5,1}.tsv
    // (scoring_compute.R) is posenrich's SOLE foreground source, SCORING is a
    // mandatory upstream dependency, never optional, so this is passed
    // unconditionally (no NO_FILE-sentinel guard): if it's absent,
    // posenrich_enrich.py hard-fails with a clear message rather than silently
    // re-deriving its own ranking.
    // --background (caastools background.output, tested positions) is equally
    // mandatory: without it the permulation null's background collapses to
    // scored positions only, invalidating the test. Also passed unconditionally;
    // posenrich_enrich.py hard-fails on any NO_FILE* sentinel or missing path.
    // In _complete runs (CT skipped), the caller must supply it via
    // params.posenrich_background_file (wired in run_single.sh.j2's
    // reuse_exploratory block from the _exploratory caastools/background.output).
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
        --p-perm-thr ${params.posenrich_p_perm_thr ?: 0.025} \
        --seed ${params.seed ?: 1998} \
        --padj-thr ${params.posenrich_padj_thr} \
        --output-dir .
    """
}

process POSENRICH_PREP_NULL {
    label 'process_posenrich_prep_null'
    publishDir path: "${params.outdir}/posenrich",
               mode: 'copy', overwrite: true,
               enabled: params.publish_intermediates

    input:
    path caas_cycle_null    // optional: perm_pos_cycle_caas.tsv.gz -> p.perm; NO_FILE* sentinel to skip

    output:
    path "caas_null_prepped.pkl", emit: prepped

    script:
    // perm_pos_cycle_caas.tsv.gz is broadcast identically to every
    // POSENRICH_RUN_BATCHED task (same null, only the GMT term sets differ
    // per batch) -- parsing it once here, instead of once per batch task,
    // is the whole point of this process. See posenrich_prep_caas_null.py.
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
    path caas_null_prepped    // caas_null_prepped.pkl from POSENRICH_PREP_NULL -> p.perm (parsed once for the whole run, not once per batch)

    output:
    path "posenrich_characterization.tsv", emit: results
    path "posenrich_leading_edge.tsv", emit: leading_edge

    script:
    def annot_arg = annot_file.name != 'NO_FILE' ? "--annot-file ${annot_file}" : ""
    def cosmic_cov_arg = !(cosmic_coverage.name =~ /^NO_FILE/) ? "--cosmic-coverage ${cosmic_coverage}" : ""
    def pai3d_cov_arg  = !(pai3d_coverage.name =~ /^NO_FILE/) ? "--pai3d-coverage ${pai3d_coverage}" : ""
    def caas_null_arg  = "--caas-null-prepped ${caas_null_prepped}"
    // The characterization layer (Pfam/UCR/FUBAR/...) is one source shared
    // across the whole run, not split per GMT file — passing it to every batch
    // would re-run and re-append it N times once POSENRICH_CONCAT merges the
    // batches. Only the designated batch (includeCharacterization=true) gets it;
    // posenrich_enrich.py's read_charset() already treats a missing/omitted
    // --characterization as "no characterization layer" for every other batch.
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
        --p-perm-thr ${params.posenrich_p_perm_thr ?: 0.025} \
        --seed ${params.seed ?: 1998} \
        --padj-thr ${params.posenrich_padj_thr} \
        --output-dir .
    """
}

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
    // SCORING's published position_lists/ dir, same channel POSENRICH_RUN
    // already consumes as its foreground source; reused here so Position
    // Characterisation's own 10/5/1% membership (test_glob_*/test_top_*/
    // test_bot_* in 14.Position_enrichment_report.Rmd) reads the identical
    // published cutoff instead of re-deriving its own quantile().
    path position_lists
    // SCORING's fcs_stats.tsv (gene + flag_* columns), same file POSENRICH_RUN
    // already consumes as --annot-file. Kept as an input for pipeline-wiring
    // compatibility; 14.Position_enrichment_report.Rmd no longer reads it
    // directly (the Overall dotplot's cross-module composite that needed it
    // was removed).
    path fcs_stats
    // cleaned_background_main.txt, staged into the Rmd's universe_file param,
    // the honest tested-gene background for the PrimateAI-3D/COSMIC join
    // sections' shared_genes_pai/shared_genes_cosmic restriction.
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
    position_lists_file     // SCORING's published position_lists/ dir, mandatory, posenrich's sole foreground source
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
    caas_cycle_null_file    // optional: perm_pos_cycle_caas.tsv.gz -> p.perm

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

    // Batch by GMT/db source: posenrich_enrich.py already scopes BH correction
    // per (direction, db) group, so splitting the GMT-file loop across
    // independent Nextflow tasks (one SLURM job each, its own memory ceiling
    // and its own retry) changes nothing statistically. batch_size=1 keeps the
    // original monolithic POSENRICH_RUN path untouched.
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
                def batchID = sprintf('posenrich_batch_%05d', idx)
                tuple(batchID, batch.size(), idx == 1, batch)
            }

        // caas_file/cleaned_background/background_output/annot_file/
        // position_lists_file are take: params (single item, but plain queue
        // channels, not genuine Nextflow value channels -- see
        // caas_permulation.nf's CAAS_PERMS_DISAMBIGUATE_BATCHED fix for the
        // full mechanism); cosmic_coverage_ch/pai3d_coverage_ch are similarly
        // one-item channels rebuilt via .ifEmpty() above. Paired positionally
        // against the many-item posenrich_batches channel, any of these would
        // silently truncate POSENRICH_RUN_BATCHED to its first batch once
        // exhausted. .first() makes each a proper reusable/broadcastable
        // channel; no-op for anything that was already a value channel.
        def caas_file_bc           = caas_file.first()
        def cleaned_background_bc  = cleaned_background.first()
        def background_output_bc   = background_output.first()
        def annot_file_bc          = annot_file.first()
        def cosmic_coverage_bc     = cosmic_coverage_ch.first()
        def pai3d_coverage_bc      = pai3d_coverage_ch.first()
        def position_lists_file_bc = position_lists_file.first()

        // Parsed once for the whole run here, instead of once per batch task
        // inside POSENRICH_RUN_BATCHED -- see POSENRICH_PREP_NULL / posenrich_prep_caas_null.py.
        POSENRICH_PREP_NULL(caas_cycle_null_file)
        def caas_null_prepped_bc = POSENRICH_PREP_NULL.out.prepped.first()

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
    // Per-position PFAM domain/clan, UCR core/flank region + per-position
    // variability, and FUBAR selection call, flattened for a direct
    // Gene/Position join (15.Comparison_report.Rmd's Interesting Positions table).
    position_characterization = POSENRICH_BUILD_GMT.out.position_char
}
