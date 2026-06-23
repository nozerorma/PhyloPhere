#!/usr/bin/env nextflow

/*
 * ACCUMULATION_GENE_LISTS
 * ───────────────────────
 * Extract ORA/FCS-ready gene lists and stats from CT_ACCUMULATION randomization output:
 *   background.txt                  — all genes tested by accumulation for this direction
 *   accumulation_${direction}_significant.txt — genes with PValueEmpirical_full_pool < threshold
 *   fcs_stats.tsv                   — ranked by -log10(PValueEmpirical_full_pool)
 *
 * Inputs
 * ──────
 *   direction : val  — 'top', 'bottom', or 'all'
 *   csv_files : path — list of all scheme CSV files for this direction
 *
 * Outputs
 * ───────
 *   direction  : val  — passed through for downstream routing
 *   gene_lists : path — both .txt files
 *   fcs_stats  : path — fcs_stats.tsv
 */

process ACCUMULATION_GENE_LISTS {
    tag "accum_gene_lists|${direction}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir path: { "${params.outdir}/accumulation/${direction}/gene_lists" },
               mode: 'copy', overwrite: true,
               pattern: '*.txt'

    input:
    val  direction
    path csv_files

    output:
    val  direction,        emit: direction
    path "*.txt",          emit: gene_lists
    path "fcs_stats.tsv",  emit: fcs_stats

    script:
    def pval_thr = params.accumulation_pval_threshold ?: 0.05
    """
    Rscript -e "
        # Find the full_pool CSV file
        files <- list.files('.', pattern = 'accumulation_${direction}_full_pool_aggregated_results\\\\.csv', full.names = TRUE)
        if (length(files) == 0) {
            # Fallback to any aggregated results CSV if full_pool is missing
            files <- list.files('.', pattern = '*_aggregated_results\\\\.csv', full.names = TRUE)
        }
        if (length(files) == 0) {
            stop('No aggregated results CSV files found for direction ${direction}')
        }
        
        csv_file <- files[1]
        df  <- read.csv(csv_file, stringsAsFactors = FALSE)
        
        # Standardise column names (Gene / gene)
        gcol <- intersect(c('Gene', 'gene'), colnames(df))[1]
        if (is.na(gcol)) stop('No Gene column found in CSV file')
        df\\\\$gene_sym <- df[[gcol]]
        
        # Find PValueEmpirical column
        pcol <- grep('PValueEmpirical', colnames(df), value = TRUE)[1]
        if (is.na(pcol)) stop('No PValueEmpirical column found in CSV file')
        
        writeLines(df\\\\$gene_sym, 'background.txt')
        
        pval_vec <- df[[pcol]]
        sig_mask <- !is.na(pval_vec) & pval_vec < ${pval_thr}
        sig_genes <- df\\\\$gene_sym[sig_mask]
        writeLines(sig_genes, 'accumulation_${direction}_significant.txt')
        
        # Build fcs_stats.tsv
        fcs_stats <- data.frame(
            gene = df\\\\$gene_sym,
            score_accumulation = -log10(pmax(pval_vec, 1e-300)),
            flag_gate_sig = sig_mask
        )
        write.table(fcs_stats, 'fcs_stats.tsv', sep = '\\t', row.names = FALSE, quote = FALSE)
        
        cat(sprintf('[ACCUMULATION_GENE_LISTS] direction=${direction}  bg=%d  sig=%d\\\\n',
            nrow(df), length(sig_genes)))
    "
    """
}
