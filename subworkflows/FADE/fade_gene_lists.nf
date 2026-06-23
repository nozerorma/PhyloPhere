#!/usr/bin/env nextflow

/*
 * FADE_GENE_LISTS
 * ───────────────
 * Extract ORA-ready gene lists from a FADE summary TSV for one direction:
 *   background.txt                  — all genes tested by FADE for this direction
 *   fade_{direction}_significant.txt — genes with max_bf >= fade_bf_threshold
 *
 * Inputs
 * ──────
 *   direction   : val  — 'top' or 'bottom'
 *   summary_tsv : path — fade_summary_{direction}.tsv from FADE_REPORT
 *
 * Outputs
 * ───────
 *   direction  : val  — passed through for downstream routing
 *   gene_lists : path — both .txt files
 */

process FADE_GENE_LISTS {
    tag "fade_gene_lists|${direction}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir path: { "${params.outdir}/selection/fade/${direction}/gene_lists" },
               mode: 'copy', overwrite: true,
               pattern: '*.txt'

    input:
    val  direction
    path summary_tsv

    output:
    val  direction,        emit: direction
    path "*.txt",          emit: gene_lists
    path "fcs_stats.tsv",  emit: fcs_stats

    script:
    def bf_thr = params.fade_bf_threshold ?: 100
    """
    Rscript -e "
        df  <- read.delim('${summary_tsv}', stringsAsFactors = FALSE)
        writeLines(df\\\$gene, 'background.txt')
        sig <- df[!is.na(df\\\$max_bf) & df\\\$max_bf >= ${bf_thr}, ]
        writeLines(sig\\\$gene, 'fade_${direction}_significant.txt')
        
        fcs_stats <- data.frame(
            gene = df\\\$gene,
            score_fade = df\\\$max_bf,
            flag_gate_sig = !is.na(df\\\$max_bf) & df\\\$max_bf >= ${bf_thr}
        )
        write.table(fcs_stats, 'fcs_stats.tsv', sep = '\\t', row.names = FALSE, quote = FALSE)
        
        cat(sprintf('[FADE_GENE_LISTS] direction=${direction}  bg=%d  sig=%d\n',
            nrow(df), nrow(sig)))
    "
    """
}
