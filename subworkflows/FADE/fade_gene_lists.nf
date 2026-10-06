#!/usr/bin/env nextflow
// fade_gene_lists.nf — Gene lists and per-gene statistics table from a FADE summary TSV.
// PhyloPhere | subworkflows/FADE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  FADE_GENE_LISTS: for one direction (top or bottom), extracts from the gene-level
 *  FADE summary the gene lists used by the AMI and FCS stages.
 *
 *  Consumes:  direction, fade_summary_<direction>.tsv (from FADE_REPORT)
 *  Produces:  background.txt (every gene tested by FADE in this direction),
 *             fade_<direction>_significant.txt (genes with max_bf >= fade_bf_threshold),
 *             fcs_stats.tsv (gene, score_fade = max_bf, flag_gate_sig = max_bf >= threshold)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Gene lists ─────────────────────────────────────────────────────────────────

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
