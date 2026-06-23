#!/usr/bin/env nextflow

/*
 * RER_GENE_LISTS
 * ──────────────
 * Extract ORA-ready gene lists from the RERconverge summary TSV:
 *   background.txt           — all genes tested by RERConverge
 *   rer_significant.txt      — p.perm < threshold & |Rho| >= rho_threshold
 *   rer_accelerating.txt     — significant & Rho > 0
 *   rer_decelerating.txt     — significant & Rho < 0
 *
 * Inputs
 * ──────
 *   summary_tsv : path — rerconverge_summary_{trait}.tsv from RER_REPORT
 *
 * Outputs
 * ───────
 *   gene_lists : path — all four .txt files (collected)
 */

process RER_GENE_LISTS {
    tag "rer_gene_lists|${params.traitname}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/RERConverge/gene_lists",
               mode: 'copy', overwrite: true,
               pattern: '*.txt'

    input:
    path summary_tsv

    output:
    path "*.txt",          emit: gene_lists
    path "fcs_stats.tsv",  emit: fcs_stats

    script:
    def pval_thr = params.rer_pval_threshold != null ? params.rer_pval_threshold : 0.05
    def pval_col_name = params.rer_pval_column != null ? params.rer_pval_column : 'p.perm'
    """
    Rscript -e '
        df  <- read.delim("${summary_tsv}", stringsAsFactors = FALSE)
        writeLines(df\$gene, "background.txt")
        
        pval_col_name <- "${pval_col_name}"
        if (pval_col_name == "p.perm" && (!"p.perm" %in% colnames(df) || all(is.na(df\$p.perm)))) {
            pval_col_name <- "p.adj"
        }
        pval_col <- df[[pval_col_name]]
        
        sig <- df[!is.na(pval_col) & pval_col < ${pval_thr}, ]
        writeLines(sig\$gene,                 "rer_significant.txt")
        writeLines(sig\$gene[sig\$Rho > 0], "rer_accelerating.txt")
        writeLines(sig\$gene[sig\$Rho < 0], "rer_decelerating.txt")
        
        flag_gate_sig <- !is.na(pval_col) & pval_col < ${pval_thr}
        flag_gate_fdr <- !is.na(df\$p.adj) & df\$p.adj < ${pval_thr}
        flag_rer_acc  <- flag_gate_sig & df\$Rho > 0
        flag_rer_decc <- flag_gate_sig & df\$Rho < 0
        
        fcs_stats <- data.frame(
            gene = df\$gene,
            score_global = sign(df\$Rho) * -log10(pmax(df\$P, 1e-300)),
            score_accelerating = ifelse(df\$Rho > 0, -log10(pmax(df\$P, 1e-300)), 0),
            score_decelerating = ifelse(df\$Rho < 0, -log10(pmax(df\$P, 1e-300)), 0),
            flag_gate_sig = flag_gate_sig,
            flag_gate_fdr = flag_gate_fdr,
            flag_rer_acc = flag_rer_acc,
            flag_rer_decc = flag_rer_decc
        )
        write.table(fcs_stats, "fcs_stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
        
        cat(sprintf("[RER_GENE_LISTS] bg=%d  sig=%d  accel=%d  decel=%d\\n",
            nrow(df), nrow(sig), sum(sig\$Rho > 0), sum(sig\$Rho < 0)))
    '
    """
}
