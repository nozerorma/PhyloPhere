#!/usr/bin/env nextflow

/*
 * RER_GENE_LISTS
 * ──────────────
 * Extract FCS/AMI-ready gene lists from the RERconverge summary TSV:
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

    publishDir path: "${params.outdir}/rerconverge/gene_lists",
               mode: 'copy', overwrite: true,
               pattern: '*.txt'

    input:
    path summary_tsv
    path gene_scores

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

        # RER own directional significance. NOTE: gate_sig/gate_fdr are CAAS concepts
        # and are deliberately NOT set here — RER significance lives in flag_rer_acc/
        # flag_rer_decc. The CAAS gates (and FADE/accum/CAAS-directional) flags are
        # joined from the scoring gene-scores file below when one is provided.
        rer_sig       <- !is.na(pval_col) & pval_col < ${pval_thr}
        flag_rer_acc  <- rer_sig & df\$Rho > 0
        flag_rer_decc <- rer_sig & df\$Rho < 0

        fcs_stats <- data.frame(
            gene = df\$gene,
            score_global = sign(df\$Rho) * -log10(pmax(df\$P, 1e-300)),
            score_accelerating = ifelse(df\$Rho > 0, -log10(pmax(df\$P, 1e-300)), 0),
            score_decelerating = ifelse(df\$Rho < 0, -log10(pmax(df\$P, 1e-300)), 0),
            flag_rer_acc = flag_rer_acc,
            flag_rer_decc = flag_rer_decc,
            stringsAsFactors = FALSE
        )

        # Optional cross-module annotation: CAAS gates (gate_sig/gate_fdr), CAAS
        # directional scores (score_top/score_bottom -> top/bottom percentiles), FADE
        # and accumulation flags — imported from a scoring gene-scores TSV when given.
        # Absent columns are simply not joined, so the report omits them gracefully.
        gs_file <- "${gene_scores}"
        if (file.exists(gs_file) && !grepl("^NO_", basename(gs_file))) {
            gs <- tryCatch(read.delim(gs_file, stringsAsFactors = FALSE), error = function(e) NULL)
            if (!is.null(gs)) {
                if (!"gene" %in% colnames(gs)) {
                    gcol <- intersect(c("Gene","GENE"), colnames(gs))[1]
                    if (!is.na(gcol)) colnames(gs)[colnames(gs) == gcol] <- "gene"
                }
                want <- c("flag_gate_sig","flag_gate_fdr","score_top","score_bottom",
                          "flag_fade","flag_fade_top","flag_fade_bottom","flag_accum")
                want <- intersect(want, colnames(gs))
                if ("gene" %in% colnames(gs) && length(want) > 0) {
                    fcs_stats <- merge(fcs_stats, gs[, c("gene", want)], by = "gene", all.x = TRUE)
                    cat(sprintf("[RER_GENE_LISTS] joined %d cross-module col(s) from gene scores: %s\\n",
                        length(want), paste(want, collapse = ", ")))
                }
            }
        }

        write.table(fcs_stats, "fcs_stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

        cat(sprintf("[RER_GENE_LISTS] bg=%d  sig=%d  accel=%d  decel=%d\\n",
            nrow(df), nrow(sig), sum(sig\$Rho > 0), sum(sig\$Rho < 0)))
    '
    """
}
