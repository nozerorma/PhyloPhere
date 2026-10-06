#!/usr/bin/env nextflow
// rer_gene_lists.nf — Gene lists and per-gene statistics table from the RERconverge summary.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_GENE_LISTS: extracts from rerconverge_summary_<trait>.tsv (RER_REPORT) the gene
 *  lists used by the FCS and AMI stages, and the fcs_stats.tsv with the RER scores and
 *  flags. The significance column is params.rer_pval_column (default p.perm, replaced
 *  by p.adj when p.perm is absent or all NA) at params.rer_pval_threshold.
 *
 *  Consumes:  summary TSV, optional SCORING gene-scores table (NO_* sentinel when absent)
 *  Produces:  background.txt (every gene tested by RERconverge),
 *             rer_significant.txt (p below the threshold),
 *             rer_accelerating.txt (significant, Rho > 0),
 *             rer_decelerating.txt (significant, Rho < 0),
 *             fcs_stats.tsv (gene, score_global, score_accelerating, score_decelerating,
 *             flag_rer_acc, flag_rer_decc, plus the cross-module columns found in the
 *             gene-scores table)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Gene lists ─────────────────────────────────────────────────────────────────

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

        # Directional RER significance lives in flag_rer_acc and flag_rer_decc; gate_sig is
        # a CAAS flag and is not set here. The CAAS directional scores and the FADE and
        # accumulation flags come from the scoring gene-scores file below, when given.
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

        # Optional cross-module columns from the scoring gene-scores TSV (CAAS directional
        # scores, FADE and accumulation flags). Columns absent from the file are not joined.
        gs_file <- "${gene_scores}"
        if (file.exists(gs_file) && !grepl("^NO_", basename(gs_file))) {
            gs <- tryCatch(read.delim(gs_file, stringsAsFactors = FALSE), error = function(e) NULL)
            if (!is.null(gs)) {
                if (!"gene" %in% colnames(gs)) {
                    gcol <- intersect(c("Gene","GENE"), colnames(gs))[1]
                    if (!is.na(gcol)) colnames(gs)[colnames(gs) == gcol] <- "gene"
                }
                want <- c("score_top","score_bottom",
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
