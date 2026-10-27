#!/usr/bin/env nextflow
// rer_gene_lists.nf — Gene lists and per-gene statistics table from the RERconverge summary.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_GENE_LISTS: extracts from rerconverge_summary_<trait>.tsv (RER_REPORT) the gene
 *  lists used by the FCS and AMI stages, and the fcs_stats.tsv with the RER scores and
 *  flags. Nominal calls use p.perm < params.evidence_p_threshold. Missing empirical
 *  evidence yields no empirical calls; there is no analytic fallback. BH is separate.
 *
 *  Consumes:  summary TSV, optional SCORING gene-scores table (NO_* sentinel when absent)
 *  Produces:  background.txt (every gene tested by RERconverge),
 *             rer_significant.txt (p below the threshold),
 *             rer_accelerating.txt (legacy filename: nominal, Rho > 0),
 *             rer_decelerating.txt (legacy filename: nominal, Rho < 0),
 *             fcs_stats.tsv (gene, score_global, flag_rer_acc, flag_rer_decc, plus the
 *             cross-module columns found in the gene-scores table). score_global is
 *             sign(Rho) * -log10(P) with the parametric P of the correlation, the gene statistic
 *             of RERconverge (getStat) that its permulations carry as corStat; FCS ranks the
 *             genes by it and builds its null from corStat.
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
    def pval_thr = params.evidence_p_threshold != null ? params.evidence_p_threshold : 0.05
    def pval_col_name = 'p.perm'
    """
    Rscript -e '
        df  <- read.delim("${summary_tsv}", stringsAsFactors = FALSE)
        writeLines(df\$gene, "background.txt")

        pval_col_name <- "${pval_col_name}"
        if (pval_col_name == "p.perm" && (!"p.perm" %in% colnames(df) || all(is.na(df\$p.perm)))) {
            df\$p.perm <- NA_real_
            pval_col_name <- "p.perm"
        }
        pval_col <- df[[pval_col_name]]

        # Legacy filenames encode the sign of rho, not CAAS top/bottom. In continuous
        # diff mode this associates relative rates with descendant-minus-ancestor changes.
        sig <- df[!is.na(pval_col) & pval_col < ${pval_thr}, ]
        writeLines(sig\$gene,                 "rer_significant.txt")
        writeLines(sig\$gene[sig\$Rho > 0], "rer_accelerating.txt")
        writeLines(sig\$gene[sig\$Rho < 0], "rer_decelerating.txt")

        # Directional RER significance lives in flag_rer_acc and flag_rer_decc; gate_sig is
        # a CAAS flag and is not set here. The CAAS directional scores and the FADE
        # flags come from the scoring gene-scores file below, when given.
        rer_sig       <- !is.na(pval_col) & pval_col < ${pval_thr}
        flag_rer_acc  <- rer_sig & df\$Rho > 0
        flag_rer_decc <- rer_sig & df\$Rho < 0

        # Gene score: sign(Rho) * -log10(P), the parametric statistic of RERconverge (getStat).
        log_p    <- -log10(df\$P)

        fcs_stats <- data.frame(
            gene = df\$gene,
            score_global = sign(df\$Rho) * log_p,
            flag_rer_acc = flag_rer_acc,
            flag_rer_decc = flag_rer_decc,
            stringsAsFactors = FALSE
        )

        # Optional cross-module columns from the scoring gene-scores TSV (CAAS directional
        # scores and FADE flags). Columns absent from the file are not joined.
        gs_file <- "${gene_scores}"
        if (file.exists(gs_file) && !grepl("^NO_", basename(gs_file))) {
            gs <- tryCatch(read.delim(gs_file, stringsAsFactors = FALSE), error = function(e) NULL)
            if (!is.null(gs)) {
                if (!"gene" %in% colnames(gs)) {
                    gcol <- intersect(c("Gene","GENE"), colnames(gs))[1]
                    if (!is.na(gcol)) colnames(gs)[colnames(gs) == gcol] <- "gene"
                }
                want <- c("score_top","score_bottom",
                          "flag_fade","flag_fade_top","flag_fade_bottom")
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
