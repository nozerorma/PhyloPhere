#!/usr/bin/env nextflow

/*
 * ACCUMULATION_GENE_LISTS
 * ───────────────────────
 * Extract ORA/FCS-ready gene lists and stats from CT_ACCUMULATION randomization output:
 *   background.txt                  — all genes tested by accumulation for this direction
 *   accumulation_${direction}_significant.txt — genes with Fisher FDR < threshold
 *   fcs_stats.tsv                   — ranked by -log10(Fisher combined p)
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
        # Load all per-group scheme CSVs and compute Fisher's combined p-value.
        # Fisher's test combines independent per-group tests without the double-
        # counting that occurs when the same position satisfies multiple groupings.
        group_schemes <- c('us', 'gs1', 'gs2', 'gs3', 'gs4')
        all_dfs <- list()
        for (scheme in group_schemes) {
            pat   <- paste0('accumulation_${direction}_', scheme, '_aggregated_results\\\\.csv')
            files <- list.files('.', pattern = pat, full.names = TRUE)
            if (length(files) == 0) next
            d     <- read.csv(files[1], stringsAsFactors = FALSE)
            gcol  <- intersect(c('Gene', 'gene'), names(d))[1]
            pcol  <- grep('PValueEmpirical', names(d), value = TRUE)[1]
            acol  <- grep('ActualCount',     names(d), value = TRUE)[1]
            if (is.na(gcol) || is.na(pcol)) next
            all_dfs[[scheme]] <- data.frame(
                gene     = d[[gcol]],
                pval     = d[[pcol]],
                actcount = if (!is.na(acol)) d[[acol]] else 0L,
                stringsAsFactors = FALSE
            )
        }

        if (length(all_dfs) == 0)
            stop('No per-group aggregated results CSV files found for direction ${direction}')

        # Merge into wide format
        df_wide <- Reduce(
            function(a, b) merge(a, b, by = 'gene', all = TRUE),
            lapply(names(all_dfs), function(s)
                setNames(all_dfs[[s]], c('gene', paste0('pval_', s), paste0('actcount_', s))))
        )

        gene_syms <- df_wide\\\\$gene
        act_total <- rowSums(df_wide[, grep('^actcount_', names(df_wide)), drop = FALSE], na.rm = TRUE)
        pval_cols <- grep('^pval_', names(df_wide), value = TRUE)

        # Fisher χ² = -2Σln(p), df = 2k; p=1 for absent groups contributes 0
        fisher_p <- apply(df_wide[, pval_cols, drop = FALSE], 1, function(ps) {
            ps[is.na(ps)] <- 1.0
            ps <- pmax(ps, 1e-300)
            pchisq(-2 * sum(log(ps)), df = 2 * length(ps), lower.tail = FALSE)
        })

        # BH FDR only on genes with at least one observed CAAS (ActualCount > 0).
        # Genes with ActualCount=0 have p=1 by construction, not by test — they
        # must not enter the FDR denominator or be reported as significant.
        tested  <- !is.na(act_total) & act_total > 0
        fdr_q   <- rep(NA_real_, length(fisher_p))
        if (any(tested))
            fdr_q[tested] <- p.adjust(fisher_p[tested], method = 'BH')

        sig_mask  <- !is.na(fdr_q) & fdr_q < ${pval_thr}
        sig_genes <- gene_syms[sig_mask]

        writeLines(gene_syms, 'background.txt')
        writeLines(sig_genes, 'accumulation_${direction}_significant.txt')

        fcs_stats <- data.frame(
            gene               = gene_syms,
            score_accumulation = -log10(pmax(fisher_p, 1e-300)),
            fisher_p           = fisher_p,
            fdr_q              = fdr_q,
            flag_gate_sig      = sig_mask
        )
        write.table(fcs_stats, 'fcs_stats.tsv', sep = '\\t', row.names = FALSE, quote = FALSE)

        cat(sprintf('[ACCUMULATION_GENE_LISTS] direction=${direction}  bg=%d  sig=%d (FDR<%g)\\\\n',
            length(gene_syms), length(sig_genes), ${pval_thr}))
    "
    """
}
