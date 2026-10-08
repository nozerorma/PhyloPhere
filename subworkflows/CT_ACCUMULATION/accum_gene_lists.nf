#!/usr/bin/env nextflow

// accum_gene_lists.nf — Gene lists and a score table from the accumulation randomization output.
// PhyloPhere | subworkflows/CT_ACCUMULATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ACCUMULATION_GENE_LISTS: combines the per-scheme empirical p-values of each gene into
 *  one Cauchy-combined p-value (CCT), applies a BH FDR over the genes with at least one
 *  observed CAAS and writes the background and significant gene lists of one direction.
 *
 *  Consumes:  direction ('top', 'bottom' or 'all') and the accumulation_<direction>_<scheme>_aggregated_results.csv
 *             files of that direction (CT_ACCUMULATION_RANDOMIZE)
 *  Produces:  background.txt (all genes with a result for this direction),
 *             accumulation_<direction>_significant.txt (genes with CCT FDR below params.accumulation_fdr),
 *             fcs_stats.tsv (gene, score_accumulation = -log10 of the CCT p, accum_cct_p, fdr_q, flag_gate_sig)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Gene lists and CCT table ─────────────────────────────────────────────────

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
    def fdr_thr = params.accumulation_fdr ?: 0.1
    """
    Rscript -e "
        # Load the per-scheme CSVs and combine their p-values with the Cauchy
        # Combination Test (CCT/ACAT), the same combiner that scoring_compute.R
        # and 10.Accumulation_report.Rmd apply to these inputs.
        #
        # CCT is used instead of a combiner that assumes independence: one
        # physical position can be a CAAS under several nested grouping schemes,
        # so the per-scheme counts, and therefore their p-values, are positively
        # correlated. CCT's Cauchy tail is heavy enough that its null
        # holds under arbitrary dependence between the combined p-values, whereas
        # an independence-based combiner inflates the type-I rate several-fold in
        # this regime.
        group_schemes <- c('us', 'gs1', 'gs2', 'gs3', 'gs4')
        all_dfs <- list()
        for (scheme in group_schemes) {
            pat   <- paste0('accumulation_${direction}_', scheme, '_aggregated_results.csv')
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

        # Merge into wide format (one row per gene, one p-value and count column per scheme)
        df_wide <- Reduce(
            function(a, b) merge(a, b, by = 'gene', all = TRUE),
            lapply(names(all_dfs), function(s)
                setNames(all_dfs[[s]], c('gene', paste0('pval_', s), paste0('actcount_', s))))
        )

        gene_syms <- df_wide[['gene']]
        act_total <- rowSums(df_wide[, grep('^actcount_', names(df_wide)), drop = FALSE], na.rm = TRUE)
        pval_cols <- grep('^pval_', names(df_wide), value = TRUE)

        # CCT: T = sum(w_i * tan((0.5 - p_i) * pi)); p = pcauchy(T, lower.tail = FALSE).
        # Absent groups (NA) are dropped rather than coerced to 1: tan((0.5-1)*pi)
        # diverges to a large negative value, so feeding in a p of exactly 1 would
        # drag T down and cancel genuine signal from the groups that did report.
        # Rows where every group reports p=1 short-circuit to 1.
        # Weights follow the position score: US 0.5 and each GS 0.125, renormalized over the groups present.
        w_scheme <- c(us = 0.5, gs1 = 0.125, gs2 = 0.125, gs3 = 0.125, gs4 = 0.125)
        w_of     <- w_scheme[substring(pval_cols, 6)]
        cct_p <- apply(df_wide[, pval_cols, drop = FALSE], 1, function(ps) {
            valid <- !is.na(ps)
            ps <- ps[valid]
            if (length(ps) == 0) return(NA_real_)
            if (all(ps >= 1)) return(1.0)
            ps   <- pmin(pmax(ps, 1e-15), 1 - 1e-15)
            w    <- w_of[valid] / sum(w_of[valid])
            stat <- sum(w * tan((0.5 - ps) * pi))
            pcauchy(stat, lower.tail = FALSE)
        })

        # BH FDR only on genes with at least one observed CAAS (ActualCount > 0).
        # Genes with ActualCount = 0 have p = 1 by construction, not by test: they
        # must not enter the FDR denominator or be reported as significant.
        tested  <- !is.na(act_total) & act_total > 0
        fdr_q   <- rep(NA_real_, length(cct_p))
        if (any(tested))
            fdr_q[tested] <- p.adjust(cct_p[tested], method = 'BH')

        sig_mask  <- !is.na(fdr_q) & fdr_q < ${fdr_thr}
        sig_genes <- gene_syms[sig_mask]

        writeLines(gene_syms, 'background.txt')
        writeLines(sig_genes, 'accumulation_${direction}_significant.txt')

        fcs_stats <- data.frame(
            gene               = gene_syms,
            score_accumulation = -log10(pmax(cct_p, 1e-300)),
            accum_cct_p        = cct_p,
            fdr_q              = fdr_q,
            flag_gate_sig      = sig_mask
        )
        write.table(fcs_stats, 'fcs_stats.tsv', sep = '\\t', row.names = FALSE, quote = FALSE)

        cat(sprintf('[ACCUMULATION_GENE_LISTS] direction=${direction}  bg=%d  sig=%d (FDR<%g)\\\\n',
            length(gene_syms), length(sig_genes), ${fdr_thr}))
    "
    """
}
