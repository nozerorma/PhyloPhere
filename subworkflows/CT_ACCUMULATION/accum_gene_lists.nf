#!/usr/bin/env nextflow

/*
 * ACCUMULATION_GENE_LISTS
 * ───────────────────────
 * Extract AMI-ready gene lists and stats from CT_ACCUMULATION randomization output:
 *   background.txt                  — all genes tested by accumulation for this direction
 *   accumulation_${direction}_significant.txt — genes with CCT-combined FDR < threshold
 *   fcs_stats.tsv                   — ranked by -log10(CCT combined p)
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
        # Load the per-group scheme CSVs and combine their p-values with the
        # Cauchy Combination Test (CCT/ACAT), the same combiner scoring_compute.R
        # and 10.Accumulation_report.Rmd use on these same inputs.
        #
        # CCT is required here rather than any combiner that assumes independence:
        # one physical position can be a CAAS under several nested grouping
        # schemes, so the per-scheme counts — and therefore their p-values — are
        # positively correlated. CCT's Cauchy tail is heavy enough that its null
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

        # Merge into wide format
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
        cct_p <- apply(df_wide[, pval_cols, drop = FALSE], 1, function(ps) {
            ps <- ps[!is.na(ps)]
            if (length(ps) == 0) return(NA_real_)
            if (all(ps >= 1)) return(1.0)
            ps   <- pmin(pmax(ps, 1e-15), 1 - 1e-15)
            stat <- sum((1 / length(ps)) * tan((0.5 - ps) * pi))
            pcauchy(stat, lower.tail = FALSE)
        })

        # BH FDR only on genes with at least one observed CAAS (ActualCount > 0).
        # Genes with ActualCount=0 have p=1 by construction, not by test — they
        # must not enter the FDR denominator or be reported as significant.
        tested  <- !is.na(act_total) & act_total > 0
        fdr_q   <- rep(NA_real_, length(cct_p))
        if (any(tested))
            fdr_q[tested] <- p.adjust(cct_p[tested], method = 'BH')

        sig_mask  <- !is.na(fdr_q) & fdr_q < ${pval_thr}
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
            length(gene_syms), length(sig_genes), ${pval_thr}))
    "
    """
}
