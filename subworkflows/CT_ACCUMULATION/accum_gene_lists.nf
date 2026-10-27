#!/usr/bin/env nextflow

// accum_gene_lists.nf — Gene lists from the accumulation test output.
// PhyloPhere | subworkflows/CT_ACCUMULATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ACCUMULATION_GENE_LISTS: takes the exact hypergeometric p-values of each gene,
 *  applies a BH FDR over the genes with at least one
 *  observed CAAS and writes the background and significant gene lists of one direction.
 *
 *  Consumes:  direction ('top', 'bottom' or 'all') and the accumulation_<direction>_results.csv
 *             file of that direction (CT_ACCUMULATION_TEST)
 *  Produces:  background.txt (all genes with a result for this direction),
 *             accumulation_<direction>_significant.txt (genes with FDR below params.accumulation_fdr)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Gene lists and table ─────────────────────────────────────────────────────

process ACCUMULATION_GENE_LISTS {
    tag "accum_gene_lists|${direction}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir path: { "${params.outdir}/posenrich/accumulation/${direction}/gene_lists" },
               mode: 'copy', overwrite: true,
               pattern: '*.txt'

    input:
    val  direction
    path csv_files

    output:
    val  direction,        emit: direction
    path "*.txt",          emit: gene_lists

    script:
    def fdr_thr = params.accumulation_fdr ?: 0.1
    """
    Rscript -e "
        # Load the aggregated results CSV and apply BH FDR over genes with >= 1 observed CAAS.
        pat   <- 'accumulation_${direction}_results.csv'
        files <- list.files('.', pattern = pat, full.names = TRUE)
        if (length(files) == 0)
            stop('No aggregated results CSV found for direction ${direction}')
        d     <- read.csv(files[1], stringsAsFactors = FALSE)
        gcol  <- intersect(c('Gene', 'gene'), names(d))[1]
        pcol  <- grep('^PValue\$', names(d), value = TRUE)[1]
        acol  <- grep('ActualCount',     names(d), value = TRUE)[1]
        if (is.na(gcol) || is.na(pcol))
            stop('Missing Gene or PValue column in the accumulation results')

        gene_syms <- d[[gcol]]
        pvals     <- as.numeric(d[[pcol]])
        act_total <- if (!is.na(acol)) as.integer(d[[acol]]) else rep(0L, length(pvals))

        # BH FDR only on genes with at least one observed CAAS (ActualCount > 0).
        # Genes with ActualCount = 0 have p = 1 by construction, not by test: they
        # must not enter the FDR denominator or be reported as significant.
        tested  <- !is.na(act_total) & act_total > 0
        fdr_q   <- rep(NA_real_, length(pvals))
        if (any(tested))
            fdr_q[tested] <- p.adjust(pvals[tested], method = 'BH')

        sig_mask  <- !is.na(fdr_q) & fdr_q < ${fdr_thr}
        sig_genes <- gene_syms[sig_mask]

        writeLines(gene_syms, 'background.txt')
        writeLines(sig_genes, 'accumulation_${direction}_significant.txt')

        cat(sprintf('[ACCUMULATION_GENE_LISTS] direction=${direction}  bg=%d  sig=%d (FDR<%g)\\\\n',
            length(gene_syms), length(sig_genes), ${fdr_thr}))
    "
    """
}
