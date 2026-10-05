"""The vectorized permulation null of the FCS computes the statistic the observed side computes.

The observed statistic of a set is the Wilcoxon AUC - 0.5 of its genes against the rest of the annotated genes
(RERconverge::fastwilcoxGMT). The null recomputes it for every permulation column in one sparse product
(fcs_null_enrichstat_vectorized). Here it is checked against the definition with base R (heavy ties, structural zeros,
NA); check_fcs_null_equivalence.R checks it against RERconverge itself where that package is installed.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "subworkflows/ENRICHMENT/local/src/fcs_enrich.R"
CHECK = Path(__file__).resolve().parent / "check_fcs_null_equivalence.R"

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")

NAMES = ("fcs_colranks", "fcs_membership_matrix", "fcs_null_enrichstat_vectorized", "fcs_permpvalenrich_vectorized")
PRELUDE = f'''
suppressPackageStartupMessages(library(Matrix))
for (e in parse("{SCRIPT}")) {{
  if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c({", ".join(repr(n).replace("'", '"') for n in NAMES)})) eval(e)
}}
'''


def _r(code):
    out = subprocess.run(["Rscript", "-e", PRELUDE + code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-1500:]
    return out.stdout.strip()


# a null matrix as the CAAS null looks: mostly structural zeros, rounded scores (ties), optionally NA cells
DEFINITION = '''
set.seed(%(seed)d)
G <- 160; N <- 24; genes <- sprintf("g%%03d", seq_len(G))
cs <- matrix(0, G, N, dimnames = list(genes, NULL))
nz <- matrix(runif(G * N) < 0.2, G, N); cs[nz] <- round(runif(sum(nz)), 1)
if (%(na)s) cs[matrix(runif(G * N) < 0.05, G, N)] <- NA
mk <- function(prefix, k, sizes) { gs <- lapply(seq_len(k), function(i) sample(genes[1:140], sample(sizes, 1)))
  names(gs) <- sprintf("%%s%%02d", prefix, seq_len(k)); list(genesets = gs, geneset.names = names(gs)) }
gmts <- list(dbA = mk("A", 14, c(3, 8, 10, 15, 30, 70)), dbB = mk("B", 9, c(10, 25, 50, 90)))
real <- lapply(gmts, function(g) data.frame(pval = rep(NA_real_, length(g$geneset.names)), stat = NA_real_, row.names = g$geneset.names))
num_g <- 10; max_g <- 60
vec <- fcs_null_enrichstat_vectorized(cs, gmts, real, num_g, max_g)
worst <- 0; cells <- 0; na_mismatch <- 0; computed <- 0
for (db in names(gmts)) {
  # the observed side hands fastwilcoxGMT the GMT without the sets over max_g, and its background is the union of the genes
  # of the sets it is given
  gs <- gmts[[db]]$genesets; kept <- gs[vapply(gs, function(set) length(intersect(set, genes)) <= max_g, logical(1))]
  genes_db <- intersect(unique(unlist(kept)), genes)
  for (s in names(gs)) for (j in seq_len(N)) {
    v <- cs[genes_db, j]; ok <- !is.na(v)
    inset <- genes_db %%in%% gs[[s]]
    x <- v[ok & inset]; y <- v[ok & !inset]
    n1 <- length(x); n2 <- length(y)
    ref <- if (!(s %%in%% names(kept)) || n1 < num_g || n1 > max_g || n2 <= 2) NA_real_ else
      unname(suppressWarnings(wilcox.test(x, y)$statistic)) / (n1 * n2) - 0.5
    got <- vec[[db]][s, j]
    cells <- cells + 1
    if (is.na(ref) != is.na(got)) na_mismatch <- na_mismatch + 1
    else if (!is.na(ref)) { computed <- computed + 1; worst <- max(worst, abs(ref - got)) }
  }
}
cat(cells, computed, na_mismatch, format(worst, digits = 3))
'''


@pytest.mark.parametrize("seed,na", [(1, "FALSE"), (2, "FALSE"), (3, "TRUE"), (4, "TRUE")], ids=["no_na_1", "no_na_2", "with_na_1", "with_na_2"])
def test_the_vectorized_statistic_is_the_wilcoxon_auc_minus_one_half_of_each_set_in_each_column(seed, na):
    cells, computed, na_mismatch, worst = _r(DEFINITION % {"seed": seed, "na": na}).split()
    assert int(cells) == 23 * 24 and int(computed) > 100        # the check is exercised: many sets pass the size gates
    assert int(na_mismatch) == 0 and float(worst) < 1e-12


def test_sets_below_the_size_gates_are_undefined_in_every_column():
    out = _r('''
cs <- matrix(c(rep(0, 20), seq(0.1, 2, by = 0.1)), 20, 2, dimnames = list(sprintf("g%02d", 1:20), NULL))
gs <- list(small = sprintf("g%02d", 1:3), fine = sprintf("g%02d", 1:8), huge = sprintf("g%02d", 1:19), pad = sprintf("g%02d", 9:20))
gmts <- list(d = list(genesets = gs, geneset.names = names(gs)))
real <- list(d = data.frame(pval = rep(NA_real_, 4), stat = NA_real_, row.names = names(gs)))
v <- fcs_null_enrichstat_vectorized(cs, gmts, real, num_g = 5, max_g = 15)$d
cat(paste(is.na(v[, 1]), collapse = ","), "|", paste(is.na(v[, 2]), collapse = ","))
''')
    assert out == "TRUE,FALSE,TRUE,FALSE | TRUE,FALSE,TRUE,FALSE"     # small: n1 < num_g; huge: over max_g, dropped; pad closes the universe


def test_the_empirical_p_counts_the_null_cells_at_least_as_extreme_with_a_pseudocount_and_skips_the_undefined_ones():
    out = _r('''
real <- list(d = data.frame(pval = rep(NA_real_, 3), stat = c(0.2, -0.2, NA), row.names = c("s1", "s2", "s3")))
null <- list(d = matrix(c(0.1, 0.3, NA, 0.2,   -0.3, 0.1, 0.2, NA,   0.0, 0.0, 0.0, 0.0), nrow = 3, byrow = TRUE,
                        dimnames = list(c("s1", "s2", "s3"), NULL)))
g <- fcs_permpvalenrich_vectorized(real, null, "greater")$d
t <- fcs_permpvalenrich_vectorized(real, null, "two.sided")$d
cat(g, "|", t)
''')
    greater, two_sided = out.split("|")
    g = [float(x) if x != "NA" else None for x in greater.split()]
    t = [float(x) if x != "NA" else None for x in two_sided.split()]
    # greater: s1 (obs 0.2) is reached by {0.3, 0.2} of the 3 defined cells; s2 (obs -0.2) by {0.1, 0.2} (-0.3 is below it)
    assert g[0] == pytest.approx((2 + 1) / (3 + 1)) and g[1] == pytest.approx((2 + 1) / (3 + 1)) and g[2] is None
    # two-sided: s1 |null| >= 0.2 -> {0.3, 0.2}; s2 |null| >= 0.2 -> {0.3, 0.2}
    assert t[0] == pytest.approx((2 + 1) / (3 + 1)) and t[1] == pytest.approx((2 + 1) / (3 + 1)) and t[2] is None


def test_the_rest_of_the_annotated_genes_needs_more_than_two_genes():
    out = _r('''
cs <- matrix(seq(0.05, 1.2, length.out = 12), 12, 1, dimnames = list(sprintf("g%02d", 1:12), NULL))
gs <- list(pad = c("g11", "g12"), tight = sprintf("g%02d", 1:10), ok = sprintf("g%02d", 1:9))   # annotated universe: 12 genes
gmts <- list(d = list(genesets = gs, geneset.names = names(gs)))
real <- list(d = data.frame(pval = rep(NA_real_, 3), stat = NA_real_, row.names = names(gs)))
v <- fcs_null_enrichstat_vectorized(cs, gmts, real, num_g = 5, max_g = 0)$d
cat(paste(is.na(v[, 1]), collapse = ","))
''')
    assert out == "TRUE,TRUE,FALSE"       # pad: n1 < num_g; tight: n2 = 2; ok: n2 = 3


def test_with_missing_cells_the_rest_of_the_genes_is_counted_per_column():
    out = _r('''
cs <- matrix(rep(seq(0.05, 1.2, length.out = 12), 2), 12, 2, dimnames = list(sprintf("g%02d", 1:12), NULL))
cs["g12", 2] <- NA                                                       # column 2 has 11 defined genes
gs <- list(pad = c("g10", "g11", "g12"), A = sprintf("g%02d", 1:9))   # annotated universe: g01..g12
gmts <- list(d = list(genesets = gs, geneset.names = names(gs)))
real <- list(d = data.frame(pval = rep(NA_real_, 2), stat = NA_real_, row.names = names(gs)))
v <- fcs_null_enrichstat_vectorized(cs, gmts, real, num_g = 5, max_g = 0)$d
cat(paste(is.na(v["A", ]), collapse = ","))
''')
    assert out == "FALSE,TRUE"     # A: 9 genes against 3 others in column 1, against 2 in column 2 (undefined)


def test_the_genes_only_a_dropped_set_holds_are_not_in_the_background():
    """fastwilcoxGMT takes as background the genes of the sets it is given; the observed side does not give it the sets over max_g."""
    out = _r('''
cs <- matrix(seq(0.05, 1.2, length.out = 12), 12, 1, dimnames = list(sprintf("g%02d", 1:12), NULL))
gs <- list(M = sprintf("g%02d", 4:9), X = c("g01", "g02", "g03", "g10"), big = sprintf("g%02d", 1:12))     # g11, g12: only in `big`
gmts <- list(d = list(genesets = gs, geneset.names = names(gs)))
real <- list(d = data.frame(pval = rep(NA_real_, 3), stat = NA_real_, row.names = names(gs)))
v <- fcs_null_enrichstat_vectorized(cs, gmts, real, num_g = 5, max_g = 10)$d
# background g01..g10 without `big`: M holds ranks 4..9 against 4 genes (g01..g03, g10): U = 39 - 21 = 18, AUC = 18 / 24 = 0.75.
# With g11 and g12 in the background it would be 18 / 36 = 0.5.
cat(v["M", 1], is.na(v["big", 1]), is.na(v["X", 1]))
''')
    assert out == "0.25 TRUE TRUE"


def _has_rerconverge():
    out = subprocess.run(["Rscript", "-e", 'cat(requireNamespace("RERconverge", quietly = TRUE))'], capture_output=True, text=True)
    return out.stdout.strip() == "TRUE"


def test_the_vectorized_statistic_equals_fastwilcoxGMTall_where_rerconverge_is_installed():
    if not _has_rerconverge():
        pytest.skip("RERconverge is not installed here; run check_fcs_null_equivalence.R in the pipeline environment")
    r = subprocess.run(["Rscript", str(CHECK), str(SCRIPT.parents[4])], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-800:]
    assert r.stdout.count("PASS") == 2 and "FAIL" not in r.stdout
