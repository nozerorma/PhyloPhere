"""FCS with a null that holds no cycle (N = 0): no private-shuffle path-sum, and no evidence label that claims a phylogenetic gate."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "subworkflows/ENRICHMENT/local/src/fcs_enrich.R"

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")

PRELUDE = f'''
suppressPackageStartupMessages(library(dplyr))
for (e in parse("{SCRIPT}")) {{
  if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("fcs_null_is_empty", "fcs_classify_evidence")) eval(e)
}}
'''


def _r(code):
    out = subprocess.run(["Rscript", "-e", PRELUDE + code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-1200:]
    return out.stdout.strip()


@pytest.mark.parametrize("perms,expected", [
    ("list(corStat_byrank = list(), caas_corStat_byrank = list(), gene_stat = 'size_adj_max')", "TRUE"),        # what scoring_caas_perms.R writes for N = 0
    ("list(caas_corStat_byrank = list(global = matrix(0, 2, 3)))", "FALSE"),
    ("list(corStat = matrix(0, 2, 3))", "FALSE"),
    ("list(gene_stat = 'size_adj_max')", "FALSE"),                                                              # not a null at all
    ("NULL", "FALSE"),
])
def test_only_a_null_with_empty_matrix_lists_is_empty(perms, expected):
    assert _r(f"cat(fcs_null_is_empty({perms}))") == expected


FRAME = '''
df <- data.frame(ranking = "global", database = "d", pathway = c("both", "one", "none"),
                 p.adj = c(0.01, 0.01, 0.9), p.perm = NA_real_, stat = 1,
                 lach_p.adj = c(0.01, 0.9, 0.9), lach_p.perm = NA_real_,
                 perm_p.adj = c(0.001, 0.001, 0.9), perm_nes = 2, stringsAsFactors = FALSE)
cl <- function(...) fcs_classify_evidence(df, fdr_wilcoxon = 0.15, fdr_lachenbruch = 0.15, fdr_permsum = 0.15, p_perm_thr = 0.025, ...)
'''


def test_with_a_usable_null_the_labels_are_the_known_ones():
    out = _r(FRAME + "r <- cl(null_empty = FALSE); cat(r$evidence_label, sep = '|')")
    assert out == "Hard evidence|Supported|Not significant"


def test_with_an_empty_null_no_row_is_called_phylogenetic_and_the_permulation_flag_cannot_pass():
    out = _r(FRAME + "r <- cl(null_empty = TRUE); cat(r$evidence_label, sep = '|'); cat('\\n'); cat(r$sig_permulation, sep = '|'); cat('\\n'); cat(r$evidence_count, sep = '|')")
    labels, perm, counts = out.splitlines()
    assert labels == "Exploratory (relative)|Exploratory (relative)|Not significant"
    assert perm == "FALSE|FALSE|FALSE"
    assert counts == "2|1|0"       # the gates that could be applied still count


def test_the_empty_null_is_detected_in_fcs_run_all_and_skips_the_private_shuffle():
    text = SCRIPT.read_text()
    assert "null_empty <- fcs_null_is_empty(corperms)" in text
    assert 'ranking %s skipped (the null holds no cycle)' in text
    assert "fcs_classify_evidence(enrich_df" in text and "null_empty = null_empty" in text
