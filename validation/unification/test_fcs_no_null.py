"""FCS for a ranking with no permulation null: no private-shuffle path-sum, and no evidence label that claims a phylogenetic gate."""
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "subworkflows/ENRICHMENT/local/src/fcs_enrich.R"

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")

PRELUDE = f'''
suppressPackageStartupMessages(library(dplyr))
for (e in parse("{SCRIPT}")) {{
  if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) == "fcs_classify_evidence") eval(e)
}}
'''

# three pathways per ranking: both gates, one gate, none; the permulation gate would pass for the first two
FRAME = '''
one <- function(rk) data.frame(ranking = rk, database = "d", pathway = c("both", "one", "none"),
                 p.adj = c(0.01, 0.01, 0.9), p.perm = NA_real_, stat = 1,
                 lach_p.adj = c(0.01, 0.9, 0.9), lach_p.perm = NA_real_,
                 perm_p.adj = c(0.001, 0.001, 0.9), perm_nes = 2, stringsAsFactors = FALSE)
df <- rbind(one("global"), one("top"))
cl <- function(...) fcs_classify_evidence(df, fdr_wilcoxon = 0.15, fdr_lachenbruch = 0.15, fdr_permsum = 0.15, p_perm_thr = 0.025, ...)
'''


def _r(code):
    out = subprocess.run(["Rscript", "-e", PRELUDE + FRAME + code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-1200:]
    return out.stdout.strip().splitlines()


def test_with_a_null_for_every_ranking_the_labels_are_the_known_ones():
    labels, perm = _r("r <- cl(); cat(r$evidence_label, sep = '|'); cat('\\n'); cat(r$sig_permulation, sep = '|')")
    assert labels == "|".join(["Hard evidence", "Supported", "Not significant"] * 2)
    assert perm == "TRUE|TRUE|FALSE|TRUE|TRUE|FALSE"


def test_a_ranking_without_a_null_has_no_phylogenetic_label_and_its_permulation_flag_cannot_pass():
    labels, perm, counts = _r("r <- cl(no_null_rankings = c('global', 'top')); cat(r$evidence_label, sep = '|'); cat('\\n'); "
                              "cat(r$sig_permulation, sep = '|'); cat('\\n'); cat(r$evidence_count, sep = '|')")
    assert labels == "|".join(["Exploratory (relative)", "Exploratory (relative)", "Not significant"] * 2)
    assert perm == "|".join(["FALSE"] * 6)
    assert counts == "2|1|0|2|1|0"       # the gates that could be applied still count


def test_only_the_rankings_without_a_null_are_capped():
    labels, = _r("r <- cl(no_null_rankings = 'global'); cat(r$evidence_label, sep = '|')")
    assert labels == "|".join(["Exploratory (relative)", "Exploratory (relative)", "Not significant",
                               "Hard evidence", "Supported", "Not significant"])


def test_the_columns_of_the_input_are_kept_and_no_helper_column_leaks():
    names, = _r("r <- cl(no_null_rankings = 'global'); cat(setdiff(names(r), names(df)), sep = ',')")
    assert names == "sig_wilcoxon,sig_lachenbruch,sig_permulation,evidence_count,evidence_label"


def test_fcs_run_all_skips_the_private_shuffle_for_a_ranking_without_a_null():
    text = SCRIPT.read_text()
    assert "private-shuffle" not in text and "fcs_null_is_empty" not in text
    assert "if (is.null(corStat_rk)) {" in text and "skipped (no permulation null for this ranking)" in text
    assert "no_null_rankings = names(corStat_byrk)[vapply(corStat_byrk, is.null, logical(1))]" in text
