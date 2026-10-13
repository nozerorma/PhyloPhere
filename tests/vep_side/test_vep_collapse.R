# test_vep_collapse.R: collapse_vep_consequence(): one consequence per (Gene, Position, side).
# The function is read from the marked block of the report, so the test exercises the code the report runs.
# Run: Rscript tests/vep_side/test_vep_collapse.R
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", ".."))
rmd  <- readLines(file.path(root, "subworkflows/ENRICHMENT/local/14.Position_enrichment_report.Rmd"))
i <- grep("^# <vep_collapse>", rmd); j <- grep("^# </vep_collapse>", rmd)
stopifnot(length(i) == 1, length(j) == 1)
eval(parse(text = rmd[i:j]))
# The copy in report 15 must be identical.
rmd15 <- readLines(file.path(root, "subworkflows/ENRICHMENT/local/15.Comparison_report.Rmd"))
stopifnot(identical(rmd[i:j], rmd15[grep("^# <vep_collapse>", rmd15):grep("^# </vep_collapse>", rmd15)]))

mk <- function(...) data.frame(..., stringsAsFactors = FALSE)
d <- mk(
  Gene = "G", Position = 1L, side = "top",
  Feature = c("T1", "T2", "T3", "T4"),
  Consequence = c("start_lost", "missense_variant", "missense_variant,splice_region_variant", "intron_variant"),
  IMPACT = c("HIGH", "MODERATE", "MODERATE", "MODIFIER"),
  Extra = c("IMPACT=HIGH", "IMPACT=MODERATE;CANONICAL=YES", "IMPACT=MODERATE;MANE_SELECT=NM_1", "IMPACT=MODIFIER")
)
r <- collapse_vep_consequence(d)
# MANE Select wins over the canonical row and over the more severe unflagged row.
stopifnot(nrow(r) == 1, r$Feature == "T3", r$transcript_rule == "MANE Select",
          r$top_consequence == "missense_variant", !r$missense_only,
          r$Consequence == "missense_variant,splice_region_variant")

# Canonical when no MANE row exists.
r <- collapse_vep_consequence(d[c(1, 2, 4), ])
stopifnot(r$Feature == "T2", r$transcript_rule == "canonical", r$missense_only)

# No flagged row: most severe consequence over all transcripts.
d3 <- d[c(1, 4), ]; d3$Extra <- c("IMPACT=HIGH", "IMPACT=MODIFIER")
r <- collapse_vep_consequence(d3)
stopifnot(r$Feature == "T1", r$transcript_rule == "most severe", r$top_consequence == "start_lost")

# Flags as columns of their own (vep 116 tab output): MANE_SELECT / MANE / CANONICAL, no Extra column.
dc <- mk(Gene = "G", Position = 1L, side = "top", Feature = c("T1", "T2", "T3"),
         Consequence = c("start_lost", "missense_variant", "missense_variant"),
         CANONICAL = c("-", "YES", "-"), MANE = c("-", "-", "MANE_Select"), MANE_SELECT = c("-", "-", "NM_000546.6"))
r <- collapse_vep_consequence(dc)
stopifnot(r$Feature == "T3", r$transcript_rule == "MANE Select")
r <- collapse_vep_consequence(dc[1:2, ])
stopifnot(r$Feature == "T2", r$transcript_rule == "canonical")

# Groups stay separate by side, and by Gene and Position; a table without Extra and Feature still works.
d4 <- rbind(
  mk(Gene = "G", Position = 1L, side = "top",    Consequence = "missense_variant"),
  mk(Gene = "G", Position = 1L, side = "bottom", Consequence = "synonymous_variant,stop_retained_variant"),
  mk(Gene = "H", Position = 1L, side = "top",    Consequence = "-"))
r <- collapse_vep_consequence(d4)
stopifnot(nrow(r) == 2, r$missense_only[r$side == "top"], !r$missense_only[r$side == "bottom"],
          r$top_consequence[r$side == "bottom"] == "stop_retained_variant")   # stop_retained ranks above synonymous

# Fallback grouping without side.
r <- collapse_vep_consequence(d4[1:2, ], by = c("Gene", "Position"))
stopifnot(nrow(r) == 1)

# Zero-row inputs keep the schema.
r <- collapse_vep_consequence(d4[3, ])
stopifnot(nrow(r) == 0, all(c("Gene", "Position", "side", "Consequence", "missense_only") %in% names(r)))
r <- collapse_vep_consequence(d4[0, ])
stopifnot(nrow(r) == 0)
cat("test_vep_collapse: OK\n")
