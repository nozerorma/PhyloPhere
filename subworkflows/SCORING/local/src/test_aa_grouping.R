#!/usr/bin/env Rscript
# Tests for aa_grouping.R — run: Rscript test_aa_grouping.R
# Asserts every one of the 20 canonical amino acids maps to the label documented
# in subworkflows/CT_DISAMBIGUATION/local/src/biochem/grouping.py, for all 5
# schemes (100 assertions), plus unknown-residue pass-through and case handling.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "aa_grouping.R"))

ok <- TRUE
check <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) ok <<- FALSE
}

# Expected maps transcribed from grouping.py (source of truth).
expected <- list(
  US = c(A="A",C="C",D="D",E="E",F="F",G="G",H="H",I="I",K="K",L="L",M="M",
         N="N",P="P",Q="Q",R="R",S="S",T="T",V="V",W="W",Y="Y"),
  GS1 = c(C="t",V="t",A="n",G="n",P="n",S="n",N="p",D="p",Q="p",E="p",
          R="b",H="b",K="b",I="h",L="h",M="h",F="h",W="h",Y="h",T="o"),
  GS2 = c(C="c",A="s",G="s",V="s",D="a",E="a",N="n",Q="n",H="n",W="n",
          R="b",K="b",I="h",L="h",F="h",P="h",Y="x",M="x",T="x",S="x"),
  GS3 = c(C="c",A="n",G="n",P="n",S="n",T="n",N="s",D="s",Q="s",E="s",
          R="b",H="b",K="b",I="l",L="l",M="l",V="l",F="g",W="g",Y="g"),
  GS4 = c(C="c",A="h",I="h",L="h",V="h",S="o",T="o",N="p",Q="p",D="a",E="a",
          R="b",H="b",G="g",P="r",K="k",M="m",F="f",W="y",Y="y")
)

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
for (sc in names(expected)) {
  for (aa in aas) {
    got <- encode_aa_r(aa, sc)
    want <- unname(expected[[sc]][[aa]])
    check(identical(got, want), sprintf("%s: %s -> %s (got %s)", sc, aa, want, got))
  }
}

# Group-count sanity (US 20, GS1 6, GS2 7, GS3 6, GS4 12).
grp_n <- function(sc) length(unique(vapply(aas, encode_aa_r, character(1), scheme = sc)))
check(grp_n("US") == 20, "US has 20 groups")
check(grp_n("GS1") == 6,  "GS1 has 6 groups")
check(grp_n("GS2") == 7,  "GS2 has 7 groups")
check(grp_n("GS3") == 6,  "GS3 has 6 groups")
check(grp_n("GS4") == 12, "GS4 has 12 groups")

# Unknown residue passes through raw (mirrors encode_aa); case-insensitive.
check(identical(encode_aa_r("Z", "GS4"), "Z"),  "unknown residue -> raw")
check(identical(encode_aa_r("x", "GS1"), "X"),  "unknown residue upper-cased")
check(identical(encode_aa_r("v", "GS3"), "l"),  "lower-case input encodes")
check(is.na(encode_aa_r("", "GS4")),            "empty -> NA")
check(is.na(encode_aa_r(NA, "US")),             "NA -> NA")
check(identical(encode_aa_r("V", "bogus"), "V"),"unknown scheme -> raw residue")

cat("\n", if (ok) "ALL TESTS PASSED" else "SOME TESTS FAILED", "\n", sep = "")
quit(status = if (ok) 0 else 1)
