# Run from the repository root with the project's Rscript.
repo <- normalizePath(".")
source(file.path(repo, "subworkflows/ENRICHMENT/local/src/term_blocks.R"))

# Complete linkage: {a,b} and {c,d} both have Jaccard 0.5 with {a,b,c,d} but share nothing,
# so they never join the same block at thr 0.5 (a greedy representative rule would merge them).
s <- list(big = c("a", "b", "c", "d"), left = c("a", "b"), right = c("c", "d"))
r <- term_blocks(s, rank = 1:3, thr = 0.5)
stopifnot(r$block[r$key == "left"] != r$block[r$key == "right"],
          sum(r$block == r$block[1]) == 2)

# Every block cut by height has a minimum pairwise Jaccard >= thr.
set.seed(1)
genes <- paste0("g", 1:60)
sets <- setNames(lapply(1:40, function(i) sample(genes, sample(5:20, 1))), paste0("t", 1:40))
for (thr in c(.1, .2, .4)) {
  b <- term_blocks(sets, rank = runif(40), thr = thr)
  stopifnot(all(b$block_min_jaccard[!is.na(b$block_min_jaccard)] >= thr - 1e-9))
}

# Lowering the threshold merges blocks (never splits them), and blocks are numbered by best member.
lo <- term_blocks(sets, rank = seq_along(sets), thr = .1)
hi <- term_blocks(sets, rank = seq_along(sets), thr = .4)
stopifnot(length(unique(lo$block)) <= length(unique(hi$block)),
          identical(hi$block[hi$is_representative], sort(hi$block[hi$is_representative])))

# The representative is the best-ranked member; ties break by name, never by row order.
rk <- c(t1 = 2, t2 = 1, t3 = 1)
same_sets <- list(t1 = c("a", "b"), t2 = c("a", "b"), t3 = c("a", "b"))
r1 <- term_blocks(same_sets, rank = rk, thr = .5)
r2 <- term_blocks(rev(same_sets), rank = rev(rk), thr = .5)
stopifnot(identical(r1$key[r1$is_representative], "t2"),
          identical(r2$key[r2$is_representative], "t2"),
          all(r1$block == 1L))

# Block 1 holds the best term; terms without genes are their own block; a lone term works.
e <- term_blocks(list(x = character(0), y = c("a", "b"), z = c("a", "b")), rank = c(3, 1, 2), thr = .2)
stopifnot(e$block[e$key == "y"] == 1L, e$block[e$key == "z"] == 1L, e$block[e$key == "x"] == 2L,
          is.na(e$block_min_jaccard[e$key == "x"]))
one <- term_blocks(list(a = c("g1")), rank = 1)
stopifnot(nrow(one) == 1L, one$is_representative, one$block == 1L, is.na(one$block_min_jaccard))
stopifnot(nrow(term_blocks(list())) == 0L)

# fcs_add_clusters: groups per ranking and effect sign over the selected rows, with the full
# gene set (gene.vals), and leaves the other rows NA.
enrich <- data.frame(
  ranking = "global", database = "db", pathway = paste0("P", 1:5),
  stat = c(.3, .25, .2, .2, -.3), q.perm = c(.01, .02, .03, .04, .05), p.perm = .001,
  gene.vals = c("a:1,b:2,c:3", "a:1,b:2,c:3,d:9", "x:1,y:2", "a:1,b:2,c:3", "a:1,b:2,c:3"),
  sel = c(TRUE, TRUE, TRUE, FALSE, TRUE), stringsAsFactors = FALSE)
out <- fcs_add_clusters(enrich, threshold = .5, selection_col = "sel")
stopifnot(identical(out$cluster_id[1:2], rep("global_positive_C001", 2)),
          out$cluster_id[3] == "global_positive_C002",
          is.na(out$cluster_id[4]), out$cluster_id[5] == "global_negative_C001",
          identical(out$is_representative[1:3], c(TRUE, FALSE, TRUE)),
          out$cluster_size[1] == 2L, abs(out$cluster_min_jaccard[1] - 0.75) < 1e-12)
cat("term_blocks tests passed\n")
