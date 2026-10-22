# Run from the repository root: Rscript tests/test_rer_permulations.R <marine_complete>
# Fixed-seed comparisons use the installed pinned RERconverge as the reference.
suppressPackageStartupMessages(library(RERconverge))
repo <- normalizePath(".")
source(file.path(repo, "subworkflows/RERCONVERGE/local/rer_permulations.R"))
run <- normalizePath(commandArgs(TRUE)[1])
trees <- readRDS(file.path(run, "rerconverge/rer_objects/gene_trees.nwk.masterTree.output"))
RER <- readRDS(file.path(run, "rerconverge/rer_objects/marine.RERmatrix.output"))[1:8, , drop = FALSE]
observed <- readRDS(file.path(run, "rerconverge/rer_results/marine.binary.output"))
fg <- attr(observed, "fg_sp")
stopifnot(length(fg) >= 3)
quiet <- function(expr) invisible(capture.output(suppressMessages(expr)))
same <- function(a, b) {
  for (n in c("corP", "corRho", "corStat")) {
    aa <- as.matrix(a[[n]]); bb <- as.matrix(b[[n]])
    bb <- bb[rownames(aa), , drop = FALSE]
    stopifnot(identical(unname(is.na(aa)), unname(is.na(bb))),
              isTRUE(all.equal(unname(aa), unname(bb), tolerance = 1e-12)))
  }
}

set.seed(20261022)
quiet(native <- getPermsBinary(3, fg, NA, NA, RER, trees, trees$masterTree, calculateenrich = FALSE))
set.seed(20261022)
quiet(adapted <- rer_perms_binary(3, fg, RER, trees, min.sp = 10))
same(native, adapted)
cat("PASS: binary CC default settings reproduce native null, fixed seed, 3 cycles\n")

# Independent scoring of the native categorical histories with non-default settings.
phen <- setNames(as.numeric(trees$masterTree$tip.label %in% fg), trees$masterTree$tip.label)
set.seed(71)
quiet(histories <- categoricalPermulations(trees, phen, rm = "ER", rp = "auto", ntrees = 2))
for (clade in c("all", "terminal", "ancestral")) {
  reference <- lapply(histories$trees, function(x) {
    tr <- trees$masterTree
    tr$edge.length <- c(x$tips, x$nodes)[tr$edge[, 2]] - 1
    names(tr$edge.length) <- NULL
    fg_j <- getFromNamespace("getForegroundsFromBinaryTree", "RERconverge")(tr)
    paths <- foreground2Paths(fg_j, trees, clade = clade)
    suppressMessages(correlateWithBinaryPhenotype(RER, paths, min.sp = 40, min.pos = 4,
                                                 winsorizeRER = 3, weighted = "auto"))
  })
  ref <- list(corP = do.call(cbind, lapply(reference, `[[`, "P")),
              corRho = do.call(cbind, lapply(reference, `[[`, "Rho")))
  ref$corStat <- sign(ref$corRho) * -log10(ref$corP)
  ref <- lapply(ref, function(x) { rownames(x) <- rownames(RER); x })
  set.seed(71)
  quiet(got <- rer_perms_binary(2, fg, RER, trees, min.sp = 40, min.pos = 4, winR = 3, clade = clade))
  same(got, ref)
}
set.seed(9)
quiet(empty <- rer_perms_binary(1, fg, RER, trees, min.sp = 10000, winR = 3))
stopifnot(all(is.na(as.matrix(empty$corStat))))
cat("PASS: binary clade all/terminal/ancestral, min.sp=40, min.pos=4, winsor=3; high cutoff keeps NA; one-cycle batch\n")

trait <- setNames(seq_along(trees$masterTree$tip.label), trees$masterTree$tip.label)
mt <- trees$masterTree
if (!ape::is.rooted(mt)) mt <- ape::root(mt, outgroup = mt$tip.label[1], resolve.root = TRUE)
mt <- ape::multi2di(mt)
mt$edge.length[mt$edge.length <= 0] <- 1e-8
set.seed(84)
quiet(native <- getPermsContinuous(2, trait, RER, NULL, trees, mt, calculateenrich = FALSE))
set.seed(84)
quiet(adapted <- rer_perms_continuous(2, trait, RER, trees, mt, min.sp = 10))
same(native, adapted)
cat("PASS: continuous simperm default settings reproduce native null, fixed seed, 2 cycles\n")

set.seed(12)
vec <- simpermvec(trait, mt)
paths <- char2Paths(stats::na.omit(vec), trees)
quiet(expected <- correlateWithContinuousPhenotype(RER, paths, min.sp = 40, winsorizeRER = 1, winsorizetrait = 2))
set.seed(12)
quiet(got <- rer_perms_continuous(1, trait, RER, trees, mt, min.sp = 40, winR = 1, winT = 2))
stopifnot(isTRUE(all.equal(got$corP[, 1], expected$P, tolerance = 1e-12)),
          isTRUE(all.equal(got$corRho[, 1], expected$Rho, tolerance = 1e-12)))
set.seed(12)
quiet(empty <- rer_perms_continuous(1, trait, RER, trees, mt, min.sp = 10000))
stopifnot(all(is.na(as.matrix(empty$corStat))), is.null(rer_winsor_limit(0)))
quiet(no_win <- rer_score_paths(RER, paths, "continuous", 10, winR = 0, winT = 0))
quiet(no_win_ref <- correlateWithContinuousPhenotype(RER, paths, min.sp = 10, winsorizeRER = NULL, winsorizetrait = NULL))
stopifnot(identical(no_win$P, no_win_ref$P))
cat("PASS: continuous min.sp and both winsorization settings match direct native scoring; cutoff and zero-winsor policies\n")
