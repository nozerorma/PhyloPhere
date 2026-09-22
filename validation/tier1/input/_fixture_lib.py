"""Shared helpers for the Tier 1 fixture build scripts (pepc/scripts/build.py,
and its respective build_cds.py).
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            if name is not None:
                out[name] = "".join(buf)
            name = ln[1:].strip()
            buf = []
        else:
            buf.append(ln.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    with path.open("w") as fh:
        for name, seq in sorted(seqs.items()):
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def date_tree(outdir: Path, outgroup: set[str], drop_tips: list[str] = ()) -> None:
    """Root outdir/tree_substitution.nwk (already-written raw phylogram, tip
    labels as used in `outgroup`/`drop_tips`) on `outgroup`, optionally
    pruning `drop_tips` first, then time-scale it to an ultrametric
    chronogram (ape::chronos, penalised likelihood, lambda=1, correlated
    rates, root age=1), writing outdir/tree.nwk.

    PhyloPhere's contrast-independence test (modified Dunn) and its OU/BM
    Phylogenetic Shift Score both assume a TIME tree -- both fixtures hit a
    concrete case of this: a fast-evolving terminal branch (Killinga in
    PEPC) inflates the diameter of any contrast
    pair containing it on the raw phylogram, pushing a real signal below the
    Dunn threshold; dating fixes it. See either fixture's README for detail.
    """
    phylo = outdir / "tree_substitution.nwk"
    drop_r = ", ".join(f'"{t}"' for t in drop_tips)
    og = ", ".join(f'"{s}"' for s in sorted(outgroup))
    rscript = (
        "suppressPackageStartupMessages(library(ape)); "
        f'p <- read.tree("{phylo}"); '
        + (f"p <- drop.tip(p, c({drop_r})); "
           f'write.tree(p, "{phylo}"); ' if drop_tips else "")
        + f"r <- root(p, outgroup = c({og}), resolve.root = TRUE); "
        "r <- multi2di(r); r$edge.length[r$edge.length <= 0] <- 1e-8; "
        'u <- chronos(r, lambda = 1, model = "correlated", quiet = TRUE); '
        'u <- ladderize(structure(unclass(u), class = "phylo")); '
        "stopifnot(is.rooted(u), is.ultrametric(u, tol = 1e-6)); "
        f'write.tree(u, "{outdir / "tree.nwk"}")'
    )
    subprocess.run(["Rscript", "-e", rscript], check=True, capture_output=True, text=True)


def rename_tree_tips(path: Path, name_map: dict[str, str]) -> None:
    """In place: swap every tip label in the newick file at `path` via
    name_map (must cover every tip currently in the file)."""
    txt = path.read_text()
    txt = re.sub(r"([(,])([A-Za-z0-9._-]+?):",
                 lambda m: f"{m.group(1)}{name_map[m.group(2)]}:", txt)
    path.write_text(txt)
