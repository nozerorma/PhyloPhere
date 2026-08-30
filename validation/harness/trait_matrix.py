"""Expand a DatasetSpec into every PhyloPhere trait-mode input file.

Covers the full cross-product of the three trait axes:

    type       : categorical | continuous
    count      : single | multi   (n_trait)
    grouping   : direct | percentilized

File formats (verified against the pipeline source):

* CAAS discovery config  — ``pindex.load_cfg`` : TSV, no header,
  ``species <tab> label(1|0) <tab> pair_id``. n_trait "multi" = a *directory*
  of such files, one per trait (``load_cfg(mode="multi")`` globs the dir).
* Phenotype / traitvalues — ``permulations.R`` : TSV. Either headerless
  ``species <tab> value`` or an annotated table with a ``species`` column and
  one or more numeric columns (``perm_pheno_col`` selects, else first numeric).

Percentilized grouping: a continuous vector -> top-q / bottom-q tails become the
fg / bg groups; the middle is dropped. Pairs are formed by rank-matching
(most-extreme fg with most-extreme bg). This is a stand-in for the
phylogenetically-informed pairing that ``contrast_selection`` does in
production; the manifest records which pairing was used.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import asdict
from pathlib import Path

from .phenotype import DatasetSpec, PhenotypeSpec

DEFAULT_QUANTILE = 0.25


# --------------------------------------------------------------------------- #
# grouping / pairing
# --------------------------------------------------------------------------- #
def percentilize(spec: PhenotypeSpec, q: float) -> tuple[list[str], list[str], list[str]]:
    """Return (foreground, background, dropped_middle) tip lists for a continuous
    trait split at the q / (1-q) quantiles."""
    if spec.kind != "continuous":
        raise TypeError("percentilize() needs a continuous trait")
    if not 0.0 < q < 0.5:
        raise ValueError(f"quantile must be in (0, 0.5), got {q}")
    ranked = sorted(spec.values.items(), key=lambda kv: kv[1])
    n = len(ranked)
    k = max(1, round(n * q))
    low = [t for t, _ in ranked[:k]]
    high = [t for t, _ in ranked[-k:]]
    middle = [t for t, _ in ranked[k:n - k]]
    if spec.higher_is_foreground:
        return high, low, middle
    return low, high, middle


def _rank_pairs(fg: list[str], bg: list[str], values: dict[str, float],
                higher_is_fg: bool) -> list[tuple[str, str]]:
    """Rank-match most-extreme fg with most-extreme bg. Truncates to the shorter
    list (unpaired tips are left out of the config, as paired mode requires)."""
    fg_sorted = sorted(fg, key=lambda t: values[t], reverse=higher_is_fg)
    bg_sorted = sorted(bg, key=lambda t: values[t], reverse=not higher_is_fg)
    return list(zip(fg_sorted, bg_sorted))


def _resolve_pairs(spec: PhenotypeSpec) -> tuple[list[tuple[str, str]], str]:
    """(pairs, method) for a categorical trait."""
    if spec.pairs:
        return spec.pairs, "explicit"
    fg, bg = spec.foreground(), spec.background()
    # order-pair deterministically; no phylo info available
    return list(zip(sorted(fg), sorted(bg))), "order(no-phylo)"


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def _write_cfg(path: Path, pairs: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for i, (fg_tip, bg_tip) in enumerate(pairs, start=1):
        lines.append(f"{fg_tip}\t1\t{i}")
        lines.append(f"{bg_tip}\t0\t{i}")
    path.write_text("\n".join(lines) + "\n")


def _write_traitvalues_single(path: Path, values: dict[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [f"{t}\t{v:g}" for t, v in sorted(values.items())]
    path.write_text("\n".join(rows) + "\n")


def _write_traitvalues_annotated(path: Path, columns: dict[str, dict[str, float]],
                                 extra: dict[str, dict[str, float]] | None = None) -> None:
    """Annotated TSV with a header: species + one column per trait (+ optional
    extra columns such as percentile ranks)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    extra = extra or {}
    tips = sorted({t for col in columns.values() for t in col})
    head = ["species", *columns.keys(), *extra.keys()]
    rows = ["\t".join(head)]
    for t in tips:
        cells = [t]
        for col in columns.values():
            cells.append(f"{col.get(t, float('nan')):g}")
        for col in extra.values():
            cells.append(f"{col.get(t, float('nan')):g}")
        rows.append("\t".join(cells))
    path.write_text("\n".join(rows) + "\n")


def _percentile_ranks(values: dict[str, float]) -> dict[str, float]:
    """percent_rank-style transform, matching phen_score: rank / (n-1) in [0,1]."""
    ranked = sorted(values, key=lambda t: values[t])
    n = len(ranked)
    if n == 1:
        return {ranked[0]: 0.0}
    return {t: i / (n - 1) for i, t in enumerate(ranked)}


# --------------------------------------------------------------------------- #
# tree pruning
# --------------------------------------------------------------------------- #
def _prune_tree(src: Path, dst: Path, keep: set[str]) -> str:
    """Prune ``src`` to the tips in ``keep`` (preserving branch lengths) and write
    to ``dst``. Prefers dendropy (what the pipeline uses); falls back to Bio.Phylo,
    then to an unpruned copy with a loud manifest warning."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    keep = set(keep)

    try:
        import dendropy  # type: ignore

        t = dendropy.Tree.get(path=str(src), schema="newick",
                              preserve_underscores=True)
        present = {l.taxon.label for l in t.leaf_node_iter() if l.taxon}
        t.retain_taxa_with_labels(sorted(keep & present))
        t.suppress_unifurcations()
        t.write(path=str(dst), schema="newick", suppress_rooting=True,
                unquoted_underscores=True)
        missing = sorted(keep - present)
        return "dendropy" + (f"; MISSING_FROM_TREE={missing}" if missing else "")
    except ImportError:
        pass
    except Exception as exc:  # parse/prune failure — do not guess, copy + flag
        shutil.copyfile(src, dst)
        return f"copied-unpruned (dendropy: {type(exc).__name__}: {exc})"

    try:
        from Bio import Phylo  # type: ignore

        t = Phylo.read(str(src), "newick")
        present = {tip.name for tip in t.get_terminals()}
        for name in present - keep:
            t.prune(name)
        Phylo.write(t, str(dst), "newick")
        missing = sorted(keep - present)
        return "Bio.Phylo" + (f"; MISSING_FROM_TREE={missing}" if missing else "")
    except Exception as exc:
        shutil.copyfile(src, dst)
        return f"copied-unpruned ({type(exc).__name__}: {exc})"


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #
def emit_trait_matrix(dataset: DatasetSpec, outdir: str | Path,
                      quantile: float = DEFAULT_QUANTILE) -> dict:
    """Write the full trait matrix for ``dataset`` under ``outdir``. Returns the
    manifest dict (also written to ``outdir/manifest.json``)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cat_traits = [t for t in dataset.traits if t.kind == "categorical"]
    con_traits = [t for t in dataset.traits if t.kind == "continuous"]
    manifest: dict = {
        "dataset": dataset.name,
        "provenance": dataset.provenance,
        "quantile": quantile,
        "emitted": {},
        "pairing": {},
        "warnings": [],
    }

    # -- tree ------------------------------------------------------------------
    tree_note = _prune_tree(dataset.tree, outdir / "tree" / "tree.nwk", dataset.tips)
    manifest["tree"] = {"source": str(dataset.tree), "note": tree_note}
    if "MISSING" in tree_note or "copied-unpruned" in tree_note:
        manifest["warnings"].append(f"tree: {tree_note}")

    # -- categorical: direct -------------------------------------------------- #
    for t in cat_traits:
        pairs, method = _resolve_pairs(t)
        manifest["pairing"][t.name] = {"method": method, "n_pairs": len(pairs)}
        _write_cfg(outdir / "categorical/direct/single" / f"{t.name}.cfg", pairs)
    if len(cat_traits) > 1:
        for t in cat_traits:
            pairs, _ = _resolve_pairs(t)
            _write_cfg(outdir / "categorical/direct/multi" / t.name, pairs)
    manifest["emitted"]["categorical/direct"] = [t.name for t in cat_traits]

    # -- categorical: percentilized (from continuous) ----------------------- #
    perc_names = []
    for t in con_traits:
        fg, bg, middle = percentilize(t, quantile)
        pairs = _rank_pairs(fg, bg, t.values, t.higher_is_foreground)
        tag = f"{t.name}.q{int(quantile * 100)}"
        perc_names.append(tag)
        manifest["pairing"][tag] = {
            "method": "rank-match(no-phylo)",
            "n_pairs": len(pairs),
            "dropped_middle": len(middle),
        }
        _write_cfg(outdir / "categorical/percentilized/single" / f"{tag}.cfg", pairs)
    if len(con_traits) > 1:
        for t in con_traits:
            fg, bg, _ = percentilize(t, quantile)
            pairs = _rank_pairs(fg, bg, t.values, t.higher_is_foreground)
            _write_cfg(outdir / "categorical/percentilized/multi" / f"{t.name}.q{int(quantile*100)}", pairs)
    manifest["emitted"]["categorical/percentilized"] = perc_names

    # -- continuous: direct ------------------------------------------------- #
    for t in con_traits:
        _write_traitvalues_single(
            outdir / "continuous/direct/single" / f"{t.name}.traitvalues.tsv", t.values
        )
    if len(con_traits) > 1:
        _write_traitvalues_annotated(
            outdir / "continuous/direct/multi" / "traitvalues.tsv",
            {t.name: t.values for t in con_traits},
        )
    manifest["emitted"]["continuous/direct"] = [t.name for t in con_traits]

    # -- continuous: percentilized ---------------------------------------- #
    for t in con_traits:
        pr = _percentile_ranks(t.values)
        _write_traitvalues_annotated(
            outdir / "continuous/percentilized/single" / f"{t.name}.traitvalues.tsv",
            {t.name: t.values},
            extra={f"{t.name}__pctile": pr},
        )
    if len(con_traits) > 1:
        _write_traitvalues_annotated(
            outdir / "continuous/percentilized/multi" / "traitvalues.tsv",
            {t.name: t.values for t in con_traits},
            extra={f"{t.name}__pctile": _percentile_ranks(t.values) for t in con_traits},
        )
    manifest["emitted"]["continuous/percentilized"] = [t.name for t in con_traits]

    # -- provenance hash --------------------------------------------------- #
    spec_bytes = json.dumps(
        {"traits": [asdict(t) for t in dataset.traits]}, sort_keys=True
    ).encode()
    manifest["spec_sha256"] = hashlib.sha256(spec_bytes).hexdigest()

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
