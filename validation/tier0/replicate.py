"""Build one Tier 0 replicate: a genome-wide phenotype + many simulated genes.

A replicate is one (tree, phenotype) draw. Within it, ``n_genes`` genes are
simulated on the SAME tree and phenotype; a fraction carry planted convergent
substitutions (the power genes), the rest are pure null. This gives:

  * a population for the gene/position ranking SCORING produces
  * a meaningful ``percent_rank(pvalue_boot)`` (needs many positions)
  * many null p-values for the KS-uniformity gate

Output (pipeline-shaped, consumed by the rendered ``run_single.sh``):

    align/g0001.fasta ... gNNNN.fasta   amino-acid alignments, tips only
    my_traits.tsv                               species <tab> <traitname>   (header)
    tree.nwk                                    the (pruned/scaled) tree
    phenotype.tsv                               species <tab> operative-fg 1|0
    truth.json                                  phenotype + per-gene planted map
                                                + quantiles + every parameter
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from . import pheno as _pheno
from .simulate import SimConfig, simulate
from .trees import PhyloTree


@dataclass
class ReplicateConfig:
    archetype: str = "bodysize"          # "echo" | "bodysize"
    n_genes: int = 60
    frac_planted_genes: float = 0.25
    traitname: str = "sim_trait"
    # per-gene sequence model (SimConfig knobs)
    n_sites: int = 400
    concentration: float = 4.0
    gamma_alpha: float = 0.7
    matrix: str = "wag"
    # planted-gene signal — identical-AA convergence only (D-T0-02). profile_shift
    # is stripped: it is neither a strict-CAAS (US) positive nor a clean grouped
    # CAAP positive (it does not target a shared amino_encoded group), so it only
    # muddied the recall denominator. A proper grouped-CAAP planted mechanism is
    # deferred (D-T0-N).
    n_planted_identical_aa: int = 24
    identical_aa_target_weight: float = 0.95
    # phenotype
    n_transitions: int = 4               # echo: number of independent origins
    bodysize_quantile: float = 0.15      # bodysize: fg = top q of the BM draw
    planted: bool = True                 # False -> null replicate (no planted signal)


def _write_fasta(path: Path, aln: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, seq in aln.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def _write_newick(tree: PhyloTree, path: Path) -> None:
    children: dict[int, list[int]] = {}
    for p, c, _ in tree.edges:
        children.setdefault(p, []).append(c)
    bl_of = {c: bl for _, c, bl in tree.edges}

    def render(node: int) -> str:
        kids = children.get(node, [])
        if not kids:
            return f"{tree.labels[node]}:{bl_of.get(node, 0.0):.8f}"
        inner = ",".join(render(k) for k in kids)
        return (f"({inner});" if node == tree.root
                else f"({inner}){tree.labels[node]}:{bl_of.get(node, 0.0):.8f}")

    path.write_text(render(tree.root) + "\n")


def build_replicate(tree: PhyloTree, rcfg: ReplicateConfig, seed: int,
                    outdir: str | Path) -> dict:
    outdir = Path(outdir)
    (outdir / "align").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    # ── phenotype ───────────────────────────────────────────────────────────
    if not rcfg.planted:
        ph = _pheno.make_null(tree, rng, rcfg.archetype,
                              n_transitions=rcfg.n_transitions,
                              quantile=rcfg.bodysize_quantile)
    elif rcfg.archetype == "echo":
        ph = _pheno.make_echo(tree, rcfg.n_transitions, rng)
    elif rcfg.archetype == "bodysize":
        ph = _pheno.make_bodysize(tree, rng, quantile=rcfg.bodysize_quantile)
    else:
        raise ValueError(f"unknown archetype {rcfg.archetype!r}")

    fg = ph.as_foreground()

    # ── genes ───────────────────────────────────────────────────────────────
    n_planted = 0 if not rcfg.planted else max(1, round(rcfg.n_genes * rcfg.frac_planted_genes))
    planted_flags = np.zeros(rcfg.n_genes, dtype=bool)
    planted_flags[:n_planted] = True
    rng.shuffle(planted_flags)

    genes: dict = {}
    ensembl_rows = ["gene\tchr\tstart\tend\tstrand\tlength\thuman_protein_id"]
    for gi in range(rcfg.n_genes):
        # NO underscore in the gene id: disambiguation_main.py derives the gene
        # name as GenePos.split("_")[0], so "gene_0002_18" -> "gene" and the gene
        # is skipped as "missing alignment". Real HGNC symbols rarely contain "_"
        # so the bug is latent; flagged for Miguel (the `Gene` column is right
        # there in the metadata and should be used instead).
        gid = f"g{gi + 1:04d}"
        # synthetic genomic coordinates: genes spread ~1 Mb apart over 3 fake
        # chromosomes. No deliberate physical clustering of planted genes yet
        # (that would be a separate positional-enrichment validation).
        chrom = f"chr{gi % 3 + 1}"
        start = (gi // 3) * 1_000_000 + 100_000
        length_bp = rcfg.n_sites * 3
        ensembl_rows.append(
            f"{gid}\t{chrom}\t{start}\t{start + length_bp}\t+\t{rcfg.n_sites}\tSYN{gi + 1:05d}"
        )
        is_planted = bool(planted_flags[gi])
        scfg = SimConfig(
            n_sites=rcfg.n_sites, concentration=rcfg.concentration,
            gamma_alpha=rcfg.gamma_alpha, matrix=rcfg.matrix,
            n_planted_profile_shift=0,
            n_planted_identical_aa=rcfg.n_planted_identical_aa if is_planted else 0,
            identical_aa_target_weight=rcfg.identical_aa_target_weight,
        )
        res = simulate(tree, fg if is_planted else None, scfg, seed=seed * 100_000 + gi)
        _write_fasta(outdir / "align" / f"{gid}.fasta", res.alignment)
        genes[gid] = {
            "planted": is_planted,
            "planted_sites": {str(k): v for k, v in res.truth.planted_sites.items()},
            "identical_aa_targets": res.truth.identical_aa_targets,
            "n_planted": len(res.truth.planted_sites),
        }

    # ── phenotype / trait files ─────────────────────────────────────────────
    _write_newick(tree, outdir / "tree.nwk")
    # species list — NAME_CURATION scans every alignment to derive this if absent
    (outdir / "ali_sp_names.txt").write_text("\n".join(sorted(tree.tips)) + "\n")
    (outdir / "gene_ensembl.tsv").write_text("\n".join(ensembl_rows) + "\n")

    sp_sorted = sorted(ph.values)
    # `family` — the reporting / CI Rmds need a `taxon_of_interest` column; the
    # taxid map needs a `family` column too. Synthetic: 5 fake families.
    fam_of = {sp: f"simfam{i % 5 + 1}" for i, sp in enumerate(sp_sorted)}

    # taxid mapping — disambiguation renames tips to numeric ids before writing
    # the PAML alignment (reconstruct.py::_write_phylip uses a single-space name
    # delimiter that breaks on species names containing letters like 'O';
    # production always supplies --tax_id so this path is what real runs use).
    tax_rows = ["tax_id\tspecies\tfamily\trank\tname_class"]
    for i, sp in enumerate(sp_sorted):
        tax_rows.append(f"{9_000_001 + i}\t{sp}\t{fam_of[sp]}\tspecies\tscientific name")
    (outdir / "taxid.tsv").write_text("\n".join(tax_rows) + "\n")

    with (outdir / "my_traits.tsv").open("w") as fh:
        fh.write(f"species\t{rcfg.traitname}\tfamily\n")
        for sp in sp_sorted:
            fh.write(f"{sp}\t{ph.values[sp]:.6g}\t{fam_of[sp]}\n")
    fg_set = set(ph.foreground_tips)
    with (outdir / "phenotype.tsv").open("w") as fh:
        for sp in sorted(ph.values):
            fh.write(f"{sp}\t{1 if sp in fg_set else 0}\n")

    # ── truth ───────────────────────────────────────────────────────────────
    truth = {
        "seed": seed,
        "archetype": rcfg.archetype,
        "trait_kind": ph.kind,
        "traitname": rcfg.traitname,
        "planted": rcfg.planted,
        "phenotype": {
            "values": {k: round(v, 6) for k, v in ph.values.items()},
            "foreground_tips": ph.foreground_tips,
            "n_origins": ph.n_transitions,
            "top_quantile": ph.top_quantile,
            "bottom_quantile": ph.bottom_quantile,
            "notes": ph.notes,
        },
        "tree": {
            "n_tips": len(tree.tips),
            "total_length": round(tree.total_length(), 6),
        },
        "genes": genes,
        "n_genes": rcfg.n_genes,
        "n_planted_genes": int(n_planted),
        "config": asdict(rcfg),
    }
    (outdir / "truth.json").write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")

    return {
        "dir": outdir.name,
        "seed": seed,
        "archetype": rcfg.archetype,
        "planted": rcfg.planted,
        "n_tips": len(tree.tips),
        "n_genes": rcfg.n_genes,
        "n_planted_genes": int(n_planted),
        "n_origins": ph.n_transitions,
        "trait_kind": ph.kind,
        "top_quantile": ph.top_quantile,
        "bottom_quantile": ph.bottom_quantile,
    }
