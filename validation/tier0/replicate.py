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
    archetype: str = "binary"            # "binary" (0/1 code) | "rate" (c/n + counts)
    lam: float = 0.5                     # Pagel's lambda for the latent trait (0 | 0.5 | 1)
    n_pairs: int = 4                     # FIXED contrast count (pipeline contrast_max_iter=3 -> 4)
    n_genes: int = 150
    frac_planted_genes: float = 0.2
    traitname: str = "sim_trait"
    # per-gene sequence model (SimConfig knobs)
    n_sites: int = 400
    concentration: float = 2.0           # lower = more profile heterogeneity (realistic)
    gamma_alpha: float = 0.7
    matrix: str = "wag"
    # planted-gene signal — two mechanisms, scored separately (D-T0-02):
    #   identical_aa  : same residue across fg lineages  -> US + all GS
    #   grouped_caap  : same GS class, residue free      -> GS only, US should not
    n_planted_identical_aa: int = 12
    n_planted_grouped_caap: int = 12
    grouped_caap_scheme: str = "GS1"
    identical_aa_target_weight: float = 0.95
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

    # ── phenotype: BM-under-lambda latent -> independent tail pairs ────────
    if rcfg.archetype not in ("binary", "rate"):
        raise ValueError(f"unknown archetype {rcfg.archetype!r}")
    ph = _pheno.make_lambda_foreground(
        tree, rcfg.n_pairs, rng,
        kind=("binary" if rcfg.archetype == "binary" else "rate"),
        planted=rcfg.planted, lam=rcfg.lam,
    )

    # ── genes ───────────────────────────────────────────────────────────────
    n_planted = 0 if not rcfg.planted else max(1, round(rcfg.n_genes * rcfg.frac_planted_genes))
    planted_flags = np.zeros(rcfg.n_genes, dtype=bool)
    planted_flags[:n_planted] = True
    rng.shuffle(planted_flags)
    # per planted gene: alternate top / bottom planting direction
    directions = ["top" if i % 2 == 0 else "bottom" for i in range(n_planted)]
    rng.shuffle(directions)
    fg_top, fg_bot = ph.foreground("top"), ph.foreground("bottom")

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
        direction = directions[sum(planted_flags[:gi])] if is_planted else None
        scfg = SimConfig(
            n_sites=rcfg.n_sites, concentration=rcfg.concentration,
            gamma_alpha=rcfg.gamma_alpha, matrix=rcfg.matrix,
            n_planted_profile_shift=0,
            n_planted_identical_aa=rcfg.n_planted_identical_aa if is_planted else 0,
            n_planted_grouped_caap=rcfg.n_planted_grouped_caap if is_planted else 0,
            grouped_caap_scheme=rcfg.grouped_caap_scheme,
            identical_aa_target_weight=rcfg.identical_aa_target_weight,
        )
        gene_fg = (fg_top if direction == "top" else fg_bot) if is_planted else None
        res = simulate(tree, gene_fg, scfg, seed=seed * 100_000 + gi)
        _write_fasta(outdir / "align" / f"{gid}.fasta", res.alignment)
        genes[gid] = {
            "planted": is_planted,
            "direction": direction,
            "planted_sites": {str(k): v for k, v in res.truth.planted_sites.items()},
            "identical_aa_targets": res.truth.identical_aa_targets,
            "grouped_caap_targets": res.truth.grouped_caap_targets,
            "n_planted": len(res.truth.planted_sites),
        }

    # ── phenotype / trait files ─────────────────────────────────────────────
    _write_newick(tree, outdir / "tree.nwk")
    # species list — NAME_CURATION scans every alignment to derive this if absent
    (outdir / "ali_sp_names.txt").write_text("\n".join(sorted(tree.tips)) + "\n")
    (outdir / "gene_ensembl.tsv").write_text("\n".join(ensembl_rows) + "\n")

    all_sp = sorted(tree.tips)
    # `family` — the reporting / CI Rmds need a `taxon_of_interest` column; the
    # taxid map needs a `family` column too. Synthetic: 5 fake families.
    fam_of = {sp: f"simfam{i % 5 + 1}" for i, sp in enumerate(all_sp)}

    # taxid mapping — disambiguation renames tips to numeric ids before writing
    # the PAML alignment (reconstruct.py::_write_phylip uses a single-space name
    # delimiter that breaks on species names containing letters like 'O';
    # production always supplies --tax_id so this path is what real runs use).
    tax_rows = ["tax_id\tspecies\tfamily\trank\tname_class"]
    for i, sp in enumerate(all_sp):
        tax_rows.append(f"{9_000_001 + i}\t{sp}\t{fam_of[sp]}\tspecies\tscientific name")
    (outdir / "taxid.tsv").write_text("\n".join(tax_rows) + "\n")

    # my_traits.tsv — binary: every tip a 0/1 code. rate: only species with data,
    # plus n_pop / n_cases count columns (CLASS 1 Jeffreys-CI path).
    data_sp = sorted(ph.values)
    with (outdir / "my_traits.tsv").open("w") as fh:
        if rcfg.archetype == "rate":
            fh.write(f"species\t{rcfg.traitname}\tn_pop\tn_cases\tfamily\n")
            for sp in data_sp:
                fh.write(f"{sp}\t{ph.values[sp]:.6g}\t{ph.n_pop[sp]}\t{ph.n_cases[sp]}\t{fam_of[sp]}\n")
        else:
            fh.write(f"species\t{rcfg.traitname}\tfamily\n")
            for sp in all_sp:
                fh.write(f"{sp}\t{int(ph.values[sp])}\t{fam_of[sp]}\n")
    fg_set = set(ph.foreground_tips)
    with (outdir / "phenotype.tsv").open("w") as fh:
        for sp in data_sp:
            fh.write(f"{sp}\t{1 if sp in fg_set else 0}\n")

    # ── truth ───────────────────────────────────────────────────────────────
    truth = {
        "seed": seed,
        "archetype": rcfg.archetype,
        "trait_kind": ph.kind,
        "traitname": rcfg.traitname,
        "planted": rcfg.planted,
        "phenotype": {
            "lambda": ph.lam,
            "values": {k: round(v, 6) for k, v in ph.values.items()},
            "foreground_tips": ph.foreground_tips,          # true top anchors
            "partner_tips": sorted(b for _, b in ph.pairs),  # true bottom partners
            "pairs": [list(pr) for pr in ph.pairs],
            "n_origins": ph.n_transitions,
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
        "lambda": rcfg.lam,
        "planted": rcfg.planted,
        "n_tips": len(tree.tips),
        "n_genes": rcfg.n_genes,
        "n_planted_genes": int(n_planted),
        "n_origins": ph.n_transitions,
        "trait_kind": ph.kind,
    }
