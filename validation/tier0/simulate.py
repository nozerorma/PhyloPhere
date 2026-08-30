"""Evolve an amino-acid alignment down a tree and inject planted convergence.

Generative model (see ``tier0/__init__.py`` for why it differs from inference):
    * WAG exchangeabilities  (``data/wag.dat``)
    * per-site equilibrium profile ~ Dirichlet(concentration * pi_WAG)
    * continuous gamma rates across sites
    * simulation by exact Gillespie along each branch (no matrix exponential;
      trivially handles branch-specific rate matrices)

Planted positives, recorded per site with their mechanism so recall can be
reported per mechanism:
    * ``identical_aa``   — on foreground edges the site's equilibrium is a near
      point mass on one target residue (the least-favoured residue under the
      site's background profile), so every foreground lineage converges to the
      *same* amino acid: the strict-CAAS (US) positive, and trivially also a
      grouped-CAAP (GS*) hit since the shared residue is in one group.
      The target is hard-set on each origin edge, then held by Q_fg.
    * ``grouped_caap``   — on foreground edges the equilibrium is ~uniform over
      one physicochemical *class* of a GS scheme (default GS1), the class chosen
      as the one most disfavoured in the site background. A fresh residue from
      the class is seeded per origin edge, so origins share the CLASS but
      usually differ in RESIDUE: a GS1-GS4 positive that US should NOT call.
    * ``profile_shift``  — DORMANT (``n_planted_profile_shift`` defaults to 0 and
      the replicate driver never sets it). On foreground edges the site would use
      a *different* Dirichlet profile — a preference shift with no shared target
      residue *or* group, so it is neither a US positive nor a clean CAAP
      positive. Kept as a hook for a future grouped-CAAP planted mechanism
      (target a shared ``amino_encoded`` physicochemical group, residue free
      within it) — see DECISIONS.md D-T0-N.

Output per replicate (``<outdir>/``):
    aln.fasta        amino-acid alignment, tips only, PAML/uppercase
    phenotype.tsv    ``species<tab>state``   state in {1,0}   (foreground = 1)
    truth.json       planted sites + mechanism, foreground tips/edges/transitions,
                     every parameter, the RNG seed
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .groups import GS_SCHEMES, SCHEMES, groups_of
from .model import AA_INDEX, N_AA, PAML_AA, RateMatrix, gamma_rates, load_rate_matrix, site_profiles
from .trees import Foreground, PhyloTree

_DATA = Path(__file__).parent / "data"


# --------------------------------------------------------------------------- #
# config
# --------------------------------------------------------------------------- #
@dataclass
class SimConfig:
    n_sites: int = 400
    concentration: float = 4.0          # Dirichlet concentration for site profiles
    gamma_alpha: float = 0.7            # across-site rate heterogeneity
    matrix: str = "wag"                 # data/<matrix>.dat  (NOT "lg" — see module doc)
    n_planted_profile_shift: int = 0    # DORMANT — see module doc
    n_planted_identical_aa: int = 0     # sites converging to the same residue (US positive)
    n_planted_grouped_caap: int = 0     # sites converging to a shared GS class (GS-only positive)
    identical_aa_target_weight: float = 0.95  # equilibrium mass on the target residue / group
    grouped_caap_scheme: str = "GS1"    # scheme whose class the grouped sites converge to
    shift_concentration: float = 2.0    # concentration of the alternate (shifted) profile

    def matrix_path(self) -> Path:
        return _DATA / f"{self.matrix}.dat"


@dataclass
class Truth:
    planted_sites: dict[int, str] = field(default_factory=dict)   # 0-based col -> mechanism
    identical_aa_targets: dict[int, str] = field(default_factory=dict)  # col -> residue
    grouped_caap_targets: dict[int, dict] = field(default_factory=dict)  # col -> {scheme, group, residues}
    foreground_tips: list[str] = field(default_factory=list)
    n_transitions: int = 0
    n_foreground_edges: int = 0
    seed: int = 0
    config: dict = field(default_factory=dict)
    tree_total_length: float = 0.0
    n_tips: int = 0

    def to_json(self) -> str:
        d = asdict(self)
        d["planted_sites"] = {str(k): v for k, v in self.planted_sites.items()}
        d["identical_aa_targets"] = {str(k): v for k, v in self.identical_aa_targets.items()}
        d["grouped_caap_targets"] = {str(k): v for k, v in self.grouped_caap_targets.items()}
        return json.dumps(d, indent=2, sort_keys=True) + "\n"


@dataclass
class SimResult:
    tree: PhyloTree
    alignment: dict[str, str]   # tip label -> aa string
    phenotype: dict[str, int]   # tip label -> 1 (fg) | 0 (bg)
    truth: Truth

    def write(self, outdir: str | Path) -> Path:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        with (outdir / "aln.fasta").open("w") as fh:
            for name, seq in self.alignment.items():
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
        with (outdir / "phenotype.tsv").open("w") as fh:
            for name, st in sorted(self.phenotype.items()):
                fh.write(f"{name}\t{st}\n")
        (outdir / "truth.json").write_text(self.truth.to_json())
        return outdir


# --------------------------------------------------------------------------- #
# per-site rate matrices
# --------------------------------------------------------------------------- #
def _point_mass_profile(target: int, weight: float) -> np.ndarray:
    pi = np.full(N_AA, (1.0 - weight) / (N_AA - 1))
    pi[target] = weight
    return pi


def _group_profile(members: list[int], weight: float) -> np.ndarray:
    """Equilibrium ~uniform over ``members`` (total mass ``weight``), tiny elsewhere."""
    pi = np.full(N_AA, (1.0 - weight) / (N_AA - len(members)))
    for m in members:
        pi[m] = weight / len(members)
    return pi


def _build_site_matrices(cfg: SimConfig, rm: RateMatrix, rng: np.random.Generator):
    """Return (Q_bg, Q_fg, targets) where Q_bg/Q_fg are (n_sites, 20, 20) stacks
    of normalised rate matrices for background / foreground edges, and ``targets``
    maps planted identical-aa site -> residue index."""
    n = cfg.n_sites
    base_prof = site_profiles(rm.pi, n, cfg.concentration, rng)

    Q_bg = np.stack([rm.q(base_prof[s]) for s in range(n)])
    Q_fg = Q_bg.copy()

    planted: dict[int, str] = {}
    targets: dict[int, int] = {}
    group_targets: dict[int, dict] = {}

    scheme = cfg.grouped_caap_scheme
    scheme_groups = groups_of(scheme, min_size=3)  # label -> [residues]

    all_sites = rng.permutation(n)
    cur = 0
    for _ in range(cfg.n_planted_profile_shift):
        s = int(all_sites[cur]); cur += 1
        alt = site_profiles(rm.pi, 1, cfg.shift_concentration, rng)[0]
        Q_fg[s] = rm.q(alt)
        planted[s] = "profile_shift"
    for _ in range(cfg.n_planted_identical_aa):
        s = int(all_sites[cur]); cur += 1
        # target: a residue currently disfavoured at this site, to make the
        # convergent signal unambiguous
        tgt = int(np.argmin(base_prof[s]))
        Q_fg[s] = rm.q(_point_mass_profile(tgt, cfg.identical_aa_target_weight))
        planted[s] = "identical_aa"
        targets[s] = tgt
    for _ in range(cfg.n_planted_grouped_caap):
        s = int(all_sites[cur]); cur += 1
        # pick the scheme class whose summed background probability at this site
        # is lowest, so foreground convergence onto it is unambiguous; the
        # residue stays free inside the class (US sees divergence, GS* sees one
        # shared class).
        label = min(
            scheme_groups,
            key=lambda g: base_prof[s][[AA_INDEX[a] for a in scheme_groups[g]]].sum(),
        )
        members = [AA_INDEX[a] for a in scheme_groups[label]]
        Q_fg[s] = rm.q(_group_profile(members, cfg.identical_aa_target_weight))
        planted[s] = "grouped_caap"
        group_targets[s] = {
            "scheme": scheme, "group": label,
            "residues": sorted(scheme_groups[label]),
        }

    return Q_bg, Q_fg, base_prof, planted, targets, group_targets


# --------------------------------------------------------------------------- #
# Gillespie
# --------------------------------------------------------------------------- #
def _evolve_branch(states: np.ndarray, Q: np.ndarray, T: np.ndarray,
                   rng: np.random.Generator) -> np.ndarray:
    """Advance every site independently along one branch.

    states : (n_sites,) int         current residue per site
    Q      : (n_sites, 20, 20)       rate matrix per site for THIS branch
    T      : (n_sites,)              effective time = branch_length * site_rate
    """
    states = states.copy()
    remaining = T.astype(float).copy()
    active = np.where(remaining > 0)[0]
    # exit rate for the current state of each site
    guard = 0
    while active.size:
        guard += 1
        if guard > 10_000:
            raise RuntimeError("Gillespie did not terminate — check Q normalisation")
        rows = Q[active, states[active]]                      # (m, 20)
        exit_rate = -rows[np.arange(active.size), states[active]]
        exit_rate = np.maximum(exit_rate, 1e-12)
        dt = rng.exponential(1.0 / exit_rate)
        jumped = dt < remaining[active]
        # sites that jump: choose destination
        j_idx = active[jumped]
        if j_idx.size:
            jrows = Q[j_idx, states[j_idx]].copy()            # (mj, 20)
            jrows[np.arange(j_idx.size), states[j_idx]] = 0.0
            jrows = np.maximum(jrows, 0.0)
            jrows /= jrows.sum(axis=1, keepdims=True)
            cdf = np.cumsum(jrows, axis=1)
            u = rng.random(j_idx.size)
            states[j_idx] = (u[:, None] < cdf).argmax(axis=1)
        remaining[active] = np.where(jumped, remaining[active] - dt, 0.0)
        active = np.where(remaining > 0)[0]
    return states


# --------------------------------------------------------------------------- #
# top level
# --------------------------------------------------------------------------- #
def simulate(tree: PhyloTree, fg: Foreground | None, cfg: SimConfig,
             seed: int) -> SimResult:
    rng = np.random.default_rng(seed)
    rm = load_rate_matrix(cfg.matrix_path())

    Q_bg, Q_fg, base_prof, planted, targets, group_targets = _build_site_matrices(cfg, rm, rng)
    rates = gamma_rates(cfg.n_sites, cfg.gamma_alpha, rng)

    fg_edge_set = fg.fg_edges if fg else set()

    # root sequence: draw each site from its background equilibrium profile
    node_states: dict[int, np.ndarray] = {}
    root_seq = np.array([
        rng.choice(N_AA, p=base_prof[s]) for s in range(cfg.n_sites)
    ], dtype=np.int64)
    node_states[tree.root] = root_seq

    # identical-aa planting: at the first foreground edge of each origin, seed the
    # target residue so every independent lineage starts converged; drift within
    # the clade is then governed by Q_fg (near point mass -> it stays).
    origin_edges: set[int] = set()
    if fg:
        parent_of = {c: p for p, c, _ in tree.edges}
        edge_index = {(p, c): i for i, (p, c, bl) in enumerate(tree.edges)}
        for node in fg.origin_nodes:
            if node in parent_of:
                origin_edges.add(edge_index[(parent_of[node], node)])

    for ei, (p, c, bl) in enumerate(tree.edges):
        parent = node_states[p]
        is_fg = ei in fg_edge_set
        Q = Q_fg if is_fg else Q_bg
        T = bl * rates
        child = parent.copy()
        if ei in origin_edges:
            for s, tgt in targets.items():
                child[s] = tgt
            # grouped-CAAP: seed a residue drawn from the target class — a fresh
            # draw per origin edge, so origins land on the SAME class but usually
            # DIFFERENT residues (US divergent, GS* convergent).
            for s, gt in group_targets.items():
                child[s] = int(rng.choice([AA_INDEX[a] for a in gt["residues"]]))
        child = _evolve_branch(child, Q, T, rng)
        node_states[c] = child

    leaves = tree._leaves()
    alignment = {
        tree.labels[i]: "".join(PAML_AA[a] for a in node_states[i]) for i in leaves
    }
    fg_tips = set(fg.tips) if fg else set()
    phenotype = {tree.labels[i]: (1 if tree.labels[i] in fg_tips else 0) for i in leaves}

    truth = Truth(
        planted_sites=planted,
        identical_aa_targets={s: PAML_AA[t] for s, t in targets.items()},
        grouped_caap_targets=group_targets,
        foreground_tips=sorted(fg_tips),
        n_transitions=fg.n_transitions if fg else 0,
        n_foreground_edges=len(fg_edge_set),
        seed=seed,
        config=asdict(cfg),
        tree_total_length=tree.total_length(),
        n_tips=len(leaves),
    )
    return SimResult(tree=tree, alignment=alignment, phenotype=phenotype, truth=truth)
