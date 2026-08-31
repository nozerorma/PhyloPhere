"""Amino-acid substitution model: PAML matrix parsing + reversible Q + site profiles.

PAML ``.dat`` format (``wag.dat``, ``lg.dat``, ``dayhoff.dat``):
    * 19 lines of lower-triangular exchangeabilities S (line k has k entries,
      k = 1..19, giving S[k][0..k-1] in PAML amino-acid order)
    * one line of 20 equilibrium frequencies pi
    * ``//`` then free-text notes

PAML amino-acid order (columns of every matrix we build):
    A R N D C Q E G H I L K M F P S T W Y V
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

PAML_AA = "ARNDCQEGHILKMFPSTWYV"
AA_INDEX = {aa: i for i, aa in enumerate(PAML_AA)}
N_AA = 20


# --------------------------------------------------------------------------- #
# PAML .dat parsing
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RateMatrix:
    """Exchangeabilities + equilibrium frequencies for a reversible AA model."""

    name: str
    S: np.ndarray   # (20, 20) symmetric, zero diagonal
    pi: np.ndarray  # (20,) sums to 1

    def q(self, pi: np.ndarray | None = None, *, normalise: bool = True) -> np.ndarray:
        """Rate matrix Q for equilibrium ``pi`` (defaults to the model's own).

        Q_ij = S_ij * pi_j  (i != j);  Q_ii = -sum_j Q_ij.
        With ``normalise`` the matrix is scaled so the expected number of
        substitutions per unit branch length is 1 (``-sum_i pi_i Q_ii = 1``),
        which makes a branch length mean the same thing at every site.
        """
        pi = self.pi if pi is None else np.asarray(pi, dtype=float)
        pi = pi / pi.sum()
        Q = self.S * pi[None, :]
        np.fill_diagonal(Q, 0.0)
        np.fill_diagonal(Q, -Q.sum(axis=1))
        if normalise:
            mu = -np.sum(pi * np.diag(Q))
            if mu <= 0:
                raise ValueError(f"{self.name}: non-positive normalisation constant {mu}")
            Q = Q / mu
        return Q


def load_rate_matrix(path: str | Path) -> RateMatrix:
    path = Path(path)
    nums: list[float] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("//", "#", "*")):
            if line.startswith("//"):
                break
            continue
        toks = line.replace(",", " ").split()
        try:
            row = [float(t) for t in toks]
        except ValueError:
            break  # hit the notes
        nums.extend(row)

    tri = N_AA * (N_AA - 1) // 2  # 190
    if len(nums) < tri + N_AA:
        raise ValueError(f"{path}: expected >= {tri + N_AA} numbers, got {len(nums)}")

    S = np.zeros((N_AA, N_AA), dtype=float)
    k = 0
    for i in range(1, N_AA):
        for j in range(i):
            S[i, j] = S[j, i] = nums[k]
            k += 1
    pi = np.array(nums[tri:tri + N_AA], dtype=float)
    if pi.min() < 0 or not np.isclose(pi.sum(), 1.0, atol=1e-3):
        raise ValueError(f"{path}: equilibrium frequencies do not look like a distribution: sum={pi.sum()}")
    return RateMatrix(name=path.stem, S=S, pi=pi / pi.sum())


# --------------------------------------------------------------------------- #
# site heterogeneity
# --------------------------------------------------------------------------- #
def site_profiles(base_pi: np.ndarray, n_sites: int, concentration: float,
                  rng: np.random.Generator) -> np.ndarray:
    """Draw ``n_sites`` equilibrium profiles ~ Dirichlet(concentration * base_pi).

    ``concentration`` large  -> profiles ~= base_pi (site-homogeneous, WAG-like)
    ``concentration`` small  -> peaked idiosyncratic profiles (conserved sites
    with site-specific preferences); this is what makes non-adaptive convergence
    occur and is the property the pipeline's single-profile ASR does not model.
    """
    if concentration <= 0:
        raise ValueError("concentration must be > 0")
    alpha = concentration * (base_pi / base_pi.sum())
    alpha = np.clip(alpha, 1e-6, None)
    prof = rng.dirichlet(alpha, size=n_sites)
    return np.clip(prof, 1e-9, None) / np.clip(prof, 1e-9, None).sum(axis=1, keepdims=True)


def gamma_rates(n_sites: int, alpha: float, rng: np.random.Generator,
                *, categories: int | None = None) -> np.ndarray:
    """Per-site relative rates with mean 1.

    Continuous Gamma(alpha, 1/alpha) by default; set ``categories`` for the
    discrete-gamma (mean of each equal-probability bin) approximation PhyloPhere's
    inference uses, so the rate model is shared but the matrix family is not.
    """
    if alpha <= 0:
        raise ValueError("gamma alpha must be > 0")
    if categories is None:
        r = rng.gamma(shape=alpha, scale=1.0 / alpha, size=n_sites)
        return r / r.mean()
    # discrete gamma: category means, then assign sites round-robin (shuffled)
    from math import gamma as _g  # noqa: F401

    try:
        from scipy.stats import gamma as _sg  # type: ignore

        edges = _sg.ppf(np.linspace(0, 1, categories + 1), a=alpha, scale=1.0 / alpha)
        means = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            xs = np.linspace(max(lo, 1e-9), hi if np.isfinite(hi) else lo + 20, 512)
            w = _sg.pdf(xs, a=alpha, scale=1.0 / alpha)
            means.append(np.trapz(xs * w, xs) / np.trapz(w, xs))
        means = np.array(means)
    except Exception:
        means = np.sort(rng.gamma(alpha, 1.0 / alpha, size=categories))
    means = means / means.mean()
    idx = np.resize(np.arange(categories), n_sites)
    rng.shuffle(idx)
    return means[idx]
