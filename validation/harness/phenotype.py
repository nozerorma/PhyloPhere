"""PhenotypeSpec / DatasetSpec — the per-tip phenotype for one validation dataset.

A ``DatasetSpec`` is a tree plus one or more ``PhenotypeSpec`` traits sharing
that tree. ``trait_matrix.emit_trait_matrix`` expands it into the concrete input
files each PhyloPhere trait mode consumes.

Spec files are JSON (``*.spec.json``) so they can live next to the fixtures they
describe and be hashed for provenance. Schema::

    {
      "name": "rh1_spectral",
      "tree": "fixtures/tier1/rh1/tree.nwk",
      "provenance": "Hauser et al. 2017; Yokoyama 2008 — see truthsets/tier1/README.md",
      "traits": [
        {
          "name": "photic_depth",
          "kind": "continuous",
          "values": {"Danio_rerio": 5.0, "Anoplopoma_fimbria": 900.0, ...},
          "higher_is_foreground": true
        },
        {
          "name": "habitat",
          "kind": "categorical",
          "labels": {"Danio_rerio": "fg", "Astyanax_mexicanus": "bg", ...},
          "pairs": [["Danio_rerio", "Astyanax_mexicanus"], ...]
        }
      ]
    }

``kind`` is one of ``categorical`` | ``continuous``. Categorical traits SHOULD
carry ``pairs`` (paired mode is mandatory downstream); if omitted the harness
rank/order-pairs and records that in the manifest.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PhenotypeSpec:
    name: str
    kind: str  # "categorical" | "continuous"
    # continuous: tip -> float ; categorical: unused (see labels)
    values: dict[str, float] = field(default_factory=dict)
    # categorical: tip -> "fg" | "bg"
    labels: dict[str, str] = field(default_factory=dict)
    # categorical: explicit (fg_tip, bg_tip) contrast pairs
    pairs: list[tuple[str, str]] = field(default_factory=list)
    # continuous: which tail is the foreground when percentilized
    higher_is_foreground: bool = True

    def __post_init__(self) -> None:
        if self.kind not in ("categorical", "continuous"):
            raise ValueError(f"{self.name}: kind must be categorical|continuous, got {self.kind!r}")
        if self.kind == "continuous" and not self.values:
            raise ValueError(f"{self.name}: continuous trait needs 'values'")
        if self.kind == "categorical":
            if not self.labels:
                raise ValueError(f"{self.name}: categorical trait needs 'labels'")
            bad = {v for v in self.labels.values() if v not in ("fg", "bg")}
            if bad:
                raise ValueError(f"{self.name}: labels must be fg|bg, got {bad}")

    @property
    def tips(self) -> set[str]:
        return set(self.values) | set(self.labels)

    def foreground(self) -> list[str]:
        if self.kind == "categorical":
            return sorted(t for t, v in self.labels.items() if v == "fg")
        raise TypeError("foreground() is only defined for categorical traits; "
                        "percentilize a continuous trait first")

    def background(self) -> list[str]:
        if self.kind == "categorical":
            return sorted(t for t, v in self.labels.items() if v == "bg")
        raise TypeError("background() is only defined for categorical traits")


@dataclass
class DatasetSpec:
    name: str
    tree: Path
    traits: list[PhenotypeSpec]
    provenance: str = ""
    source_path: Path | None = None

    @property
    def tips(self) -> set[str]:
        out: set[str] = set()
        for t in self.traits:
            out |= t.tips
        return out

    def trait(self, name: str) -> PhenotypeSpec:
        for t in self.traits:
            if t.name == name:
                return t
        raise KeyError(name)


def load_dataset_spec(path: str | Path) -> DatasetSpec:
    path = Path(path)
    raw = json.loads(path.read_text())
    traits = []
    for t in raw["traits"]:
        traits.append(
            PhenotypeSpec(
                name=t["name"],
                kind=t["kind"],
                values={k: float(v) for k, v in t.get("values", {}).items()},
                labels=dict(t.get("labels", {})),
                pairs=[tuple(p) for p in t.get("pairs", [])],
                higher_is_foreground=bool(t.get("higher_is_foreground", True)),
            )
        )
    tree = Path(raw["tree"])
    if not tree.is_absolute():
        tree = (path.parent / tree).resolve()
    return DatasetSpec(
        name=raw["name"],
        tree=tree,
        traits=traits,
        provenance=raw.get("provenance", ""),
        source_path=path,
    )
