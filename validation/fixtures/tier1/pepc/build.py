"""Sanitise the raw PCOC/ConDor PEPC C3/C4 sedge dataset into a Tier 1 fixture.

Raw inputs (committed, fetched from github.com/evolbioinfo/condor/test_data,
originally Besnard et al. 2009, MBE 26:1909, courtesy of the authors via PCOC):

    pepc.aa.coor_mays.fa            79 Cyperaceae PEPC (ppc) amino-acid sequences,
                                    970 columns == maize PEPC1 (UniProt P04711,
                                    exactly 970 aa) coordinates, 1:1.
    pepc.phyml_tree.txt             PhyML tree, 79 tips, aLRT support on nodes.
    besnard2009_convergent_species  23 tips with the "genotypic" C4 annotation
                                    (presence of the A780S determinant).
    outgroup.txt                    Chrysithrix (root; excluded from fg/bg).

Output (pipeline-shaped, tips = sanitised species names, "." -> "_"):

    align/PEPC.fasta
    tree.nwk
    my_traits.tsv        species <tab> c4 (1|0) <tab> family
    phenotype.tsv        species <tab> 1|0     (C4 foreground, outgroup dropped)
    ali_sp_names.txt / gene_ensembl.tsv / taxid.tsv    synthetic support files

Numbering note: the truth set (validation/truthsets/tier1/pepc_c4.sites.tsv) is
in maize P04711 coordinates, which equal alignment columns 1:1 here.
"""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent
RAW_ALN = HERE / "pepc.aa.coor_mays.fa"
RAW_TREE = HERE / "pepc.phyml_tree.txt"
C4_FILE = HERE / "besnard2009_convergent_species.txt"
OG_FILE = HERE / "outgroup.txt"

# Pinned raw inputs (gitignored — fetched on demand). github.com/evolbioinfo/
# condor test_data, commit-agnostic master; originally Besnard et al. 2009.
_BASE = "https://raw.githubusercontent.com/evolbioinfo/condor/master/test_data/"
_RAW = {
    RAW_ALN: "cyp_coding.aa.coor_mays.fa",
    RAW_TREE: "cyp_coding.phy_phyml_tree.txt",
    C4_FILE: "besnard2009_convergent_species.txt",
    OG_FILE: "outgroup.txt",
}


def _fetch() -> None:
    for dst, name in _RAW.items():
        if dst.exists():
            continue
        print(f"fetch {name}")
        with urllib.request.urlopen(_BASE + name, timeout=60) as r:
            dst.write_bytes(r.read())


def _san(name: str) -> str:
    return name.replace(".", "_")


def _read_fasta(path: Path) -> dict[str, str]:
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


def main() -> None:
    _fetch()
    aln = {_san(k): v for k, v in _read_fasta(RAW_ALN).items()}
    tree = RAW_TREE.read_text().strip()
    # sanitise tip names in the newick (they appear as "(<name>:" or ",<name>:")
    tree = re.sub(r"([(,])([A-Za-z0-9._-]+?):",
                  lambda m: f"{m.group(1)}{_san(m.group(2))}:", tree)

    c4 = {_san(x) for x in C4_FILE.read_text().split()}
    outgroup = {_san(x) for x in OG_FILE.read_text().split()}
    assert c4 <= set(aln), c4 - set(aln)
    assert outgroup <= set(aln), outgroup - set(aln)

    L = len(next(iter(aln.values())))
    assert L == 970, f"expected 970 maize-coordinate columns, got {L}"
    assert all(len(s) == L for s in aln.values())

    (HERE / "align").mkdir(exist_ok=True)
    with (HERE / "align" / "PEPC.fasta").open("w") as fh:
        for name, seq in sorted(aln.items()):
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    (HERE / "tree.nwk").write_text(tree + "\n")

    species = sorted(aln)
    fam_of = {sp: "Cyperaceae" for sp in species}
    with (HERE / "my_traits.tsv").open("w") as fh:
        fh.write("species\tc4\tfamily\n")
        for sp in species:
            if sp in outgroup:
                continue
            fh.write(f"{sp}\t{1 if sp in c4 else 0}\t{fam_of[sp]}\n")
    with (HERE / "phenotype.tsv").open("w") as fh:
        for sp in species:
            if sp in outgroup:
                continue
            fh.write(f"{sp}\t{1 if sp in c4 else 0}\n")

    (HERE / "ali_sp_names.txt").write_text("\n".join(species) + "\n")
    (HERE / "gene_ensembl.tsv").write_text(
        "gene\tchr\tstart\tend\tstrand\tlength\thuman_protein_id\n"
        f"PEPC\tchr1\t100000\t{100000 + L * 3}\t+\t{L}\tP04711\n"
    )
    with (HERE / "taxid.tsv").open("w") as fh:
        fh.write("tax_id\tspecies\tfamily\trank\tname_class\n")
        for i, sp in enumerate(species):
            fh.write(f"{9100001 + i}\t{sp}\tCyperaceae\tspecies\tscientific name\n")

    n_c4 = sum(1 for sp in species if sp in c4)
    print(f"PEPC fixture: {len(species)} tips  "
          f"({n_c4} C4, {len(species) - n_c4 - len(outgroup)} C3, "
          f"{len(outgroup)} outgroup)  x {L} maize-coord columns")


if __name__ == "__main__":
    main()
