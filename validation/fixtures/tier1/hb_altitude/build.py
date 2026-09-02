"""Build the Tier 1 high-altitude haemoglobin fixture from Zhu et al. 2018 PNAS.

Zhu X, Guan Y, Signore AV, Natarajan C, ... Storz JF (2018). "Divergent and
parallel routes of biochemical adaptation in high-altitude passerine birds from
the Qinghai-Tibet Plateau." PNAS 115(8):1865-1870. doi:10.1073/pnas.1720487115

Phenotype: altitude (Qinghai-Tibet Plateau tits + long-tailed tits). The paper
analyses it categorically (high / low); the fixture ALSO carries a continuous
elevation trait (per-species range midpoint) - that continuous encoding is ours,
see README.md.

Truth (validation/truthsets/tier1/hb_altitude.sites.tsv):
  - alphaA A34T : parallel, mutagenesis-confirmed, in Pseudopodoces humilis
                  (= Parus humilis) and Lophophanes dichrous. The alphaA GENE
                  tree (Fig S2A) shows the two Thr34 states arose independently.
  - alphaA P119A: mutagenesis-confirmed affinity increase in Aegithalos bonvaloti
                  (parallel with the bar-headed goose, which is outside this set).
  - divergent tier: the AncParidae -> P. humilis / L. dichrous substitution lists
                    (Fig S3); not individually functionally validated.

Raw input (gitignored, fetched on demand):
  hb_genbank.gb   GenBank records MG772099-MG772439 (341 seqs: alphaA/alphaD/betaA
                  globin, ~3-16 population isolates per species per gene).

Output (pipeline-shaped, tips = binomials, "." / " " -> "_"):
  align/HBA.fasta  align/HBD.fasta  align/HBB.fasta   per-gene AA alignments
                   (initiator Met stripped; for HBA/HBB, column N == mature-chain
                   residue N, 1:1, so the truth set needs no ref row - verified
                   in _assert_landmarks below)
  tree.nwk         species tree: Fig S1 topology (Johansson 2013 + Li 2016),
                   ML branch lengths on the concatenation, then rooted on the
                   Aegithalidae outgroup and time-scaled to an ultrametric
                   chronogram (ape::chronos, penalised likelihood) — PhyloPhere
                   contrast selection assumes a time tree (see _datetree below)
  tree_substitution.nwk   the raw ML phylogram, kept for provenance
  gene_trees/HBA.nwk HBD.nwk HBB.nwk   unconstrained per-gene ML trees
                   (reproduce the Fig S2 genealogical discordance)
  my_traits.tsv    species <tab> elev_mid <tab> altitude(H|L) <tab> family
  phenotype.tsv    species <tab> 1|0   (1 = high-altitude)
  ali_sp_names.txt / gene_ensembl.tsv / taxid.tsv     synthetic support files

Needs the fixture-build env (mafft, iqtree):
  micromamba env create -f validation/fixtures/fixture-tools.yml
  micromamba run -n phylophere-fixtures python3 validation/fixtures/tier1/hb_altitude/build.py
"""

from __future__ import annotations

import subprocess
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).parent
RAW_GB = HERE / "hb_genbank.gb"

# MG772099-MG772439 : Zhu et al. 2018 globin sequences (NCBI Nucleotide).
_ACC = [f"MG772{n:03d}" for n in range(99, 440)]
_EFETCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&rettype=gb&retmode=text&id="
)

# ---------------------------------------------------------------------------- #
# focal taxa  ==  Fig S1 (16 species).  Elevational range (m): the 7 values in
# bold are quoted verbatim in the main text; the rest are read off the Fig S1
# bars (+/- ~150 m, checked by MR 2026-08-31).  trait = midpoint.
#   name in fixture              low   high   text-exact?
FOCAL: dict[str, tuple[int, int, bool, str]] = {
    "Pseudopodoces_humilis":     (3100, 5500, True,  "Paridae"),   # = Parus humilis
    "Lophophanes_dichrous":      (2300, 4600, True,  "Paridae"),
    "Periparus_rubidiventris":   (2400, 4300, True,  "Paridae"),
    "Periparus_ater_aemodius":   (2100, 4600, True,  "Paridae"),
    "Aegithalos_bonvaloti":      (1500, 4600, False, "Aegithalidae"),
    "Poecile_davidi":            (2300, 3500, False, "Paridae"),
    "Poecile_montanus":          (1500, 4000, False, "Paridae"),
    "Parus_monticolus":          (1500, 3300, False, "Paridae"),
    "Sylviparus_modestus":       ( 700, 3300, False, "Paridae"),
    "Pardaliparus_venustulus":   ( 400, 3400, False, "Paridae"),
    "Machlolophus_spilonotus":   ( 900, 2800, False, "Paridae"),   # = Parus spilonotus
    "Aegithalos_fuliginosus":    ( 900, 2500, False, "Aegithalidae"),
    "Cyanistes_cyanus":          ( 600, 2400, False, "Paridae"),
    "Parus_minor":               (   0, 2000, True,  "Paridae"),
    "Poecile_palustris":         (   0, 2100, True,  "Paridae"),
    "Periparus_ater_pekinensis": (   0, 1800, True,  "Paridae"),
}
# GenBank organism string -> fixture name
_ORG2NAME = {n.replace("_", " "): n for n in FOCAL}

# high / low split as used by Zhu et al. (Table S2 (H)/(L); here by midpoint,
# which reproduces their assignment: cut at 2500 m).
_HIGH_CUT = 2500

# Fig S1 topology (Johansson et al. 2013 MPE 69:852; Li et al. 2016 MPE 104:14).
# Sylviparus basal in Paridae; Pseudopodoces nested in the Parus group; the two
# Aegithalidae form the outgroup.  Branch lengths are added by IQ-TREE below.
_TOPOLOGY = (
    "((Aegithalos_bonvaloti,Aegithalos_fuliginosus),"
    "(Sylviparus_modestus,"
    "((Pardaliparus_venustulus,(Periparus_rubidiventris,"
    "(Periparus_ater_aemodius,Periparus_ater_pekinensis))),"
    "((Lophophanes_dichrous,(Poecile_montanus,(Poecile_palustris,Poecile_davidi))),"
    "(Cyanistes_cyanus,(Pseudopodoces_humilis,"
    "(Machlolophus_spilonotus,(Parus_minor,Parus_monticolus))))))));"
)

_GENES = ("HBA", "HBD", "HBB")
_MATURE_LEN = {"HBA": 141, "HBD": 141, "HBB": 146}


# ---------------------------------------------------------------------------- #
def _fetch() -> None:
    if RAW_GB.exists():
        return
    print(f"fetch {len(_ACC)} GenBank records -> {RAW_GB.name}")
    buf = []
    for i in range(0, len(_ACC), 100):
        chunk = ",".join(_ACC[i : i + 100])
        with urllib.request.urlopen(_EFETCH + chunk, timeout=180) as r:
            buf.append(r.read().decode())
    RAW_GB.write_text("".join(buf))


def _run(env_cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(env_cmd, check=True, capture_output=True, text=True, **kw)


def _rmtree(p: Path) -> None:
    import shutil

    if p.exists():
        shutil.rmtree(p)


def _gene_of(desc: str) -> str:
    for g in _GENES:
        if f"({g})" in desc:
            return g
    return "?"


def _consensus(seqs: list[str]) -> str:
    """Per-column majority over a species' allele set (ties -> first seen)."""
    L = max(len(s) for s in seqs)
    out = []
    for i in range(L):
        col = [s[i] for s in seqs if i < len(s)]
        out.append(Counter(col).most_common(1)[0][0])
    return "".join(out)


def _load_proteins() -> dict[str, dict[str, str]]:
    from Bio import SeqIO

    raw: dict[str, dict[str, list[str]]] = {g: defaultdict(list) for g in _GENES}
    for rec in SeqIO.parse(RAW_GB, "genbank"):
        org = rec.annotations.get("organism", "")
        name = _ORG2NAME.get(org)
        if name is None:
            continue  # non-focal taxon (extra Aegithalos, pooled "Periparus ater")
        g = _gene_of(rec.description)
        if g == "?":
            continue
        for f in rec.features:
            if f.type == "CDS" and "translation" in f.qualifiers:
                raw[g][name].append(f.qualifiers["translation"][0])
                break

    out: dict[str, dict[str, str]] = {}
    for g in _GENES:
        out[g] = {}
        for name in FOCAL:
            seqs = raw[g].get(name, [])
            if not seqs:
                raise SystemExit(f"{g}: no sequence for {name}")
            cons = _consensus(seqs)
            # alphaA and betaA: initiator Met is cleaved from the mature chain,
            # so strip it -> column N == mature residue N (the truth-set frame).
            # alphaD: the N-terminal Met is retained in the mature chain -> keep.
            if g in ("HBA", "HBB") and cons[0] == "M":
                cons = cons[1:]
            out[g][name] = cons
        print(f"  {g}: {len(out[g])} species  "
              f"(len {sorted({len(s) for s in out[g].values()})})")
    return out


def _align(seqs: dict[str, str], gene: str, outdir: Path) -> dict[str, str]:
    unaln = outdir / f"{gene}.unaln.fasta"
    with unaln.open("w") as fh:
        for name, s in sorted(seqs.items()):
            fh.write(f">{name}\n{s}\n")
    res = _run(["mafft", "--localpair", "--maxiterate", "1000",
                "--quiet", str(unaln)])
    aln: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for ln in res.stdout.splitlines():
        if ln.startswith(">"):
            if name:
                aln[name] = "".join(buf)
            name, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if name:
        aln[name] = "".join(buf)
    unaln.unlink()
    return aln


def _assert_landmarks(aln: dict[str, dict[str, str]]) -> None:
    """The two mutagenesis-confirmed truth sites must land where the truth set
    says, in raw alignment-column coordinates (== mature residue, 1-based)."""
    hba = aln["HBA"]
    # no gaps expected in the tit HBA block -> column == residue
    if any("-" in s for s in hba.values()):
        raise SystemExit("HBA alignment has gaps - truth-set column mapping is "
                         "no longer identity; add a ref row / mapping step")
    col34 = {n: s[33] for n, s in hba.items()}       # 1-based residue 34
    col119 = {n: s[118] for n, s in hba.items()}
    thr34 = {n for n, a in col34.items() if a == "T"}
    ala119 = {n for n, a in col119.items() if a == "A"}
    # humilis + dichrous carry the mutagenesis-confirmed derived Thr34.
    # Sylviparus modestus (low-altitude, basal Paridae) ALSO shows Thr34 - but
    # that is retention of the ancestral avian state (chicken alphaA is Thr34
    # too; AncParidae is the derived Ala34).  It is a genuine background species
    # carrying the focal residue; the fixture keeps it (README "Known
    # complications") because resolving retention-vs-reversal is exactly the
    # CT_DISAMBIGUATION / ASR job.
    truth34 = {"Pseudopodoces_humilis", "Lophophanes_dichrous"}
    assert truth34 <= thr34, f"missing Thr34: {truth34 - thr34}"
    assert thr34 <= truth34 | {"Sylviparus_modestus"}, f"unexpected Thr34: {thr34}"
    assert ala119 == {"Aegithalos_bonvaloti"}, ala119
    print(f"  landmark OK: alphaA34 Thr in {sorted(thr34)} "
          f"(truth: {sorted(truth34)}); alphaA119 Ala in {sorted(ala119)}")


_OUTGROUP = ("Aegithalos_bonvaloti", "Aegithalos_fuliginosus")  # Aegithalidae


def _iqtree(concat: Path, per_gene: dict[str, Path], outdir: Path) -> None:
    work = outdir / "_iqtree"
    work.mkdir(exist_ok=True)
    # species tree: fixed Fig S1 topology, ML branch lengths on the concatenation
    topo = work / "topology.nwk"
    topo.write_text(_TOPOLOGY + "\n")
    _run(["iqtree", "-s", str(concat), "-te", str(topo), "-m", "LG+G",
          "-pre", str(work / "sp"), "-redo", "-quiet"])
    (outdir / "tree_substitution.nwk").write_text((work / "sp.treefile").read_text())
    _datetree(outdir)
    # per-gene unconstrained ML trees (Fig S2 discordance) — left as substitution
    # phylograms on purpose: these exist to reproduce the genealogical discordance.
    gt = outdir / "gene_trees"
    gt.mkdir(exist_ok=True)
    for g, fa in per_gene.items():
        _run(["iqtree", "-s", str(fa), "-m", "LG+G", "-pre", str(work / g),
              "-redo", "-quiet"])
        (gt / f"{g}.nwk").write_text((work / f"{g}.treefile").read_text())


def _datetree(outdir: Path) -> None:
    """tree.nwk = the ML phylogram rooted on Aegithalidae and time-scaled to an
    ultrametric chronogram (ape::chronos, penalised likelihood, lambda=1,
    correlated rates; root age = 1, relative time).

    PhyloPhere's contrast-independence test (modified Dunn index) and its OU/BM
    Phylogenetic Shift Score both assume a TIME tree. On the raw ML phylogram the
    Sino-Himalayan ground tit *Pseudopodoces humilis* carries a terminal branch
    ~4x its sister's — that is molecular *rate*, not divergence *time*. It
    inflates the diameter of any contrast pair containing P. humilis, so the
    Pseudopodoces~Parus_minor contrast (truth pair for alphaA A34T) scores a
    sub-1 Dunn on the phylogram and is dropped. Dating fixes this: on the
    chronogram it clears the Dunn threshold. The unmodified phylogram is kept as
    tree_substitution.nwk for provenance / methods that want substitution lengths.
    """
    phylo = outdir / "tree_substitution.nwk"
    og = ", ".join(f'"{s}"' for s in _OUTGROUP)
    rscript = (
        "suppressPackageStartupMessages(library(ape)); "
        f'p <- read.tree("{phylo}"); '
        f"r <- root(p, outgroup = c({og}), resolve.root = TRUE); "
        "r <- multi2di(r); r$edge.length[r$edge.length <= 0] <- 1e-8; "
        'u <- chronos(r, lambda = 1, model = "correlated", quiet = TRUE); '
        'u <- ladderize(structure(unclass(u), class = "phylo")); '
        "stopifnot(is.rooted(u), is.ultrametric(u, tol = 1e-6)); "
        f'write.tree(u, "{outdir / "tree.nwk"}")'
    )
    _run(["Rscript", "-e", rscript])


def _write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, s in sorted(seqs.items()):
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")


def main() -> None:
    _fetch()
    prot = _load_proteins()

    aln_dir = HERE / "align"
    aln_dir.mkdir(exist_ok=True)
    aln = {g: _align(prot[g], g, aln_dir) for g in _GENES}
    _assert_landmarks(aln)
    for g in _GENES:
        _write_fasta(aln_dir / f"{g}.fasta", aln[g])

    # concatenated alignment for species-tree branch lengths
    names = sorted(FOCAL)
    concat = {n: "".join(aln[g][n] for g in _GENES) for n in names}
    cat_fa = HERE / "_concat.fasta"
    _write_fasta(cat_fa, concat)
    per_gene_fa = {g: aln_dir / f"{g}.fasta" for g in _GENES}
    _iqtree(cat_fa, per_gene_fa, HERE)
    cat_fa.unlink()
    _rmtree(HERE / "_iqtree")

    # traits.  `altitude` is a 0/1 fg/bg code (1 = high, cut at _HIGH_CUT m) so
    # the pipeline's ordinal path picks it up directly; `altitude_hl` keeps the
    # human-readable H/L for eyeballing.
    mids = {n: (lo + hi) // 2 for n, (lo, hi, _, _) in FOCAL.items()}
    with (HERE / "my_traits.tsv").open("w") as fh:
        fh.write("species\telev_mid\taltitude\taltitude_hl\tfamily\n")
        for n in names:
            hi = 1 if mids[n] >= _HIGH_CUT else 0
            fh.write(f"{n}\t{mids[n]}\t{hi}\t{'H' if hi else 'L'}\t{FOCAL[n][3]}\n")
    with (HERE / "phenotype.tsv").open("w") as fh:
        for n in names:
            fh.write(f"{n}\t{1 if mids[n] >= _HIGH_CUT else 0}\n")

    # synthetic support files (mirrors the PEPC fixture)
    (HERE / "ali_sp_names.txt").write_text("\n".join(names) + "\n")
    with (HERE / "gene_ensembl.tsv").open("w") as fh:
        fh.write("gene\tchr\tstart\tend\tstrand\tlength\thuman_protein_id\n")
        for i, g in enumerate(_GENES):
            L = _MATURE_LEN[g]
            fh.write(f"{g}\tchr{i+1}\t100000\t{100000 + L*3}\t+\t{L}\t.\n")
    with (HERE / "taxid.tsv").open("w") as fh:
        fh.write("tax_id\tspecies\tfamily\trank\tname_class\n")
        for i, n in enumerate(names):
            fh.write(f"{9200001+i}\t{n}\t{FOCAL[n][3]}\tspecies\tscientific name\n")

    nH = sum(1 for n in names if mids[n] >= _HIGH_CUT)
    print(f"\nHb-altitude fixture: {len(names)} tips ({nH} high / {len(names)-nH} low)"
          f" x 3 globin genes  (HBA {len(aln['HBA']['Parus_minor'])}, "
          f"HBD {len(aln['HBD']['Parus_minor'])}, HBB {len(aln['HBB']['Parus_minor'])} cols)")


if __name__ == "__main__":
    sys.exit(main())
