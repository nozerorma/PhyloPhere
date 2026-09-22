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
    tree.nwk             PhyML tree rooted on the Chrysithrix outgroup and
                         time-scaled to an ultrametric chronogram (ape::chronos)
                         - PhyloPhere contrast selection assumes a time tree
    tree_substitution.nwk   the raw PhyML phylogram, kept for provenance
    my_traits.tsv        species <tab> c4 (1|0) <tab> family
    phenotype.tsv        species <tab> 1|0     (C4 foreground, outgroup dropped)
    ali_sp_names.txt / gene_ensembl.tsv / taxid.tsv    synthetic support files

Numbering note: the truth set (validation/truthsets/tier1/pepc_c4.sites.tsv) is
in maize P04711 coordinates, which equal alignment columns 1:1 here.
"""

from __future__ import annotations

import re
import sys
import urllib.request
from collections import Counter
from pathlib import Path

from code_accession import CODE_ACCESSION

HERE = Path(__file__).resolve().parent.parent  # fixture root (this script lives in scripts/)
REPO_ROOT = HERE.parents[3]  # .../PhyloPhere -- pepc -> input -> tier1 -> validation -> PhyloPhere
sys.path.insert(0, str(REPO_ROOT / "bin"))
sys.path.insert(0, str(HERE.parent))  # validation/tier1/input -- shared _fixture_lib
from generate_taxid_map import resolve_taxids  # noqa: E402
from _fixture_lib import date_tree, read_fasta, rename_tree_tips, write_fasta  # noqa: E402

RAW_DIR = HERE / "besnard2009"
RAW_ALN = RAW_DIR / "pepc.aa.coor_mays.fa"
RAW_TREE = RAW_DIR / "pepc.phyml_tree.txt"
RAW_GB = RAW_DIR / "ppc1_genbank.gb"
C4_FILE = RAW_DIR / "besnard2009_convergent_species.txt"
OG_FILE = RAW_DIR / "outgroup.txt"

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
    RAW_DIR.mkdir(exist_ok=True)
    for dst, name in _RAW.items():
        if dst.exists():
            continue
        print(f"fetch {name}")
        with urllib.request.urlopen(_BASE + name, timeout=60) as r:
            dst.write_bytes(r.read())


def _san(name: str) -> str:
    return name.replace(".", "_")


def _species_names(codes: set[str]) -> tuple[dict[str, str], dict[str, str]]:
    """Returns (tip_name, base_binomial), both keyed by code, for every code
    with a resolved GenBank accession in code_accession.py's CODE_ACCESSION
    (Abildgaar has none and is omitted here -- dropped from the fixture
    entirely, see README). The binomial is GenBank's own current ORGANISM
    annotation for that accession (i.e. present-day accepted taxonomy, not
    necessarily Besnard 2009's own label -- e.g. Baumea -> Machaerina
    articulata, Killinga -> Rhynchospora colorata). tip_name is base_binomial
    plus an accession suffix for the 12 species with 2-4 accessions (e.g.
    Cyperus_eragrostis_FM208065), so tips stay unique; base_binomial is the
    same for every tip of a shared species (needed for a single, shared
    NCBI tax_id lookup per species rather than per tip -- see _ncbi_tax_ids).
    """
    if not RAW_GB.exists():
        raise SystemExit(
            f"missing {RAW_GB} -- this is a one-off manual efetch of the "
            "Besnard 2009 ppc-1 GenBank records (see README's reverse-check "
            "section), not fetched automatically by this script."
        )
    from Bio import SeqIO

    org_by_acc = {
        rec.id.split(".")[0]: rec.annotations["organism"]
        for rec in SeqIO.parse(RAW_GB, "genbank")
    }
    code_org = {c: org_by_acc[CODE_ACCESSION[c]] for c in codes if CODE_ACCESSION.get(c)}
    counts = Counter(code_org.values())
    base_binomial = {code: _san(org.replace(" ", "_")) for code, org in code_org.items()}
    tip_name = {
        code: f"{base}_{CODE_ACCESSION[code]}" if counts[code_org[code]] > 1 else base
        for code, base in base_binomial.items()
    }
    assert len(set(tip_name.values())) == len(tip_name), "species-name collision after renaming"
    return tip_name, base_binomial


def _ncbi_tax_ids(binomials: set[str]) -> dict[str, int]:
    """species binomial (underscored) -> NCBI taxonomy id, via
    bin/generate_taxid_map.py's resolve_taxids -- the same live-NCBI-first,
    local-ete3-fallback resolution PhyloPhere's own tax_id auto-generation
    uses. That auto-generation resolves directly against tree TIP labels,
    though, which fails outright for every multi-accession tip here (e.g.
    "Cyperus eragrostis FM208065" is not a taxon -- resolve_taxids' own
    genus+species fallback handles that specific shape, but this fixture
    still resolves once per BASE species up front, shared across a
    multi-accession group's sibling tips, and supplies a ready-made
    --tax_id file rather than leaving it for PhyloPhere to generate).
    """
    resolved, unresolved = resolve_taxids(sorted(binomials))
    if unresolved:
        raise SystemExit(f"NCBI taxonomy lookup failed for: {unresolved}")
    return resolved


def main() -> None:
    _fetch()
    aln_by_code = {_san(k): v for k, v in read_fasta(RAW_ALN).items()}

    # Rename tips to real species binomials (GenBank's current ORGANISM
    # annotation for each accession -- see _species_names). Abildgaar has no
    # resolved accession and is dropped here, from every output, for the
    # first time putting the AA-alignment track and the CDS track (which has
    # always excluded it, see build_cds.py) in agreement.
    name_map, base_map = _species_names(set(aln_by_code))
    dropped = sorted(set(aln_by_code) - set(name_map))
    if dropped:
        print(f"Dropping (no real GenBank accession): {', '.join(dropped)}")
    aln = {name_map[k]: v for k, v in aln_by_code.items() if k in name_map}

    tree = RAW_TREE.read_text().strip()
    # sanitise tip names in the newick (they appear as "(<name>:" or ",<name>:")
    tree = re.sub(r"([(,])([A-Za-z0-9._-]+?):",
                  lambda m: f"{m.group(1)}{_san(m.group(2))}:", tree)

    c4_codes = {_san(x) for x in C4_FILE.read_text().split()}
    outgroup_codes = {_san(x) for x in OG_FILE.read_text().split()}
    assert c4_codes <= set(aln_by_code), c4_codes - set(aln_by_code)
    assert outgroup_codes <= set(aln_by_code), outgroup_codes - set(aln_by_code)
    assert not (c4_codes & set(dropped)), "a dropped tip is listed as C4 -- fixture is broken"
    assert not (outgroup_codes & set(dropped)), "the outgroup tip was dropped -- fixture is broken"
    c4 = {name_map[x] for x in c4_codes}
    outgroup = {name_map[x] for x in outgroup_codes}

    L = len(next(iter(aln.values())))
    assert L == 970, f"expected 970 maize-coordinate columns, got {L}"
    assert all(len(s) == L for s in aln.values())

    (HERE / "align").mkdir(exist_ok=True)
    write_fasta(HERE / "align" / "PEPC.fasta", aln)
    (HERE / "tree_substitution.nwk").write_text(tree + "\n")
    date_tree(HERE, outgroup_codes, dropped)
    # Final pass: swap the short PCOC/ConDor codes for real species names in
    # both committed tree files (date_tree only ever sees the original,
    # already-verified tip identifiers it was given for outgroup/drop_tips).
    for f in (HERE / "tree_substitution.nwk", HERE / "tree.nwk"):
        rename_tree_tips(f, name_map)

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
    # Real NCBI tax_ids, not a synthetic placeholder: PhyloPhere's own
    # auto-generation (bin/generate_taxid_map.py) does exact-name lookup
    # against tip labels directly and fails on every accession-suffixed
    # multi-accession tip, so this fixture supplies --tax_id itself, one
    # shared lookup per base species rather than per tip.
    tip_to_base = {name_map[code]: base_map[code] for code in name_map}
    tax_id_of_base = _ncbi_tax_ids(set(base_map.values()))
    with (HERE / "taxid.tsv").open("w") as fh:
        fh.write("tax_id\tspecies\tfamily\trank\tname_class\n")
        for sp in species:
            fh.write(f"{tax_id_of_base[tip_to_base[sp]]}\t{sp}\tCyperaceae\tspecies\tscientific name\n")

    n_c4 = sum(1 for sp in species if sp in c4)
    print(f"PEPC fixture: {len(species)} tips  "
          f"({n_c4} C4, {len(species) - n_c4 - len(outgroup)} C3, "
          f"{len(outgroup)} outgroup)  x {L} maize-coord columns")


if __name__ == "__main__":
    main()
