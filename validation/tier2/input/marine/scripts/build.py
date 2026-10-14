#!/usr/bin/env python3
"""Build the Tier 2 fixture: marine lifestyle in 58 placental mammals, UCSC 100-way (hg19) alignments.

Inputs (fetched by hand, see ../README.md):
    --ucsc-dir   knownCanonical.exonAA.fa.gz   exon-wise amino acid alignments, 100 species, one canonical
                                               transcript per gene
                 knownCanonical.txt.gz         canonical transcripts of the current hg19 gene set (coordinates)
                 kgXref.txt.gz                 transcript -> gene symbol and UniProt accession
                 kgAlias.txt.gz                transcript -> aliases, which carry the symbols of the RER set
                                               that the current gene set has renamed
    --rds        UCSC62mammaltrees.rds         RERconverge trees object: 19,149 gene trees, master tree

Outputs, in the fixture directory (the alignments go to --align-out):
    tree.nwk, tree_substitution.nwk, gene_trees.nwk, rer_genes.tsv   (rds_to_inputs.R)
    <align-out>/<GENE>.fasta   stitched exons of the 58 placentals; tip names are Genus_species; columns that are
                               gaps in all 58 species are removed; a final stop column is removed
    my_traits.tsv              species, marine, control, order
    prune/marine_tips.txt      the marine tips, pruned in the terrestrial-control run
    ali_sp_names.txt, taxid.tsv, gene_ensembl.tsv, unmapped_genes.txt

A UCSC transcript is assigned to a gene by coordinates: the identifiers of the 2016 alignments (uc0xxxxx) are
not in the current gene tables, so the hg19 exons of each transcript are located inside the canonical
transcripts of knownCanonical.txt.gz, and the symbol comes from kgXref.
"""
from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

FIXTURE = Path(__file__).resolve().parents[1]
VALIDATION = FIXTURE.parents[2]
REPO = VALIDATION.parent
sys.path.insert(0, str(VALIDATION / "tier1" / "input"))
sys.path.insert(0, str(REPO / "bin"))
from _fixture_lib import write_fasta  # noqa: E402

# UCSC assembly code in the 100-way alignment -> (species, order). The RER object names the walrus odoRosDi.
SPECIES = {
    "hg19": ("Homo_sapiens", "Primates"), "panTro4": ("Pan_troglodytes", "Primates"),
    "gorGor3": ("Gorilla_gorilla", "Primates"), "ponAbe2": ("Pongo_abelii", "Primates"),
    "nomLeu3": ("Nomascus_leucogenys", "Primates"), "rheMac3": ("Macaca_mulatta", "Primates"),
    "macFas5": ("Macaca_fascicularis", "Primates"), "papHam1": ("Papio_hamadryas", "Primates"),
    "chlSab1": ("Chlorocebus_sabaeus", "Primates"), "calJac3": ("Callithrix_jacchus", "Primates"),
    "saiBol1": ("Saimiri_boliviensis", "Primates"), "otoGar3": ("Otolemur_garnettii", "Primates"),
    # NCBI has no species Tupaia chinensis; the Chinese tree shrew is the subspecies T. belangeri chinensis
    "tupChi1": ("Tupaia_belangeri", "Scandentia"),
    "speTri2": ("Ictidomys_tridecemlineatus", "Rodentia"), "jacJac1": ("Jaculus_jaculus", "Rodentia"),
    "micOch1": ("Microtus_ochrogaster", "Rodentia"), "criGri1": ("Cricetulus_griseus", "Rodentia"),
    "mesAur1": ("Mesocricetus_auratus", "Rodentia"), "mm10": ("Mus_musculus", "Rodentia"),
    "rn5": ("Rattus_norvegicus", "Rodentia"), "hetGla2": ("Heterocephalus_glaber", "Rodentia"),
    "cavPor3": ("Cavia_porcellus", "Rodentia"), "chiLan1": ("Chinchilla_lanigera", "Rodentia"),
    "octDeg1": ("Octodon_degus", "Rodentia"),
    "oryCun2": ("Oryctolagus_cuniculus", "Lagomorpha"), "ochPri3": ("Ochotona_princeps", "Lagomorpha"),
    "susScr3": ("Sus_scrofa", "Artiodactyla"), "vicPac2": ("Vicugna_pacos", "Artiodactyla"),
    "camFer1": ("Camelus_ferus", "Artiodactyla"), "turTru2": ("Tursiops_truncatus", "Artiodactyla"),
    "orcOrc1": ("Orcinus_orca", "Artiodactyla"), "panHod1": ("Pantholops_hodgsonii", "Artiodactyla"),
    "bosTau7": ("Bos_taurus", "Artiodactyla"), "oviAri3": ("Ovis_aries", "Artiodactyla"),
    "capHir1": ("Capra_hircus", "Artiodactyla"),
    "equCab2": ("Equus_caballus", "Perissodactyla"), "cerSim1": ("Ceratotherium_simum", "Perissodactyla"),
    "felCat5": ("Felis_catus", "Carnivora"), "canFam3": ("Canis_lupus", "Carnivora"),
    "musFur1": ("Mustela_putorius", "Carnivora"), "ailMel1": ("Ailuropoda_melanoleuca", "Carnivora"),
    "odoRosDiv1": ("Odobenus_rosmarus", "Carnivora"), "lepWed1": ("Leptonychotes_weddellii", "Carnivora"),
    "pteAle1": ("Pteropus_alecto", "Chiroptera"), "pteVam1": ("Pteropus_vampyrus", "Chiroptera"),
    "myoDav1": ("Myotis_davidii", "Chiroptera"), "myoLuc2": ("Myotis_lucifugus", "Chiroptera"),
    "eptFus1": ("Eptesicus_fuscus", "Chiroptera"),
    "eriEur2": ("Erinaceus_europaeus", "Eulipotyphla"), "sorAra2": ("Sorex_araneus", "Eulipotyphla"),
    "conCri1": ("Condylura_cristata", "Eulipotyphla"),
    "loxAfr3": ("Loxodonta_africana", "Proboscidea"), "eleEdw1": ("Elephantulus_edwardii", "Macroscelidea"),
    "triMan1": ("Trichechus_manatus", "Sirenia"), "chrAsi1": ("Chrysochloris_asiatica", "Afrosoricida"),
    "echTel2": ("Echinops_telfairi", "Afrosoricida"), "oryAfe1": ("Orycteropus_afer", "Tubulidentata"),
    "dasNov3": ("Dasypus_novemcinctus", "Cingulata"),
}
RER_CODE = {"odoRosDiv1": "odoRosDi"}  # tip label in the RER object when it differs from the alignment's
MARINE = ["turTru2", "orcOrc1", "odoRosDiv1", "lepWed1", "triMan1"]
CONTROL = ["oryAfe1", "vicPac2", "camFer1", "myoLuc2", "myoDav1"]  # terrestrial control foreground (Chikina et al. 2016)
MIN_SPECIES = 15  # as rer_minsp: genes with fewer species are not written
VALID = set("ACDEFGHIKLMNPQRSTVWY-")
STOP = ("Z", "*")
MIN_COVER = 0.5  # minimum share of a UCSC transcript's span covered by the canonical transcript it is assigned to


def iter_blocks(path: Path, species: set[str], wanted_tx: set[str] | None = None):
    """Yield (transcript, [(species, exon, coordinates, sequence)]) for consecutive records of one transcript.

    Only records of `species` are kept; with `wanted_tx`, transcripts outside it are skipped without parsing.
    """
    cur, recs, hdr, seq = None, [], None, []

    def flush():
        if hdr is not None:
            name, *rest = hdr.split()
            tx, sp, ex, _ = name[1:].rsplit("_", 3)
            if sp in species and (wanted_tx is None or tx in wanted_tx):
                recs.append((sp, int(ex), rest[3] if len(rest) > 3 else "", "".join(seq)))

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                hdr, seq = line.rstrip("\n"), []
                tx = line.split()[0][1:].rsplit("_", 3)[0]
                if cur is not None and tx != cur:
                    if recs:
                        yield cur, recs
                    recs = []
                cur = tx
            else:
                seq.append(line.strip())
        flush()
        if cur is not None and recs:
            yield cur, recs


def parse_coord(c: str):
    """'chr3:129247577-129247937+' -> (chrom, start, end, strand)."""
    m = re.match(r"(\w+):(\d+)-(\d+)([+-])", c)
    return (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)) if m else None


def load_canonical(ucsc: Path):
    """Canonical transcripts: bins of 100 kb -> intervals, ENST -> (symbol, UniProt), ENST -> (chrom, start, end)."""
    xref = {}
    with gzip.open(ucsc / "kgXref.txt.gz", "rt") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            xref[f[0]] = (f[4], f[2])
    bins, span, n = defaultdict(list), {}, 0
    with gzip.open(ucsc / "knownCanonical.txt.gz", "rt") as fh:
        for ln in fh:
            chrom, s, e, _, enst, _ = ln.rstrip("\n").split("\t")[:6]
            if enst not in xref or not xref[enst][0]:
                continue
            s, e, n = int(s), int(e), n + 1
            span[enst] = (chrom, s, e)
            for b in range(s // 100000, e // 100000 + 1):
                bins[(chrom, b)].append((s, e, enst))
    print(f"canonical transcripts with a symbol: {n}")
    return bins, xref, span


def load_aliases(ucsc: Path, wanted: set[str], xref: dict) -> dict[str, set[str]]:
    """Old symbol -> set of current gene symbols, from kgAlias (alias of any isoform) and kgXref (isoform symbol)."""
    out = defaultdict(set)
    with gzip.open(ucsc / "kgAlias.txt.gz", "rt") as fh:
        for ln in fh:
            enst, alias = ln.rstrip("\n").split("\t")[:2]
            if alias in wanted and enst in xref and xref[enst][0]:
                out[alias].add(xref[enst][0])
    return out


def assign_genes(ucsc: Path, rer_genes: set[str]):
    """Pass 1: map every UCSC transcript to a gene symbol by exon overlap, then every RER gene to a transcript.

    A UCSC transcript goes to the canonical transcript that covers the largest share of its coding exon bases (the
    smallest one on a tie, so a nested gene wins over its host). Full containment is too strict (the 2016 transcript can
    extend past the current one by a few bases or a UTR), the span of the exons is not reliable (a transcript can have
    exons in a paralog) and the span of the canonical transcript includes UTRs, so the share is taken over the exons.
    Canonical transcripts of named genes, then those with a UniProt accession, are preferred: unnamed ENSG genes
    (readthroughs, lncRNAs) overlap named ones.

    A RER gene is matched by its symbol; when the current gene set does not know that symbol, by the current
    symbol that kgAlias gives it, provided there is exactly one and no other RER gene already holds it.
    Returns ({RER gene: (human protein length, transcript, chrom, start, end, strand, UniProt)}, {RER gene: current
    symbol} for the genes matched by alias, [ambiguous RER genes]).
    """
    bins, xref, span = load_canonical(ucsc)
    by_symbol = {}  # current symbol -> (cover, human protein length, transcript, chrom, start, end, strand, UniProt)
    for tx, recs in iter_blocks(ucsc / "knownCanonical.exonAA.fa.gz", {"hg19"}):
        coords = [c for c in (parse_coord(r[2]) for r in recs if r[2]) if c]
        if not coords:
            continue
        total = sum(e - s + 1 for _, s, e, _ in coords)
        overlap = defaultdict(int)
        for chrom, s, e, _ in coords:
            seen = set()
            for bi in range(s // 100000, e // 100000 + 1):
                for cs, ce, enst in bins.get((chrom, bi), ()):
                    if enst not in seen:
                        seen.add(enst)
                        overlap[enst] += max(0, min(e, ce) - max(s, cs) + 1)
        cands = [(round(bp / total, 2), enst) for enst, bp in overlap.items() if bp / total >= MIN_COVER]
        cands = [c for c in cands if not xref[c[1]][0].startswith("ENSG")] or cands  # named genes first
        cands = [c for c in cands if xref[c[1]][1]] or cands  # then those with a UniProt accession
        if not cands:
            continue
        cover, enst = max(cands, key=lambda c: (c[0], -(span[c[1]][2] - span[c[1]][1])))
        chrom, cs, ce = span[enst]
        symbol, uniprot = xref[enst]
        hum = "".join(r[3] for r in sorted(recs, key=lambda r: r[1]))
        plen = len(hum.replace("-", "")) - hum.endswith(STOP)  # residues of the human protein, stop excluded
        if symbol not in by_symbol or (cover, plen) > by_symbol[symbol][:2]:
            by_symbol[symbol] = (cover, plen, tx, chrom, cs, ce, coords[0][3], uniprot)

    by_symbol = {k: v[1:] for k, v in by_symbol.items()}  # drop the cover share
    best = {g: by_symbol[g] for g in rer_genes if g in by_symbol}
    held = set(best)  # current symbols already taken by a RER gene of the same name
    unknown = rer_genes - set(best)
    aliases = load_aliases(ucsc, unknown, xref)
    via_alias, ambiguous = {}, []
    for g in sorted(unknown):
        current = {c for c in aliases.get(g, ()) if c in by_symbol}
        if len(current) == 1:
            (c,) = current
            if c not in held:
                best[g], via_alias[g] = by_symbol[c], c
                held.add(c)
                continue
        if current:
            ambiguous.append(g)
    return best, via_alias, ambiguous


def clean_alignment(rows: dict[str, str]) -> dict[str, str]:
    """Drop a final stop column, map non-standard residues to X, drop species and columns without data."""
    human = rows["Homo_sapiens"]
    if human.endswith(STOP):  # UCSC writes the stop codon as Z
        rows = {k: v[:-1] for k, v in rows.items()}
    rows = {k: "".join(ch if ch in VALID else "X" for ch in v) for k, v in rows.items()}
    rows = {k: v for k, v in rows.items() if v.strip("-X")}
    if "Homo_sapiens" not in rows:
        return {}
    keep = [i for i in range(len(rows["Homo_sapiens"])) if any(v[i] != "-" for v in rows.values())]
    return {k: "".join(v[i] for i in keep) for k, v in rows.items()}


def write_alignments(ucsc: Path, best: dict, out: Path) -> set[str]:
    """Pass 2: stitch the exons of the chosen transcripts and write one FASTA per gene."""
    out.mkdir(parents=True, exist_ok=True)
    by_tx = {v[1]: sym for sym, v in best.items()}
    written = set()
    for tx, recs in iter_blocks(ucsc / "knownCanonical.exonAA.fa.gz", set(SPECIES), set(by_tx)):
        parts = defaultdict(list)
        for sp, ex, _, seq in sorted(recs, key=lambda r: (r[0], r[1])):
            parts[sp].append(seq)
        rows = {SPECIES[sp][0]: "".join(p) for sp, p in parts.items()}
        assert len({len(v) for v in rows.values()}) == 1, f"unequal exon lengths in {tx}"
        rows = clean_alignment(rows)
        if len(rows) < MIN_SPECIES:
            continue
        write_fasta(out / (re.sub(r"[^A-Za-z0-9_.-]", "_", by_tx[tx]) + ".fasta"), rows)
        written.add(by_tx[tx])
    return written


def write_tables(best: dict, align_n: int) -> None:
    from generate_taxid_map import resolve_taxids  # live NCBI first, as the pipeline's own tax_id generation

    species = sorted(v[0] for v in SPECIES.values())
    order = {v[0]: v[1] for v in SPECIES.values()}
    marine = {SPECIES[c][0] for c in MARINE}
    control = {SPECIES[c][0] for c in CONTROL}
    with (FIXTURE / "my_traits.tsv").open("w") as fh:
        fh.write("species\tmarine\tcontrol\torder\n")
        for sp in species:
            fh.write(f"{sp}\t{int(sp in marine)}\t{int(sp in control)}\t{order[sp]}\n")
    (FIXTURE / "prune").mkdir(exist_ok=True)
    (FIXTURE / "prune" / "marine_tips.txt").write_text("\n".join(sorted(marine)) + "\n")
    (FIXTURE / "ali_sp_names.txt").write_text("\n".join(species) + "\n")
    resolved, unresolved = resolve_taxids(species)
    if unresolved:
        raise SystemExit(f"NCBI taxonomy lookup failed for: {unresolved}")
    with (FIXTURE / "taxid.tsv").open("w") as fh:
        fh.write("tax_id\tspecies\torder\trank\tname_class\n")
        for sp in species:
            fh.write(f"{resolved[sp]}\t{sp}\t{order[sp]}\tspecies\tscientific name\n")
    with (FIXTURE / "gene_ensembl.tsv").open("w") as fh:
        fh.write("gene\tchr\tstart\tend\tstrand\tlength\thuman_protein_id\n")
        for sym in sorted(best):
            plen, _, chrom, cs, ce, strand, uniprot = best[sym]
            fh.write(f"{sym}\t{chrom}\t{cs}\t{ce}\t{strand}\t{plen}\t{uniprot or 'NA'}\n")
    print(f"tables written: {len(species)} species, {len(best)} genes mapped, {align_n} alignments")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ucsc-dir", type=Path, required=True)
    ap.add_argument("--rds", type=Path, required=True)
    ap.add_argument("--align-out", type=Path, required=True)
    ap.add_argument("--stage", choices=["all", "trees", "align"], default="all")
    a = ap.parse_args()

    if a.stage in ("all", "trees"):
        map_file = FIXTURE / "species_map.tsv"
        map_file.write_text("rer_code\tspecies\n" + "".join(
            f"{RER_CODE.get(c, c)}\t{v[0]}\n" for c, v in sorted(SPECIES.items())))
        subprocess.run(["Rscript", "--vanilla", str(Path(__file__).with_name("rds_to_inputs.R")),
                        str(a.rds), str(map_file), str(FIXTURE)], check=True)
    if a.stage in ("all", "align"):
        rer_genes = {ln.split("\t")[0] for ln in (FIXTURE / "rer_genes.tsv").read_text().splitlines()
                     if int(ln.split("\t")[1]) >= MIN_SPECIES}
        best, via_alias, ambiguous = assign_genes(a.ucsc_dir, rer_genes)
        missing = sorted(rer_genes - set(best))
        (FIXTURE / "unmapped_genes.txt").write_text("\n".join(missing) + "\n")
        (FIXTURE / "gene_aliases.tsv").write_text("rer_gene\tcurrent_symbol\n" + "".join(
            f"{g}\t{c}\n" for g, c in sorted(via_alias.items())))
        print(f"RER genes with at least {MIN_SPECIES} placentals: {len(rer_genes)}; mapped: {len(best)} "
              f"({len(via_alias)} through kgAlias); unmapped: {len(missing)} (unmapped_genes.txt), "
              f"of which {len(ambiguous)} with an ambiguous or already held alias")
        written = write_alignments(a.ucsc_dir, best, a.align_out)
        write_tables({g: v for g, v in best.items() if g in written}, len(written))


if __name__ == "__main__":
    main()
