"""Build codon-level CDS alignments for the Hb fixture, for use as
ortholog_characterizator's --cds_dir (quality -> translation/BMGE -> phylogeny
-> positive_selection/FUBAR), alongside (not replacing) build.py's own
protein-level align/*.fasta used directly by the PhyloPhere GUI templates.

Same GenBank records as build.py (hb_genbank.gb, MG772099-MG772439), same
per-species-per-gene majority consensus, but taken at the nucleotide level
from each CDS feature's own coordinates (Bio.SeqFeature.extract) instead of
the record's annotated /translation. The codon alignment is then built by
back-translating each species' nucleotide consensus into the existing,
already-validated align/{GENE}.fasta amino-acid alignment's gap pattern
(one codon per non-gap AA column, "---" per gap column) -- this guarantees
the same topology/columns as the fixture already in use, rather than risking
a second independent alignment drifting from it.

Output: align_cds/{HBA,HBD,HBB}.fasta (codon nucleotide alignments, tips
matching align/*.fasta exactly).

Needs the fixture-build env (mafft, iqtree) only insofar as build.py must
have already been run once (align/*.fasta must exist).
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).parent
RAW_GB = HERE / "hb_genbank.gb"

sys.path.insert(0, str(HERE))
from build import FOCAL, _ORG2NAME, _GENES, _gene_of  # noqa: E402


def _consensus(seqs: list[str]) -> str:
    """Per-column majority over a species' allele set (ties -> first seen).
    Tolerates differing isolate lengths (partial-cds boundary calls differ
    per GenBank submission) the same way build.py's protein-level version
    does: majority over whichever isolates cover that position."""
    L = max(len(s) for s in seqs)
    out = []
    for i in range(L):
        col = [s[i] for s in seqs if i < len(s)]
        out.append(Counter(col).most_common(1)[0][0])
    return "".join(out)


def _load_cds() -> dict[str, dict[str, str]]:
    raw: dict[str, dict[str, list[str]]] = {g: defaultdict(list) for g in _GENES}
    for rec in SeqIO.parse(RAW_GB, "genbank"):
        org = rec.annotations.get("organism", "")
        name = _ORG2NAME.get(org)
        if name is None:
            continue
        g = _gene_of(rec.description)
        if g == "?":
            continue
        for f in rec.features:
            if f.type == "CDS" and "translation" in f.qualifiers:
                nt = str(f.extract(rec.seq)).upper()
                # Drop a trailing stop codon if present so length is a clean
                # multiple of 3 matching the annotated (stop-excluded) AA.
                aa_len = len(f.qualifiers["translation"][0])
                nt = nt[: aa_len * 3]
                raw[g][name].append(nt)
                break

    out: dict[str, dict[str, str]] = {}
    for g in _GENES:
        out[g] = {}
        for name in FOCAL:
            seqs = raw[g].get(name, [])
            if not seqs:
                raise SystemExit(f"{g}: no CDS for {name}")
            out[g][name] = _consensus(seqs)
    return out


def _back_translate_into_alignment(nt_cds: str, aa_aligned: str, strip_met: bool) -> str:
    """Walk aa_aligned column by column; emit one codon per non-gap column
    (consuming nt_cds 3 at a time) and '---' per gap column."""
    codons = [nt_cds[i : i + 3] for i in range(0, len(nt_cds), 3)]
    if strip_met:
        codons = codons[1:]
    out = []
    ci = 0
    for col in aa_aligned:
        if col == "-":
            out.append("---")
        else:
            out.append(codons[ci])
            ci += 1
    if ci != len(codons):
        raise SystemExit(
            f"codon/AA-column count mismatch: consumed {ci} of {len(codons)} codons"
        )
    return "".join(out)


def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def _write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, s in sorted(seqs.items()):
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")


def main() -> None:
    aln_dir = HERE / "align"
    for g in _GENES:
        if not (aln_dir / f"{g}.fasta").exists():
            raise SystemExit(f"missing {aln_dir / f'{g}.fasta'} -- run build.py first")

    cds = _load_cds()
    out_dir = HERE / "align_cds"
    out_dir.mkdir(exist_ok=True)

    for g in _GENES:
        aa_aligned = _read_fasta(aln_dir / f"{g}.fasta")
        strip_met = g in ("HBA", "HBB")
        codon_aln: dict[str, str] = {}
        for name, aa_row in aa_aligned.items():
            nt = cds[g][name]
            # Sanity check: nucleotide consensus translates to the same
            # ungapped AA sequence already used in the trusted fixture.
            translated = str(Seq(nt).translate(to_stop=True))
            if strip_met and translated[0] == "M":
                translated = translated[1:]
            expected = aa_row.replace("-", "")
            if translated != expected:
                raise SystemExit(
                    f"{g}/{name}: nucleotide-consensus translation != "
                    f"protein-consensus fixture\n  nt->aa:  {translated}\n"
                    f"  fixture: {expected}"
                )
            codon_aln[name] = _back_translate_into_alignment(nt, aa_row, strip_met)
        _write_fasta(out_dir / f"{g}.fasta", codon_aln)
        ncol = len(next(iter(codon_aln.values())))
        print(f"  {g}: {len(codon_aln)} species x {ncol} nt columns "
              f"({ncol // 3} codons)")

    print(f"\nCDS alignments written to {out_dir}")


if __name__ == "__main__":
    sys.exit(main())
