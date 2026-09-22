"""Build a codon-level CDS alignment for the PEPC fixture, for use as
ortholog_characterizator's --cds_dir (quality -> translation/BMGE ->
phylogeny -> positive_selection/FUBAR), alongside (not replacing)
build.py's own protein-level align/PEPC.fasta used directly by the
PhyloPhere GUI templates.

Real per-accession GenBank CDS (besnard2009/ppc1_genbank.gb, one
CDS `join()` feature per accession, intron-free once extracted via
Bio.SeqFeature.extract) mapped to each of align/PEPC.fasta's 78 tips 1:1
(code_accession.py -- single-accession species from Table
S1 directly; the 12 multi-accession species disambiguated by best-offset
local protein identity, see that file's docstring). Every one of the 78
tips' back-translated codon row is checked here against the fixture's
existing amino-acid row and written to align_cds/verify.tsv -- 100%
coverage for all 78 (exact for most; a handful of accessions cover a
longer or differently-registered stretch of the gene than the fixture's
455-aa fragment, in which case the fixture's row is an exact substring of
the accession's own translation -- see the coverage check below).

`Abildgaar` has no real accession (see README's reverse-check section) so
it is dropped from the fixture entirely by build.py, before this script
ever runs -- align/PEPC.fasta already has only 78 tips.

Same back-translation technique as ../hb/build_cds.py: walk the existing,
already-validated align/PEPC.fasta amino-acid alignment column by column,
emit one codon per non-gap column (consuming the real CDS 3 nt at a time)
and "---" per gap column, so the codon alignment shares the exact same
970-column (2910 nt) coordinate system as the trusted fixture already in
use, rather than risking a second independent alignment drifting from it.

Output: align_cds/PEPC.fasta (78 tips x 2910 nt columns).
"""

from __future__ import annotations

import sys
from pathlib import Path

import difflib

from Bio import Seq, SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
HERE = SCRIPT_DIR.parent  # fixture root (this script lives in scripts/)
RAW_CDS_DIR = HERE / "besnard2009"

sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(HERE.parent))  # validation/tier1/input -- shared _fixture_lib
from code_accession import CODE_ACCESSION  # noqa: E402
from build import _species_names  # noqa: E402
from _fixture_lib import read_fasta, write_fasta  # noqa: E402

# align/PEPC.fasta's tips are real species names (build.py's own rename, from
# the same CODE_ACCESSION table by way of GenBank's ORGANISM annotation) --
# not the original PCOC/ConDor short codes CODE_ACCESSION is keyed by. Derive
# the same code->name map build.py used, then re-key by name, so a lookup by
# tip name here always agrees with what's actually in align/PEPC.fasta.
_NAME_MAP, _ = _species_names(set(CODE_ACCESSION))
CODE_ACCESSION = {_NAME_MAP[code]: acc for code, acc in CODE_ACCESSION.items()
                  if code in _NAME_MAP}


def _load_cds_by_accession() -> dict[str, str]:
    out: dict[str, str] = {}
    for rec in SeqIO.parse(RAW_CDS_DIR / "ppc1_genbank.gb", "genbank"):
        acc = rec.id.split(".")[0]
        cds = [f for f in rec.features if f.type == "CDS"][0]
        out[acc] = str(cds.extract(rec.seq)).upper()
    return out


def _back_translate_into_alignment(nt_cds: str, aa_aligned: str, codon_start: int,
                                    translated: str, expected: str) -> str:
    """Walk aa_aligned column by column; emit one codon per non-gap column
    (consuming nt_cds 3 at a time, from codon_start) and '---' per gap column.

    `translated` (this accession's own conceptual translation) and
    `expected` (the fixture's ungapped row) are usually identical, in which
    case codon i simply belongs to non-gap column i. The one exception
    (Chrysithr/FM208000, see the docstring above) has one extra residue in
    `translated` that isn't in the fixture -- a small real allelic
    difference from whatever exact individual PCOC's own alignment used.
    difflib locates any such extra/missing residues so the matching codons
    still land in the right columns instead of a naive frame-shift.
    """
    nt_cds = nt_cds[codon_start - 1:]
    usable_len = (len(nt_cds) // 3) * 3
    codons = [nt_cds[i:i + 3] for i in range(0, usable_len, 3)]
    if translated == expected:
        codon_of_expected_pos = list(range(len(expected)))
    else:
        codon_of_expected_pos = [None] * len(expected)
        for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, translated, expected, autojunk=False).get_opcodes():
            if tag == "equal":
                for k in range(i2 - i1):
                    codon_of_expected_pos[j1 + k] = i1 + k
        if any(c is None for c in codon_of_expected_pos):
            raise SystemExit("could not align translated CDS to fixture row (unmapped residues)")
    out = []
    ei = 0  # index into `expected` / codon_of_expected_pos
    for col in aa_aligned:
        if col in ("-", "X"):
            out.append("---")
        else:
            out.append(codons[codon_of_expected_pos[ei]])
            ei += 1
    return "".join(out)


def main() -> None:
    aln_path = HERE / "align" / "PEPC.fasta"
    if not aln_path.exists():
        raise SystemExit(f"missing {aln_path} -- run build.py first")

    aa_aligned = read_fasta(aln_path)
    cds_by_acc = _load_cds_by_accession()
    codon_start_by_acc = {
        rec.id.split(".")[0]: int(
            [f for f in rec.features if f.type == "CDS"][0].qualifiers.get("codon_start", ["1"])[0]
        )
        for rec in SeqIO.parse(RAW_CDS_DIR / "ppc1_genbank.gb", "genbank")
    }

    out_dir = HERE / "align_cds"
    out_dir.mkdir(exist_ok=True)

    codon_aln: dict[str, str] = {}
    skipped = []
    verify_rows = ["code\taccession\tcoverage\texact"]
    for name, aa_row in aa_aligned.items():
        acc = CODE_ACCESSION.get(name)
        if acc is None:
            skipped.append(name)
            continue
        nt = cds_by_acc[acc]
        codon_start = codon_start_by_acc[acc]
        usable = nt[codon_start - 1:]
        usable = usable[: (len(usable) // 3) * 3]
        translated = str(Seq.Seq(usable).translate())
        expected = aa_row.replace("-", "").replace("X", "")
        # Exact match for most tips. A few accessions cover a longer or
        # differently-registered stretch of the gene than the fixture's own
        # 455-aa fragment (e.g. Cyp_era1/FM208065, or Ele_vivA/AB085948
        # which is a complete mRNA) -- .ratio() penalises those for the
        # length difference even when every fixture residue is an exact
        # substring match, so check *coverage of `expected`* by "equal"
        # opcodes instead (this also covers Chrysithr/FM208000, which has
        # one extra residue relative to whatever exact allele/consensus
        # PCOC's own alignment used).
        sm = difflib.SequenceMatcher(None, translated, expected, autojunk=False)
        covered = sum(j2 - j1 for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag == "equal")
        coverage = covered / len(expected)
        if coverage < 0.95:
            raise SystemExit(
                f"{name} ({acc}): nucleotide-consensus translation != "
                f"protein-alignment fixture (coverage={coverage:.1%})\n"
                f"  nt->aa:  {translated}\n  fixture: {expected}"
            )
        if translated != expected:
            print(f"  NOTE: {name} ({acc}) matches at {coverage:.1%} coverage, not "
                  "identical/co-terminal (see build_cds.py docstring)")
        verify_rows.append(f"{name}\t{acc}\t{coverage:.4f}\t{translated == expected}")
        codon_aln[name] = _back_translate_into_alignment(nt, aa_row, codon_start, translated, expected)

    write_fasta(out_dir / "PEPC.fasta", codon_aln)
    (out_dir / "verify.tsv").write_text("\n".join(verify_rows) + "\n")
    ncol = len(next(iter(codon_aln.values())))
    print(f"PEPC: {len(codon_aln)} tips x {ncol} nt columns ({ncol // 3} codons)")
    if skipped:
        print(f"Skipped (no real accession): {', '.join(skipped)}")
    print(f"\nCDS alignment written to {out_dir}")
    print(f"Verification report written to {out_dir / 'verify.tsv'}")


if __name__ == "__main__":
    main()
