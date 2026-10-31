"""Validated protein alignment IO shared by selection preparation helpers."""
from pathlib import Path
from Bio import AlignIO, SeqIO


def read_alignment(path):
    path = Path(path)
    # Detect FASTA from content, including extensionless inputs.
    with path.open() as fh:
        first = next((line.strip() for line in fh if line.strip()), "")
    if first.startswith(">"):
        with path.open() as fh:
            records = list(SeqIO.parse(fh, "fasta"))
    else:
        records = None
        for fmt in ("phylip-relaxed", "phylip-sequential", "phylip"):
            try:
                records = list(AlignIO.read(path, fmt))
                break
            except ValueError:
                continue
        if records is None:
            raise ValueError(f"Cannot parse alignment as FASTA or PHYLIP: {path}")
    if not records:
        raise ValueError(f"Empty alignment: {path}")
    seen, lengths = set(), set()
    for rec in records:
        if not rec.id or rec.id in seen:
            raise ValueError(f"Empty or duplicate taxon {rec.id!r}: {path}")
        seen.add(rec.id)
        seq = str(rec.seq)
        if not seq or any(c.isspace() for c in seq):
            raise ValueError(f"Empty sequence or embedded whitespace for {rec.id}: {path}")
        if set(seq.upper()) - set("ABCDEFGHIJKLMNOPQRSTUVWXYZ-?*."):
            raise ValueError(f"Invalid protein sequence for {rec.id}: {path}")
        lengths.add(len(seq))
        rec.description = ""
    if len(lengths) != 1:
        raise ValueError(f"Unequal alignment sequence lengths: {path}")
    return records


def informative(record):
    # Preserve the historical PHYLIP rule, now applied to FASTA too.
    return bool(set(str(record.seq).upper()) - set("-?X"))


def tree_taxa(path):
    import dendropy
    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    labels = [leaf.taxon.label if leaf.taxon else None for leaf in tree.leaf_node_iter()]
    if not labels or any(not x for x in labels) or len(labels) != len(set(labels)):
        raise ValueError(f"Empty, unnamed or duplicate tree tips: {path}")
    return set(labels)
