"""The files the rest of the pipeline reads about the observed labeling, written from the b_0 slice of the core.

A perm-replay batch leaves, per gene, `<id>.b0.discovery.tsv` (the discovery.tab rows of b_0, only when it has a
hit), `<id>.b0.background` (`<gene>\\t<positions tested | NULL>`) and, from the ASR scoring, `<gene>.master.csv.gz`
(the gene's rows of the master table). This module turns the files of all the batches into

  discovery.tab, background.output, background_genes.output,
  meta_caas/global_meta_caas.tsv and meta_caas/<group>_meta_caas.tsv, caas_convergence_master.csv.

Every table is ordered by gene (byte order, as `LC_ALL=C sort` does), so it does not depend on the order in which
the batches finished; rows of one gene keep the order they were written in. The CAAS ids of the meta tables are
the content ids of `core.meta`, the same ones `tag_support` of the master holds. The meta tables have the columns
and the text format of the tables the pattern-annotation report writes (readr conventions: an empty cell is `NA`,
logicals are `TRUE`/`FALSE`).
"""
import csv
import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from src.core.master import write_master_csv
from src.core.meta import row_id

DISCOVERY_SUFFIX = ".b0.discovery.tsv"
BACKGROUND_SUFFIX = ".b0.background"
MASTER_SUFFIX = ".master.csv.gz"

# the columns of discovery.tab, with the two conserved-pair columns: what an empty discovery.tab carries
EMPTY_DISCOVERY_HEADER = ["gene", "mode", "caap_group", "trait", "position", "caas", "amino_encoded", "pattern",
                          "ffgn", "fbgn", "gfg", "gbg", "mfg", "mbg", "ffg", "fbg", "ms", "is_conserved_meta", "conserved_pair"]

META_COLUMNS = ["tag", "GenePos", "Gene", "Position", "pattern", "caap_group", "is_conserved_meta", "conserved_pair",
                "caas", "amino_encoded"]
_HYP = re.compile(r"H[0-9]+")
_HYP_LAST = re.compile(r".*(H[0-9]+).*")


def batch_files(dirs: Iterable, suffix: str) -> List[Path]:
    """The files of the batch directories that end in `suffix`, in file-name order."""
    return sorted((f for d in dirs for f in Path(d).glob(f"*{suffix}")), key=lambda f: (f.name, str(f)))


def write_discovery(files: Sequence[Path], out: Path) -> int:
    """discovery.tab: one header (the files', or the empty-table header), then every file's rows by gene.
    Returns the number of rows."""
    header: Optional[str] = None
    rows: List[tuple] = []
    for f in files:
        with open(f) as fh:
            head = fh.readline()
            if header is None:
                header = head
            for line in fh:
                rows.append((line.split("\t", 1)[0], line))
    rows.sort(key=lambda r: r[0])  # stable
    with open(out, "w") as fh:
        fh.write(header if header is not None else "\t".join(EMPTY_DISCOVERY_HEADER) + "\n")
        fh.writelines(line for _, line in rows)
    return len(rows)


def write_background(files: Sequence[Path], out: Path, out_genes: Path) -> int:
    """background.output (`<gene>\\t<positions>` lines by gene; `Gene\\tPosition` alone when there is none) and
    background_genes.output (the genes with at least one tested position, sorted, once each)."""
    lines: List[tuple] = []
    for f in files:
        with open(f) as fh:
            lines += [(line.split("\t", 1)[0].rstrip("\n"), line if line.endswith("\n") else line + "\n") for line in fh if line.strip()]
    lines.sort(key=lambda r: r[0])
    with open(out, "w") as fh:
        fh.write("".join(line for _, line in lines) if lines else "Gene\tPosition\n")
    genes = set()
    for _, line in lines:
        cols = line.rstrip("\n").split("\t")
        if len(cols) >= 2 and cols[1] not in ("", "Position", "NULL"):
            genes.add(cols[0])
    with open(out_genes, "w") as fh:
        fh.writelines(g + "\n" for g in sorted(genes))
    return len(lines)


def _na(value: Optional[str]) -> str:
    return "NA" if value is None or value == "" else value


def write_meta(discovery: Path, out_dir: Path) -> Dict[str, int]:
    """global_meta_caas.tsv and <group>_meta_caas.tsv from discovery.tab, row by row in its order.
    Returns the number of rows of each file by group ('global' included)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    counts: Dict[str, int] = {}
    handles: Dict[str, object] = {}
    with open(discovery, newline="") as src:
        reader = csv.DictReader(src, delimiter="\t", quoting=csv.QUOTE_NONE)
        has_trait = "trait" in (reader.fieldnames or [])
        group_cols = META_COLUMNS + (["trait"] if has_trait else [])
        global_cols = group_cols + ["hyp_id"] if has_trait else META_COLUMNS + ["hyp_id"]
        glob = open(out_dir / "global_meta_caas.tsv", "w")
        glob.write("\t".join(global_cols) + "\n")
        counts["global"] = 0
        try:
            for row in reader:
                gene = row["gene"]
                group = row.get("caap_group") or "US"
                trait = row.get("trait")
                pair = row.get("conserved_pair") or ""
                fields = [row_id(gene, row), f"{gene}_{row['position']}", gene, row["position"], _na(row.get("pattern")), group,
                          "TRUE" if (row.get("is_conserved_meta") or "").strip().upper() == "TRUE" else "FALSE",
                          re.sub(r"^\d+:", "", pair) if pair else "NA", _na(row.get("caas")), _na(row.get("amino_encoded"))]
                if has_trait:
                    fields.append(_na(trait))
                    hyp = _HYP_LAST.sub(r"\1", trait) if trait and _HYP.search(trait) else trait
                    glob.write("\t".join(fields + [_na(hyp)]) + "\n")
                else:
                    glob.write("\t".join(fields + ["NA"]) + "\n")
                counts["global"] += 1
                if group not in handles:
                    handles[group] = open(out_dir / f"{group}_meta_caas.tsv", "w")
                    handles[group].write("\t".join(group_cols) + "\n")
                    counts[group] = 0
                handles[group].write("\t".join(fields) + "\n")
                counts[group] += 1
        finally:
            glob.close()
            for h in handles.values():
                h.close()
    return counts


def write_master(shards: Sequence[Path], out: Path, fallback_fields: Sequence[str]) -> int:
    """caas_convergence_master.csv: the shards' rows ordered by (gene, msa_pos), header from the first shard
    (`fallback_fields` when there is none). Returns the number of rows."""
    rows = []
    fields: Optional[List[str]] = None
    for f in shards:
        with gzip.open(f, "rt", newline="") as fh:
            reader = csv.DictReader(fh)
            fields = fields or list(reader.fieldnames or [])
            for row in reader:
                rows.append((row["gene"], int(row["msa_pos"]) if row["msa_pos"] else None, row))
    return write_master_csv(rows, out, fields or list(fallback_fields))
