# contract.py — Write the observed-labeling tables (discovery, background, meta CAAS, master) from the b_0 batches.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/core/

"""
The files the rest of the pipeline reads about the observed labeling, written from the b_0 slice of the core.

Imported by: contract_main.py
Inputs: per-batch `<id>.b0.discovery.tsv`, `<id>.b0.background`, `<id>.b0.universe` and `<gene>.master.csv.gz` files
Outputs: discovery.tab, background.output, background_genes.output, meta_caas/meta_caas.tsv, caas_convergence_master.csv

A perm-replay batch leaves, per gene, `<id>.b0.discovery.tsv` (the discovery.tab rows of b_0, only when it has a
hit), `<id>.b0.background` (`<gene>\\t<positions tested | NULL>`) and, from the ASR scoring, `<gene>.master.csv.gz`
(the gene's rows of the master table). This module turns the files of all the batches into

  discovery.tab, background.output, background_genes.output,
  meta_caas/meta_caas.tsv, caas_convergence_master.csv.

Every table is ordered by gene (byte order, as `LC_ALL=C sort` does), so it does not depend on the order in which
the batches finished; rows of one gene keep the order they were written in. The CAAS ids of the meta tables are
the content ids of `core.meta`, the same ones `tag_support` of the master holds. The meta tables have the columns
and the text format of the tables the pattern-annotation report writes (readr conventions: an empty cell is `NA`,
logicals are `TRUE`/`FALSE`).
"""
import csv
import gzip
import hashlib
import math
import tempfile
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from src.core.master import write_master_csv
from src.core.meta import row_id

DISCOVERY_SUFFIX = ".b0.discovery.tsv"
BACKGROUND_SUFFIX = ".b0.background"
from src.core.universe import UNIVERSE_SUFFIX  # noqa: E402  (`*.b0.universe`: the positional universe shards)
MASTER_SUFFIX = ".master.csv.gz"

# Header of an empty discovery.tab: the columns of discovery.tab, including the two conserved-pair columns.
EMPTY_DISCOVERY_HEADER = ["gene", "trait", "position", "caas", "pattern",
                          "ffgn", "fbgn", "gfg", "gbg", "mfg", "mbg", "ffg", "fbg", "ms", "is_conserved_meta", "conserved_pair", "shared_residue_count"]

META_COLUMNS = ["tag", "GenePos", "Gene", "Position", "pattern", "is_conserved_meta", "conserved_pair",
                "caas", "shared_residue_count"]
_HYP = re.compile(r"H[0-9]+")
_HYP_LAST = re.compile(r".*(H[0-9]+).*")


def batch_files(dirs: Iterable, suffix: str) -> List[Path]:
    """The files of the batch directories that end in `suffix`, in file-name order."""
    return sorted((f for d in dirs for f in Path(d).glob(f"*{suffix}")), key=lambda f: (f.name, str(f)))


def write_discovery(files: Sequence[Path], out: Path) -> int:
    """discovery.tab: one header (the files', or the empty-table header), then every file's rows by gene.
    Returns the number of rows."""
    # SQLite spills the global sort to disk instead of holding the discovery in RAM.
    import sqlite3
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.discovery-merge-', dir=out.parent) as scratch:
        db = sqlite3.connect(str(Path(scratch) / 'sort.sqlite'))
        db.execute('CREATE TABLE rows (gene TEXT, line TEXT)')
        header = None
        n = 0
        for source in files:
            with source.open() as fh:
                columns = fh.readline().rstrip('\r\n')
                if header is not None and header != columns:
                    raise ValueError(f'Discovery shard schemas differ: {source}')
                header = columns
                if not {'gene','position','caas','pattern'} <= set(columns.split('\t')):
                    raise ValueError(f'Invalid discovery header: {source}')
                for line in fh:
                    if len(line.rstrip('\r\n').split('\t')) != len(columns.split('\t')):
                        raise ValueError(f'Discovery column count differs: {source}')
                    db.execute('INSERT INTO rows VALUES (?,?)', (line.split('\t',1)[0], line if line.endswith('\n') else line+'\n'))
                    n += 1
        db.commit()
        with out.open('w') as fh:
            fh.write((header if header is not None else '\t'.join(EMPTY_DISCOVERY_HEADER)) + '\n')
            for line, in db.execute('SELECT line FROM rows ORDER BY gene, rowid'):
                fh.write(line)
        db.close()
    return n


def write_background(files: Sequence[Path], universe_files: Sequence[Path], out: Path, out_genes: Path) -> int:
    """background.output: `<gene>\\t<positions>` lines by gene (`NULL` for no position; `Gene\\tPosition` alone when
    there is no gene), the positions of the positional universe of each gene; and background_genes.output, the genes
    with at least one universe position, sorted, once each.

    `files` are the `*.b0.background` shards (the positions b_0 tested) and `universe_files` the `*.b0.universe` ones
    (those of them that pass the CAAS test under some orientation of the observed pairs, src.core.universe). The gene
    family is the genes that have a universe position, so the cleaned gene background (this list minus the genes the
    gene filter removes) has no gene without positions to test. A gene with none keeps its `NULL` line here and is not
    in background_genes.output."""
    from src.core.universe import parse_line, read_files
    tested = {}
    for path in files:
        for line in Path(path).read_text().splitlines():
            if not line.strip():
                continue
            gene, positions = parse_line(line)
            if gene in tested:
                raise ValueError(f'Duplicate gene in the tested positions: {gene}')
            tested[gene] = set(positions)
    universe = read_files(universe_files)
    if set(universe) != set(tested):
        raise ValueError('The positional universe and the tested positions cover different genes: '
                         f'{sorted(set(universe) ^ set(tested))[:10]}')
    for gene, positions in universe.items():
        if not set(positions) <= tested[gene]:
            raise ValueError(f'Universe position not tested by b_0 in {gene}')
    genes = sorted(universe)
    with out.open('w') as fh:
        if not genes:
            fh.write('Gene\tPosition\n')
        for gene in genes:
            fh.write(f"{gene}\t{','.join(map(str, universe[gene])) if universe[gene] else 'NULL'}\n")
    out_genes.write_text(''.join(gene + '\n' for gene in genes if universe[gene]))
    return len(genes)


def _na(value: Optional[str]) -> str:
    # readr convention: an empty cell is NA
    return "NA" if value is None or value == "" else value


def write_meta(discovery: Path, out_dir: Path) -> int:
    """meta_caas.tsv from discovery.tab, row by row in its order. Returns the number of rows."""
    from src.core.observed import is_strict_row
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(discovery, newline="") as src:
        reader = csv.DictReader(src, delimiter="\t", quoting=csv.QUOTE_NONE)
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not {'gene','position','caas','pattern'} <= set(fields):
            raise ValueError(f'Invalid discovery header: {discovery}')
        has_trait = 'trait' in fields
        cols = META_COLUMNS + (["trait"] if has_trait else []) + ["hyp_id"]
        with open(out_dir / "meta_caas.tsv", "w") as out:
            out.write("\t".join(cols) + "\n")
            for row in reader:
                from src.core.discovery import validate_row
                from src.utils.amino import normalize_amino_list
                validate_row(row)
                gene = row["gene"]
                trait = row.get("trait")
                pair = row.get("conserved_pair") or ""
                fields = [row_id(gene, row), f"{gene}_{row['position']}", gene, row["position"], _na(row.get("pattern")),
                          "TRUE" if not is_strict_row(row) else "FALSE",
                          re.sub(r"^\d+:", "", pair) if pair else "NA", _na(row.get("caas"))]
                parts = row['caas'].split('/')
                fields.append(str(len(set(normalize_amino_list(parts[0])) & set(normalize_amino_list(parts[1])))))
                if has_trait:
                    hyp = _HYP_LAST.sub(r"\1", trait) if trait and _HYP.search(trait) else trait
                    fields += [_na(trait), _na(hyp)]
                else:
                    fields.append("NA")
                out.write("\t".join(fields) + "\n")
                n += 1
    return n


def write_master(shards: Sequence[Path], out: Path, fallback_fields: Sequence[str], audit=None) -> int:
    """Merge sorted per-gene shards with bounded memory and exact schema/key checks."""
    fields = list(fallback_fields)
    def read(path):
        previous = None
        with gzip.open(path, 'rt', newline='') as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames != fields:
                raise ValueError(f'Master shard schema differs: {path}')
            for row in reader:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError(f'Malformed master row: {path}')
                gene, position = row['gene'], int(row['msa_pos'])
                if gene != path.name[:-len(MASTER_SUFFIX)]:
                    raise ValueError(f'Master shard must contain its named gene: {path}')
                key = (gene, position)
                if position < 0 or (previous is not None and key < previous):
                    raise ValueError(f'Unsorted/invalid master positions: {path}')
                if row['side'] not in {'top','bottom','none'}:
                    raise ValueError(f'Invalid master side: {path}')
                score = float(row['asr_path_score'])
                if not math.isfinite(score) or not 0 <= score <= 1:
                    raise ValueError(f'Invalid master score: {path}')
                previous = key
                yield gene, position, row

    def checked():
        gene = None
        last_position = None
        sides = set()
        count = 0
        digest = hashlib.sha256()
        def finish():
            if audit is not None and gene is not None:
                audit[gene] = dict(positions=count, sha256=digest.hexdigest())
        # One file per gene: opening one at a time avoids both genome-sized
        # lists and a file descriptor for every gene.
        merged = (row for path in sorted(shards, key=lambda p: p.name[:-len(MASTER_SUFFIX)]) for row in read(path))
        previous = None
        for g, pos, row in merged:
            if previous is not None and (g, pos) < previous:
                raise ValueError('Observed shards overlap or are not ordered')
            previous = (g, pos)
            if g != gene:
                finish()
                gene, last_position, count, digest = g, None, 0, hashlib.sha256()
            if pos != last_position:
                last_position, sides = pos, set()
                digest.update(f'{pos}\n'.encode())
                count += 1
            side = row['side']
            if side in sides or ('none' in sides) or (side == 'none' and sides):
                raise ValueError(f'Duplicate/conflicting master position/side: {g}:{pos}:{side}')
            sides.add(side)
            yield g, pos, row
        finish()
    return write_master_csv(checked(), out, fields, already_sorted=True)
