"""core.contract: the observed contract files written from the b_0 slice.

Oracles: the scripts of CONCAT_DISCOVERY and CONCAT_BACKGROUND (run on the same files, in every arrival order), the
code of the pattern-annotation report (7.CAAS_pattern_annotation.Rmd, extracted with knitr::purl and run without
rendering) for the meta tables, and the frozen PEPC discovery.tab and master for the byte-exact reconstructions.
"""
import csv
import gzip
import io
import itertools
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
SRC = ROOT / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(HERE))
from src.core import contract  # noqa: E402
from src.core.labelings import design_max_pairs  # noqa: E402
from src.core.meta import caas_id  # noqa: E402
from src.reporting.disambiguation_writers import _generate_dynamic_fields  # noqa: E402
from test_concat import HEADER, _run, _script, F1, F2, F3, BG  # noqa: E402

GOLD = HERE / "golden/pepc_c4_complete"
RMD = ROOT / "subworkflows/CT_META_CAAS/local/7.CAAS_pattern_annotation.Rmd"
H19 = "\t".join(contract.EMPTY_DISCOVERY_HEADER)


def _batch_dirs(tmp_path, files):
    """One directory per file, as the batches leave them: {name: text} -> [dir]"""
    dirs = []
    for i, (name, text) in enumerate(files):
        d = tmp_path / f"batch_{i}"
        d.mkdir()
        (d / name).write_text(text)
        dirs.append(d)
    return dirs


# ── discovery.tab and background: the former concatenation processes are the oracle ───────────────

@pytest.mark.parametrize("order", list(itertools.permutations([("A.b0.discovery.tsv", F1), ("B.b0.discovery.tsv", F2), ("C.b0.discovery.tsv", F3)])))
def test_discovery_equals_the_former_concatenation_in_every_arrival_order(tmp_path, order):
    base = tmp_path / "in"
    base.mkdir()
    dirs = _batch_dirs(base, order)
    out = tmp_path / "discovery.tab"
    contract.write_discovery(contract.batch_files(dirs, contract.DISCOVERY_SUFFIX), out)
    work = _run("CONCAT_DISCOVERY", [t for _, t in order], tmp_path, "discovery")
    assert out.read_text() == (work / "discovery.tab").read_text() and out.read_text().count(HEADER) == 1


def test_an_empty_discovery_has_the_header_of_the_former_empty_table(tmp_path):
    out = tmp_path / "discovery.tab"
    assert contract.write_discovery([], out) == 0
    work = tmp_path / "w"
    work.mkdir()
    (work / "run.sh").write_text(_script("CONCAT_DISCOVERY"))
    assert subprocess.run(["bash", "run.sh"], cwd=work, capture_output=True).returncode == 0
    assert out.read_text() == (work / "discovery.tab").read_text() == H19 + "\n"


@pytest.mark.parametrize("order", list(itertools.permutations(list(BG.items()))))
def test_background_equals_the_former_concatenation_in_every_arrival_order(tmp_path, order):
    base = tmp_path / "in"
    base.mkdir()
    dirs = _batch_dirs(base, [(f"{k}.b0.background", v) for k, v in order])
    out, genes = tmp_path / "background.output", tmp_path / "background_genes.output"
    contract.write_background(contract.batch_files(dirs, contract.BACKGROUND_SUFFIX), out, genes)
    work = _run("CONCAT_BACKGROUND", [v for _, v in order], tmp_path, "background")
    assert out.read_text() == (work / "background.output").read_text()
    assert genes.read_text() == (work / "background_genes.output").read_text() == "SHANK2\nSPN\n"


def test_an_empty_background_is_the_header_and_no_genes(tmp_path):
    out, genes = tmp_path / "b.output", tmp_path / "g.output"
    contract.write_background([], out, genes)
    assert out.read_text() == "Gene\tPosition\n" and genes.read_text() == ""


# ── meta tables: the report's own code is the oracle ─────────────────────────────────────────────

def _r_ready():
    return shutil.which("Rscript") and subprocess.run(
        ["Rscript", "-e", "for (p in c('tidyverse','knitr','scales','data.table','rmarkdown')) library(p, character.only = TRUE)"],
        capture_output=True).returncode == 0


needs_r = pytest.mark.skipif(not _r_ready(), reason="R with tidyverse and knitr is not available")


def _rmd_meta(work, discovery_text):
    """Run the code of the pattern-annotation report on a discovery.tab and return the meta tables it writes."""
    work.mkdir()
    (work / "discovery.tab").write_text(discovery_text)
    (work / "background.output").write_text("Gene\tPosition\n")
    (work / "run.R").write_text(f'''
knitr::purl("{RMD}", output = "rmd7.R", documentation = 0, quiet = TRUE)
code <- readLines("rmd7.R")
code <- code[!grepl("library\\\\(paletteer\\\\)", code)]
i <- grep("^filtered_discovery_path <- params", code)
code <- c(code[1:(i - 1)], 'params <- list(discovery_input = "discovery.tab", background_input = "background.output", output_dir = "out", caap_mode = TRUE, seed = "1998")', code[i:length(code)])
paletteer_d <- function(...) setNames(rep("#000000", 10), NULL)   # colors only: not used by the export
pdf(NULL)
eval(parse(text = code))
''')
    p = subprocess.run(["Rscript", "run.R"], cwd=work, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-1500:]
    return {f.name: f.read_text() for f in (work / "meta_caas").glob("*.tsv")}


def _compare_meta(tmp_path, discovery_text):
    ref = _rmd_meta(tmp_path / "rmd", discovery_text)
    (tmp_path / "d.tab").write_text(discovery_text)
    counts = contract.write_meta(tmp_path / "d.tab", tmp_path / "meta")
    got = {f.name: f.read_text() for f in (tmp_path / "meta").glob("*.tsv")}
    assert set(got) == set(ref) and counts["global"] == len(ref["global_meta_caas.tsv"].splitlines()) - 1
    rows = list(csv.DictReader(io.StringIO(discovery_text), delimiter="\t", quoting=csv.QUOTE_NONE))
    for name in ref:
        a, b = got[name].splitlines(), ref[name].splitlines()
        assert a[0] == b[0] and len(a) == len(b), name
        # every column but the id is the report's own text; the id is the content id of the row
        assert [l.split("\t")[1:] for l in a] == [l.split("\t")[1:] for l in b], name
    grp_rows = {}
    for r in rows:
        grp_rows.setdefault(r["caap_group"], []).append(r)
    for name, lines in got.items():
        sel = rows if name == "global_meta_caas.tsv" else grp_rows[name.split("_meta_caas")[0]]
        expect = [caas_id(r["gene"], r["position"], r["trait"], r["caap_group"], r["caas"], r["amino_encoded"], r["pattern"]) for r in sel]
        assert [l.split("\t")[0] for l in lines.splitlines()[1:]] == expect, name


@needs_r
def test_meta_tables_equal_those_of_the_report_on_the_frozen_pepc_discovery(tmp_path):
    _compare_meta(tmp_path, gzip.open(GOLD / "discovery.tab.gz", "rt").read())


@needs_r
def test_meta_tables_equal_those_of_the_report_with_missing_cells_and_unusual_traits(tmp_path):
    cols = contract.EMPTY_DISCOVERY_HEADER
    def row(gene, grp, trait, pos, caas, amino, pat, cons, pair):
        d = dict.fromkeys(cols, "1")
        d.update(gene=gene, mode="CAAP", caap_group=grp, trait=trait, position=str(pos), caas=caas, amino_encoded=amino, pattern=str(pat),
                 is_conserved_meta=cons, conserved_pair=pair, ms="")
        return "\t".join(d[c] for c in cols)
    rows = [row("G2", "US", "traitfile_H3.tab", 7, "AB/C", "AB/C", 1, "TRUE", "2:1,2"),
            row("G1", "GS1", "traitfile.tab", 0, "", "ab/c", 2, "FALSE", ""),
            row("G1", "US", "x_H1_H22.tab", 12, "D/E", "", 3, "", "10:4"),
            row("G1", "US", "traitfile_H10.tab", 12, "D/E", "D/E", 3, "true", "")]
    _compare_meta(tmp_path, "\n".join([H19] + rows) + "\n")


def test_an_empty_discovery_gives_header_only_meta_tables(tmp_path):
    (tmp_path / "d.tab").write_text(H19 + "\n")
    counts = contract.write_meta(tmp_path / "d.tab", tmp_path / "meta")
    assert counts == {"global": 0} and [f.name for f in (tmp_path / "meta").iterdir()] == ["global_meta_caas.tsv"]
    assert (tmp_path / "meta/global_meta_caas.tsv").read_text().split("\n")[0].endswith("\ttrait\thyp_id")


# ── master ───────────────────────────────────────────────────────────────────────────────────────

def _shard(path, text):
    with gzip.open(path, "wt", newline="") as fh:
        fh.write(text)


GOLD_MASTER = (GOLD / "caas_convergence_master.csv").read_text()


def test_the_master_of_one_gene_is_reproduced_byte_for_byte(tmp_path):
    _shard(tmp_path / "PEPC.master.csv.gz", GOLD_MASTER)
    assert contract.write_master([tmp_path / "PEPC.master.csv.gz"], tmp_path / "m.csv", []) == 217
    assert (tmp_path / "m.csv").read_text() == GOLD_MASTER


@pytest.mark.parametrize("reverse", [False, True])
def test_the_master_is_ordered_by_gene_whatever_the_shard_order(tmp_path, reverse):
    rows = list(csv.reader(io.StringIO(GOLD_MASTER)))
    other = io.StringIO()
    w = csv.writer(other, lineterminator="\r\n")
    w.writerow(rows[0])
    w.writerows([["PEPD"] + r[1:] for r in rows[1:]])
    _shard(tmp_path / "PEPC.master.csv.gz", GOLD_MASTER)
    _shard(tmp_path / "PEPD.master.csv.gz", other.getvalue())
    files = [tmp_path / "PEPD.master.csv.gz", tmp_path / "PEPC.master.csv.gz"]
    contract.write_master(files[::-1] if reverse else files, tmp_path / "m.csv", [])
    expect = GOLD_MASTER + other.getvalue().split("\r\n", 1)[1].replace("\r\n", "\n")  # read_text() normalizes line ends
    assert (tmp_path / "m.csv").read_text() == expect


def test_a_master_without_shards_has_the_columns_of_the_design(tmp_path):
    fields = _generate_dynamic_fields(4)
    assert contract.write_master([], tmp_path / "m.csv", fields) == 0
    assert (tmp_path / "m.csv").read_text().splitlines()[0] == GOLD_MASTER.splitlines()[0] == ",".join(fields)


# ── the command line ─────────────────────────────────────────────────────────────────────────────

def test_the_command_line_writes_every_contract_file_and_ignores_a_sentinel_among_the_directories(tmp_path):
    d = tmp_path / "b0"
    d.mkdir()
    (d / "PEPC.b0.discovery.tsv").write_text(gzip.open(GOLD / "discovery.tab.gz", "rt").read())
    (d / "PEPC.b0.background").write_text("PEPC\t1,2,3\n")
    _shard(d / "PEPC.master.csv.gz", GOLD_MASTER)
    (tmp_path / "NO_B0_OBSERVED").write_text("")
    design = tmp_path / "design"
    design.mkdir()
    for h in range(1, 3):  # a four-pair design: the master schema comes from here
        (design / f"traitfile_H{h}.tab").write_text("".join(f"s{i}a\t1\t{i}\ns{i}b\t0\t{i}\n" for i in range(1, 5)))
    assert design_max_pairs(design) == 4
    p = subprocess.run([sys.executable, str(SRC / "contract_main.py"), "--b0-dirs", str(tmp_path / "NO_B0_OBSERVED"), str(d),
                        "--design", str(design), "--output-dir", str(tmp_path / "out")], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-1500:]
    out = tmp_path / "out"
    assert (out / "discovery.tab").read_text() == gzip.open(GOLD / "discovery.tab.gz", "rt").read()
    assert (out / "caas_convergence_master.csv").read_text() == GOLD_MASTER
    assert (out / "background.output").read_text() == "PEPC\t1,2,3\n" and (out / "background_genes.output").read_text() == "PEPC\n"
    assert sorted(f.name for f in (out / "meta_caas").iterdir()) == ["GS1_meta_caas.tsv", "GS2_meta_caas.tsv", "GS3_meta_caas.tsv",
                                                                     "GS4_meta_caas.tsv", "US_meta_caas.tsv", "global_meta_caas.tsv"]


def test_the_command_line_with_no_directory_writes_the_empty_tables(tmp_path):
    design = tmp_path / "design"
    design.mkdir()
    (design / "traitfile_H1.tab").write_text("a\t1\t1\nb\t0\t1\n")
    p = subprocess.run([sys.executable, str(SRC / "contract_main.py"), "--b0-dirs", str(tmp_path / "NO_B0_OBSERVED"),
                        "--design", str(design), "--output-dir", str(tmp_path / "out")], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-1500:]
    out = tmp_path / "out"
    assert (out / "discovery.tab").read_text() == H19 + "\n" and (out / "background.output").read_text() == "Gene\tPosition\n"
    assert (out / "background_genes.output").read_text() == ""
    assert (out / "caas_convergence_master.csv").read_text().strip() == ",".join(_generate_dynamic_fields(1))
