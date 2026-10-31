"""Regression checks for NAME_CURATION identifiers, validation and Nextflow staging."""

import csv
import itertools
import shutil
import subprocess
import sys
from pathlib import Path

import dendropy
import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "subworkflows/TRAIT_ANALYSIS/local/src"
sys.path.insert(0, str(SOURCE))

from alignment_species import alignment_species  # noqa: E402
from tree_cleanup import curate_traits, load_ali_sp_names, load_tax_id_map  # noqa: E402


@pytest.fixture
def inputs(tmp_path):
    (tmp_path / "tree.nwk").write_text("((A_sp:1,B_sp:2):3,C_sp:4);\n")
    (tmp_path / "names.txt").write_text("A_sp\nB_sp\nC_sp\n")
    (tmp_path / "tax.tsv").write_text("tax_id\tspecies\tfamily\n1\tA_sp\tFam\n2\tB_sp\tFam\n3\tC_sp\tFam\n")
    (tmp_path / "traits.csv").write_text("species,trait\nA_sp,1\nB_sp,2\nC_sp,3\n")
    return tmp_path


def run_cleanup(directory, *, bundle=True, traits=True, names=True, taxonomy=True):
    args = [sys.executable, str(SOURCE / "tree_cleanup.py"),
            "--tree", str(directory / "tree.nwk"),
            "--ali-sp-names", str(directory / "names.txt") if names else "NO_FILE",
            "--tax-id", str(directory / "tax.tsv") if taxonomy else "NO_FILE",
            "--output", str(directory / "out.nwk"), "--report", str(directory / "report.tsv")]
    if traits:
        args += ["--traits", str(directory / "traits.csv"), "--traits-out", str(directory / "out.csv")]
    if bundle:
        args += ["--species-table", str(directory / "species.tsv"),
                 "--species-report", str(directory / "species.txt"),
                 "--taxid-map", str(directory / "map.tsv")]
    return subprocess.run(args, capture_output=True, text=True, timeout=30)


@pytest.mark.parametrize("bundle", [False, True])
def test_cli_always_curates_requested_traits(inputs, bundle):
    result = run_cleanup(inputs, bundle=bundle)
    assert result.returncode == 0, result.stderr
    with (inputs / "out.csv").open() as fh:
        assert list(csv.reader(fh)) == [["species", "trait"], ["A_sp", "1"], ["B_sp", "2"], ["C_sp", "3"]]
    tree = dendropy.Tree.get(path=str(inputs / "out.nwk"), schema="newick", preserve_underscores=True)
    leaves = {n.taxon.label: n for n in tree.leaf_node_iter()}
    assert set(leaves) == {"A_sp", "B_sp", "C_sp"}
    assert leaves["A_sp"].distance_from_root() == 4
    assert leaves["B_sp"].distance_from_root() == 5
    assert leaves["C_sp"].distance_from_root() == 4


@pytest.mark.parametrize(("filename", "contents", "message"), [
    ("tree.nwk", "(('A sp':1,A_sp:2):3,C_sp:4);", "normalize to the same species"),
    ("names.txt", "Unknown\n", "no tree tips matched"),
    ("names.txt", "", "no alignment species names"),
    ("traits.csv", "", "header"),
    ("traits.csv", "species,trait\n", "no trait data rows"),
    ("traits.csv", "trait,species\n1\n", "row 2: expected 2 columns"),
    ("traits.csv", "species,trait\nA_sp,1,extra\n", "row 2: expected 2 columns"),
    ("traits.csv", "species,species\nA_sp,1\n", "duplicate column"),
    ("traits.csv", "wrong,trait\nA_sp,1\n", "column 'species' not found"),
    ("traits.csv", "species,tax_id,trait\nA_sp,bad,1\n", "positive integer"),
    ("traits.csv", "species,trait\nA_sp,1\nA_sp,999\n", "conflicting trait values"),
    ("tax.tsv", "tax_id\tspecies\tfamily\n1\tA_sp\n", "row 2: expected 3 columns"),
    ("tax.tsv", "id\tname\n1\tA_sp\n", "requires 'tax_id' and 'species'"),
    ("tax.tsv", "tax_id\tspecies\ninvalid\tA_sp\n", "positive integer"),
    ("tax.tsv", "tax_id\tspecies\n1\tA_sp\n2\tA_sp\n", "conflicting tax_ids"),
])
def test_bad_inputs_fail_before_writing_outputs(inputs, filename, contents, message):
    (inputs / filename).write_text(contents)
    result = run_cleanup(inputs)
    assert result.returncode != 0
    assert message in result.stderr
    assert "Traceback" not in result.stderr
    for name in ["out.nwk", "report.tsv", "out.csv", "species.tsv", "map.tsv"]:
        assert not (inputs / name).exists()


def test_missing_taxonomy_is_an_error_but_sentinel_is_optional(inputs):
    (inputs / "tax.tsv").unlink()
    failed = run_cleanup(inputs)
    assert failed.returncode != 0 and "tax.tsv" in failed.stderr
    assert run_cleanup(inputs, taxonomy=False).returncode == 0


def test_unmatched_trait_rows_are_reported_without_inventing_matches(inputs):
    (inputs / "traits.csv").write_text("species,trait\nUnknown,1\n")
    result = run_cleanup(inputs)
    assert result.returncode == 0, result.stderr
    assert (inputs / "out.csv").read_text().strip() == "species,trait"
    assert "Unknown\t\tremoved" in (inputs / "species.txt").read_text()


@pytest.mark.parametrize("reverse", [False, True])
def test_conflicting_synonym_values_are_rejected(reverse):
    rows = [["Alias", "1", "1"], ["A_sp", "1", "999"]]
    with pytest.raises(ValueError, match="conflicting trait values"):
        curate_traits(["species", "tax_id", "trait"], rows[::-1] if reverse else rows,
                      "species", {"A_sp"}, {"A_sp": "1", "Alias": "1"}, {"A_sp": ("1", "1", "")})


def test_ambiguous_trait_taxid_requires_canonical_name():
    mapping = {"A_sp": "1", "B_sp": "1"}
    resolved = {"A_sp": ("1", "1", ""), "B_sp": ("1", "2", "")}
    with pytest.raises(ValueError, match="ambiguous tax_id 1"):
        curate_traits(["species", "tax_id", "trait"], [["Alias", "1", "99"]],
                      "species", set(mapping), mapping, resolved)
    kept, _ = curate_traits(["species", "tax_id", "trait"], [["B_sp", "1", "99"]],
                            "species", set(mapping), mapping, resolved)
    assert kept == [["B_sp", "2", "99"]]


def test_ambiguous_tree_alias_is_rejected(inputs):
    (inputs / "tree.nwk").write_text("(Alias:1,C_sp:1);")
    (inputs / "tax.tsv").write_text("tax_id\tspecies\n1\tAlias\n1\tA_sp\n1\tB_sp\n3\tC_sp\n")
    result = run_cleanup(inputs)
    assert result.returncode != 0 and "ambiguous tax_id 1" in result.stderr
    assert not (inputs / "out.nwk").exists()


def test_identical_duplicates_collapse_to_exact_name(inputs):
    (inputs / "tax.tsv").write_text("tax_id\tspecies\n1\tAlias\n1\tA_sp\n")
    (inputs / "traits.csv").write_text("species,trait\nAlias,1\nA_sp,1\nA_sp,1\n")
    result = run_cleanup(inputs)
    assert result.returncode == 0, result.stderr
    assert (inputs / "out.csv").read_text().splitlines() == ["species,trait", "A_sp,1"]


def test_taxonomy_repeated_rows_and_txt_delimiter(inputs):
    path = inputs / "taxonomy.txt"
    path.write_text("tax_id\tspecies\tfamily\n1\tA_sp\t\n1\tA_sp\tFam\n1\tA_sp\tFam\n")
    assert load_tax_id_map(path) == ({"A_sp": "1"}, {"1": ["A_sp"]}, {"A_sp": "Fam"})


def test_precomputed_names_reject_normalization_collisions(inputs):
    path = inputs / "names.txt"
    path.write_text("A sp\nA_sp\n")
    with pytest.raises(ValueError, match="normalize to"):
        load_ali_sp_names(path)


ALIGNMENTS = {
    "fasta": ">A_sp description\nAA\n>B_sp description\nAC\n>C_sp\nCC\n",
    "phylip": "3 2\nA_sp      AA\nB_sp      AC\nC_sp      CC\n",
    "clustal": "CLUSTAL W\n\nA_sp AA\nB_sp AC\nC_sp CC\n\n",
    "stockholm": "# STOCKHOLM 1.0\nA_sp AA\nB_sp AC\nC_sp CC\n//\n",
    "nexus": "#NEXUS\nbegin data;\ndimensions ntax=3 nchar=2;\nformat datatype=protein;\nmatrix\nA_sp AA\nB_sp AC\nC_sp CC\n;\nend;\n",
}


@pytest.mark.parametrize("fmt", list(ALIGNMENTS))
def test_scan_follows_file_symlinks_and_uses_record_ids(tmp_path, fmt):
    directory = tmp_path / "alignments"
    directory.mkdir()
    source = tmp_path / "source.aln"
    source.write_text(ALIGNMENTS[fmt])
    (directory / "gene.aln").symlink_to(source)
    (directory / "notes.txt").write_text(">Not_a_species\n")
    assert alignment_species(directory, fmt) == ["A_sp", "B_sp", "C_sp"]


@pytest.mark.parametrize("kind", ["empty_directory", "empty_file", "broken_symlink", "malformed_alignment"])
def test_scan_reports_invalid_inputs(tmp_path, kind):
    if kind == "empty_file":
        (tmp_path / "gene.fa").write_text("")
    elif kind == "broken_symlink":
        (tmp_path / "gene.fa").symlink_to(tmp_path / "missing.fa")
    elif kind == "malformed_alignment":
        (tmp_path / "gene.phy").write_text("not an alignment")
    with pytest.raises(ValueError):
        alignment_species(tmp_path, "phylip" if kind == "malformed_alignment" else "fasta")


def nextflow_run(inputs, *, taxonomy, traits, names, fmt=None):
    if not shutil.which("nextflow"):
        pytest.skip("Nextflow is required for process-staging regression checks")
    (inputs / "subworkflows").symlink_to(ROOT / "subworkflows", target_is_directory=True)
    (inputs / "audit.nf").write_text(
        "include { NAME_CURATION } from './subworkflows/TRAIT_ANALYSIS/ta_name_curation.nf'\n"
        "workflow { NAME_CURATION(Channel.value(file(params.tree)), "
        "Channel.value(file(params.tax_id)), Channel.value(file(params.my_traits))) }\n")
    (inputs / "audit.config").write_text(
        "process.executor = 'local'\nprocess.shell = ['/bin/bash', '-euo', 'pipefail']\n"
        "params.sp_colname = 'species'\nparams.ali_sp_names = ''\nparams.alignment = ''\nparams.ali_format = 'fasta'\n")
    args = ["nextflow", "-log", str(inputs / "nextflow.log"), "run", str(inputs / "audit.nf"),
            "-c", str(inputs / "audit.config"), "-w", str(inputs / "work"), "-ansi-log", "false",
            "--tree", str(inputs / "tree.nwk"), "--outdir", str(inputs / "out"),
            "--tax_id", str(inputs / "tax.tsv") if taxonomy else "NO_FILE",
            "--my_traits", str(inputs / "traits.csv") if traits else "NO_FILE"]
    if names:
        args += ["--ali_sp_names", str(inputs / "names.txt")]
    if fmt:
        directory = inputs / "alignments"
        directory.mkdir()
        source = inputs / "source.aln"
        source.write_text(ALIGNMENTS[fmt])
        (directory / "gene.aln").symlink_to(source)
        args += ["--alignment", str(directory), "--ali_format", fmt]
    result = subprocess.run(args, cwd=inputs, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    output = inputs / "out/name_curation"
    tree = dendropy.Tree.get(path=str(output / "curated_tree.nwk"), schema="newick", preserve_underscores=True)
    assert {n.taxon.label for n in tree.leaf_node_iter()} == {"A_sp", "B_sp", "C_sp"}
    assert (output / "curated_traits.csv").exists() == traits
    assert (output / "species_table.tsv").exists()
    assert (output / "species_taxid_map.tsv").exists()


@pytest.mark.parametrize(("taxonomy", "traits", "names"), list(itertools.product([False, True], repeat=3)))
def test_nextflow_all_optional_input_combinations(inputs, taxonomy, traits, names):
    nextflow_run(inputs, taxonomy=taxonomy, traits=traits, names=names)


@pytest.mark.parametrize("fmt", ["fasta", "phylip"])
def test_nextflow_scan_with_missing_taxonomy_and_symlinks(inputs, fmt):
    nextflow_run(inputs, taxonomy=False, traits=True, names=False, fmt=fmt)
