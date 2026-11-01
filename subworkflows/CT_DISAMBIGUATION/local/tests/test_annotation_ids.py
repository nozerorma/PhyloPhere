"""Annotation compatibility must preserve scoring IDs and observed/null lengths."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data.loaders import annotation_keys, load_ensembl_genes
from src.core.postproc import load_gene_lengths


@pytest.fixture
def references(tmp_path):
    annotation = tmp_path / 'genes.tsv'
    annotation.write_text('gene\tlength\nPPP3R2\t172\nAC1.2\t20\nOTHER\t10\n')
    maps = tmp_path / 'maps'
    maps.mkdir()
    (maps / 'PPP3R2.Homo_sapiens.map.tsv').write_text(
        'ori_codon_col\tstatus\tprot_ali_col\n'
        + ''.join(f'{i}\tselected\t{i}\n' for i in range(1, 176)))
    return annotation, maps


def test_annotation_whitelist_keeps_complete_alignment_ids(references):
    annotation, _ = references
    genes = ['PPP3R2.Homo_sapiens', 'AC1.2.Papio_anubis', 'OTHER.prot', 'UNKNOWN.Homo_sapiens']
    assert load_ensembl_genes(annotation, genes=genes) == set(genes[:-1])
    assert load_ensembl_genes(annotation) == {'PPP3R2', 'AC1.2', 'OTHER'}


def test_literal_dotted_ids_and_exact_annotation_priority():
    assert annotation_keys(['AC1.2', 'G.1'], {'AC1', 'G'}) == {}
    assert annotation_keys(['G.Homo_sapiens'], {'G', 'G.Homo_sapiens'}) == {
        'G.Homo_sapiens': 'G.Homo_sapiens'}


def test_ambiguous_annotation_aliases_fail():
    with pytest.raises(ValueError, match='Ambiguous annotated gene alias'):
        annotation_keys(['G.Homo_sapiens.prot'], {'G', 'G.Homo_sapiens'})


def test_observed_and_null_use_same_map_length_and_ids(references):
    annotation, maps = references
    gene = 'PPP3R2.Homo_sapiens'
    observed_sources, null_sources = {}, {}
    observed = load_gene_lengths(str(annotation), str(maps), genes=[gene], sources=observed_sources)
    null = load_gene_lengths(str(annotation), str(maps), genes=[gene, 'OTHER'], sources=null_sources)
    assert observed == {gene: 175.0}
    assert observed[gene] == null[gene]
    assert observed_sources[gene] == null_sources[gene] == 'untrimmed_alignment_columns'
    assert 'PPP3R2' not in observed and 'PPP3R2' not in null


def test_annotation_fallback_uses_literal_alignment_id(references):
    annotation, maps = references
    sources = {}
    assert load_gene_lengths(str(annotation), str(maps), genes=['AC1.2.Papio_anubis'], sources=sources) == {
        'AC1.2.Papio_anubis': 20.0}
    assert sources == {'AC1.2.Papio_anubis': 'annotation_units'}


def test_empty_observed_does_not_record_unused_annotation_lengths(references):
    annotation, maps = references
    sources = {}
    assert load_gene_lengths(str(annotation), str(maps), genes=[], sources=sources) == {}
    assert sources == {}
