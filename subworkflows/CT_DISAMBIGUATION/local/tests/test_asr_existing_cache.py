"""Existing flat ASRs must be reused without a manifest or hashed directory."""
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.asr import asr_single
from src.asr.cache import find_entry, prepare_config


def entry(path):
    path.mkdir(parents=True)
    (path / 'rst').write_text('cached reconstruction')
    (path / 'tree_paml.nwk').write_text('(A:1,B:1);')
    return path


@pytest.mark.parametrize('gene', ['G.Homo_sapiens', 'G.Mus_musculus', 'G.Papio_anubis.prot'])
def test_existing_asr_does_not_launch_paml(tmp_path, monkeypatch, gene):
    saved = entry(tmp_path / 'asr_G')
    (saved / 'provenance.json').write_text('obsolete metadata is not required')
    config = asr_single.SingleGeneASRConfig(Path('unused.fa'), Path('unused.nwk'), output_dir=tmp_path)
    alignment = SimpleNamespace(alignment=SimpleNamespace(get_alignment_length=lambda: 2), taxid_to_species={})
    monkeypatch.setattr(asr_single, 'parse_paml_rst_node_level', lambda *a, **kw: ({3: {1: {'A': 1.0}, 2: {'A': 1.0}}}, {}))
    monkeypatch.setattr(asr_single, '_write_node_id_map', lambda *a: None)

    def forbidden(*a, **kw):
        raise AssertionError('Existing cache must not run PAML')

    monkeypatch.setattr(asr_single, '_run_asr_pipeline', forbidden)
    result = asr_single.run_asr_pipeline(gene, config, alignment_data=alignment, tree_data=object())
    assert result.rst_file == saved / 'rst'
    assert prepare_config(config).output_dir == tmp_path


def test_exact_complete_id_wins_and_internal_dots_stay_literal(tmp_path):
    entry(tmp_path / 'asr_AC1')
    with pytest.raises(FileNotFoundError):
        find_entry(tmp_path, 'AC1.2.Homo_sapiens')
    versioned = entry(tmp_path / 'asr_AC1.2')
    assert find_entry(tmp_path, 'AC1.2.Homo_sapiens') == versioned
    exact = entry(tmp_path / 'asr_AC1.2.Homo_sapiens')
    assert find_entry(tmp_path, 'AC1.2.Homo_sapiens') == exact


def test_incomplete_cache_is_not_returned(tmp_path):
    p = entry(tmp_path / 'asr_G')
    (p / 'tree_paml.nwk').unlink()
    with pytest.raises(FileNotFoundError):
        find_entry(tmp_path, 'G')


def test_nested_existing_layout(tmp_path):
    saved = entry(tmp_path / 'G' / 'asr_G')
    assert find_entry(tmp_path, 'G.Pan_troglodytes') == saved
