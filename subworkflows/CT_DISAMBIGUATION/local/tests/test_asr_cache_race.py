"""Concurrent cache publication must preserve verified ASR entries and failures."""
import errno
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.asr.cache import find_entry
from src.core import driver


@pytest.fixture
def cache(tmp_path, monkeypatch):
    config = SimpleNamespace(output_dir=tmp_path)
    target = tmp_path / 'asr_G'
    result = SimpleNamespace(posteriors_node={1: {1: {'A': 1.0}}})

    def load(gene, used, alignment):
        entry = find_entry(used.output_dir, gene)
        if (entry / 'rst').read_text() != 'complete reconstruction\n':
            raise ValueError('Unreadable reconstruction')
        return result

    monkeypatch.setattr(driver, 'load_precomputed_asr', load)
    monkeypatch.setattr('time.sleep', lambda delay: None)
    return config, target, result


def populate(path):
    path.mkdir()
    for name in ('rst', 'rst1', 'tree_paml.nwk'):
        (path / name).write_text('complete reconstruction\n')


def test_real_nonempty_directory_race_reuses_valid_winner(cache, monkeypatch):
    config, target, result = cache

    def compute(*args, **kwargs):
        source = config.output_dir / 'stage'
        populate(source)
        assert not target.exists()
        populate(target)  # another writer wins after exists()
        os.replace(source, target)

    monkeypatch.setattr(driver, 'run_asr_pipeline', compute)
    assert driver._compute_or_read_asr('G', config, None, None) is result
    assert find_entry(config.output_dir, 'G') == target


@pytest.mark.parametrize('damage', ['unreadable', 'missing_rst', 'absent'])
def test_unverifiable_winner_does_not_hide_publication_error(cache, monkeypatch, damage):
    config, target, _ = cache
    if damage != 'absent':
        populate(target)
        if damage == 'unreadable':
            (target / 'rst').write_text('partial')
        else:
            (target / 'rst').unlink()
    error = FileExistsError(errno.EEXIST, 'concurrent publication')

    def compute(*args, **kwargs):
        raise error

    monkeypatch.setattr(driver, 'run_asr_pipeline', compute)
    with pytest.raises(FileExistsError) as caught:
        driver._compute_or_read_asr('G', config, None, None)
    assert caught.value is error


def test_delayed_winner_visibility(cache, monkeypatch):
    config, target, result = cache
    populate(target)
    read = driver._read_cached
    attempts = []

    def delayed(*args):
        attempts.append(1)
        return (None, 'miss') if len(attempts) == 1 else read(*args)

    def compute(*args, **kwargs):
        raise FileNotFoundError(errno.ENOENT, 'concurrent rename')

    monkeypatch.setattr(driver, '_read_cached', delayed)
    monkeypatch.setattr(driver, 'run_asr_pipeline', compute)
    assert driver._compute_or_read_asr('G', config, None, None) is result
    assert len(attempts) == 2


@pytest.mark.parametrize('error', [PermissionError(errno.EACCES, 'denied'), RuntimeError('PAML failed')])
def test_other_errors_are_not_hidden_even_with_valid_cache(cache, monkeypatch, error):
    config, target, _ = cache
    populate(target)

    def compute(*args, **kwargs):
        raise error

    monkeypatch.setattr(driver, 'run_asr_pipeline', compute)
    with pytest.raises(type(error)) as caught:
        driver._compute_or_read_asr('G', config, None, None)
    assert caught.value is error
