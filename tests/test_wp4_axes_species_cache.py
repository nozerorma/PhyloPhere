"""Cached domain scores must not borrow the species support of another hypothesis."""
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'subworkflows/CT_DISAMBIGUATION/local'))
from src.convergence import disambiguate_single as ds
from src.convergence.fop_pool import pool_domains


@pytest.mark.parametrize('explicit_signature', [False, True])
def test_score_cache_keeps_current_hypothesis_species(monkeypatch, explicit_signature):
    calls = []

    def scorer(pair_details, per_node_dist, node_index):
        calls.append(pair_details)
        pair = pair_details[0]
        return {
            'top': {'domain_scores': {1: 0.6}, 'domain_der': {1: 'D'}, 'domain_anc': {1: 'A'},
                    'domain_der_species': {1: pair['top_species'][0]}},
            'bottom': {'domain_scores': {}, 'domain_der': {}, 'domain_anc': {}, 'domain_der_species': {}},
            'domain_meta': {1: {'mrca_id': 7, 'state': 'A', 'posterior': 1.0}},
        }

    monkeypatch.setattr(ds, 'compute_domain_scores', scorer)
    cache, site_cache = {}, {}
    signature = ((1, 7, 'A', 'D', 'A'),) if explicit_signature else None

    def run(top, bottom, use_cache=True):
        pair = {'pair_id': 1, 'node_id': 7, 'focal_state': 'A', 'top_tip_mode': 'D',
                'bottom_tip_mode': 'A', 'top_species': [top], 'bottom_species': [bottom]}
        return ds._position_axes(SimpleNamespace(position_one_based=10), None, {}, {}, [pair],
                                 site_cache, cache if use_cache else None, signature)

    first = run('species_A', 'species_C')
    second = run('species_B', 'species_D')
    assert len(calls) == 1  # scores are still reused for identical ancestral and tip states
    assert first['top']['domain_der_species'] == {1: 'species_A'}
    assert second['top']['domain_der_species'] == {1: 'species_B'}
    assert second == run('species_B', 'species_D', use_cache=False)
    assert second['top']['domain_scores'] == first['top']['domain_scores']

    design = {('H1', 1): ('species_A', 'species_C', 1.0),
              ('H2', 1): ('species_B', 'species_D', 1.0)}
    pooled = pool_domains([{'hyp': 'H1', 'sides': first}, {'hyp': 'H2', 'sides': second}], design)
    assert pooled['top']['n_sp'] == 2
    assert pooled['top']['asr_path_score'] == pytest.approx(0.6)
