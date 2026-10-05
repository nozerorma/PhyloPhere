"""analyze_gene_disambiguation scores the entries and the contrasts it is given: no metadata file, trait file, diagnostics dump or ASR mode."""
import inspect
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.convergence import disambiguate_single as ds  # noqa: E402

RETIRED_PARAMETERS = ("caas_positions", "caas_metadata_path", "trait_file_path", "diagnostics_dir", "asr_mode")


def test_the_scorer_takes_none_of_the_retired_parameters():
    params = inspect.signature(ds.analyze_gene_disambiguation).parameters
    assert not set(RETIRED_PARAMETERS) & set(params), set(RETIRED_PARAMETERS) & set(params)


def test_the_entries_and_the_contrasts_are_required():
    params = inspect.signature(ds.analyze_gene_disambiguation).parameters
    for name in ("caas_entries", "trait_pairs"):
        assert params[name].default is inspect.Parameter.empty, name


@pytest.mark.parametrize("name", RETIRED_PARAMETERS)
def test_a_retired_parameter_is_rejected(name):
    with pytest.raises(TypeError, match=name):
        ds.analyze_gene_disambiguation("G", None, None, caas_entries=[], trait_pairs={}, **{name: None})


def test_the_historical_aliases_are_gone():
    assert not hasattr(ds, "analyze_gene_biochemistry") and not hasattr(ds, "analyze_caas_position_biochemistry")
