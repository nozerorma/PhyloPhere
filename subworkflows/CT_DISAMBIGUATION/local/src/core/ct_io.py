"""Locate CT's canonical strict readers in a checkout or a staged Nextflow task."""
import importlib.util
from pathlib import Path


def _load():
    candidates = [Path(__file__).with_name('_ct_design_io.py'),
                  Path(__file__).resolve().parents[4] / 'CT/local/modules/design_io.py']
    for path in candidates:
        if path.is_file():
            spec = importlib.util.spec_from_file_location('_phylophere_design_io', path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise ImportError('CT design_io.py is required; stage it as src/core/_ct_design_io.py')


strict_io = _load()
