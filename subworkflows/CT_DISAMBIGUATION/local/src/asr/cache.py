"""ASR cache lookup in the existing shared directory, without mandatory manifests."""
import hashlib
import re
from dataclasses import replace
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def prepare_config(config, alignment_data=None, tree_data=None):
    model = str(config.model).lower()
    model = {'poison': 'poisson', 'jones': 'jtt'}.get(model, model)
    from src.asr.reconstruct import MODEL_SPECS
    if model not in MODEL_SPECS:
        raise ValueError(f'Unsupported ASR model: {config.model}')
    return replace(config, model=model, output_dir=Path(config.output_dir))


def find_entry(root, gene):
    """Exact gene ID first, then producer .prot/reference-species suffixes.

    Existing flat and per-gene directory layouts are supported. Arbitrary dots
    and versions are preserved. No provenance sidecar is required.
    """
    root = Path(root)
    names = [gene]
    stem = gene
    if stem.endswith('.prot'):
        stem = stem[:-5]
        names.append(stem)
    base, sep, species = stem.rpartition('.')
    if sep and re.fullmatch(r'[A-Z][a-z]+_[a-z]+(?:_[a-z]+)*', species):
        names.append(base)
    for name in names:
        for directory in (root / f'asr_{name}', root / name / f'asr_{name}'):
            if all((directory / f).is_file() for f in ('rst', 'tree_paml.nwk')):
                return directory
    raise FileNotFoundError(f'ASR cache entry not found for {gene} in {root}')
