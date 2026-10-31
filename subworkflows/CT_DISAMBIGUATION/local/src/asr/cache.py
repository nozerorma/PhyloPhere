"""ASR cache provenance. Legacy/unverifiable entries are never reused."""
import hashlib
import json
from dataclasses import replace
from functools import lru_cache
from io import StringIO
from pathlib import Path
from Bio import Phylo


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


@lru_cache(maxsize=1)
def implementation():
    here = Path(__file__).parent
    return {name: digest(here / name) for name in
            ('asr_single.py', 'reconstruct.py', 'posterior.py', 'tree_parser.py', 'cache.py')}


def prepare_config(config, alignment_data, tree_data):
    model = str(config.model).lower()
    aliases = {'poison': 'poisson', 'jones': 'jtt'}
    model = aliases.get(model, model)
    from src.asr.reconstruct import MODEL_SPECS
    if model not in MODEL_SPECS:
        raise ValueError(f'Unsupported ASR model: {config.model}')
    buffer = StringIO()
    Phylo.write(tree_data.tree, buffer, 'newick', format_branch_length='%1.17g')
    rate = MODEL_SPECS[model]['aa_rate_file']
    rate_file = Path(__file__).parents[2] / 'dat' / rate if rate else None
    alignment_hash = hashlib.sha256()
    for record in alignment_data.alignment:
        alignment_hash.update(str(record.id).encode()); alignment_hash.update(b'\0')
        alignment_hash.update(str(record.seq).upper().encode()); alignment_hash.update(b'\0')
    provenance = dict(version=1, model=model, fix_blength=2,
                      alignment_sha256=alignment_hash.hexdigest(),
                      tree_sha256=hashlib.sha256(buffer.getvalue().encode()).hexdigest(),
                      mapping_sha256=hashlib.sha256(json.dumps(alignment_data.species_to_taxid, sort_keys=True).encode()).hexdigest(),
                      rate_sha256=digest(rate_file) if rate_file else None, implementation=implementation())
    key = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    return replace(config, model=model, provenance=provenance, output_dir=Path(config.output_dir) / 'validated_v1' / key)


def validate_entry(directory, expected):
    stamp = directory / 'provenance.json'
    if expected is None or not stamp.is_file():
        raise FileNotFoundError(f'ASR cache provenance unavailable: {directory}')
    try:
        saved = json.loads(stamp.read_text())
    except (ValueError, OSError) as exc:
        raise ValueError(f'Invalid ASR provenance: {directory}') from exc
    if saved.get('inputs') != json.loads(json.dumps(expected)):
        raise ValueError(f'ASR cache inputs differ: {directory}')
    artifacts = saved.get('artifacts', {})
    if not {'rst', 'rst1', 'tree_paml.nwk'} <= artifacts.keys():
        raise ValueError(f'Incomplete ASR entry: {directory}')
    for name, checksum in artifacts.items():
        if name not in {'rst', 'rst1', 'tree_paml.nwk'} or digest(directory / name) != checksum:
            raise ValueError(f'ASR artifact differs: {directory / name}')


def write_entry(directory, inputs):
    files = {name: digest(directory / name) for name in ('rst', 'rst1', 'tree_paml.nwk')}
    (directory / 'provenance.json').write_text(json.dumps(dict(inputs=inputs, artifacts=files), sort_keys=True, indent=2))
