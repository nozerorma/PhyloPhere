"""Local audit only: replay fixed inputs without running PAML or modifying production caches."""
import argparse
import csv
import json
import logging
import pickle
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--code', required=True, type=Path)
p.add_argument('--inputs', required=True, type=Path)
p.add_argument('--gene', default='A1BG')
p.add_argument('--context', required=True, type=Path)
p.add_argument('--build-cache', type=Path)
p.add_argument('--output', required=True, type=Path)
a = p.parse_args()
sys.path.insert(0, str(a.code.resolve()))
logging.basicConfig(level=logging.ERROR)
from src.core import driver
from src.core.labelings import read_trait_pairs, observed_pss
from src.core.observed import observed_entries, analyze_observed

gene = a.gene + '.Homo_sapiens'
inp = a.inputs.resolve()
if a.build_cache:
    from src.asr.asr_single import ASRResults
    from src.asr.posterior import parse_paml_rst_node_level
    cache = a.build_cache.resolve()
    posteriors, mapping = parse_paml_rst_node_level(cache/'rst', cache/'tree_paml.nwk',
                                                  threshold=0.1, return_mapping=True)
    def fixed_cache(_gene, config, alignment_data):
        return ASRResults(gene=gene, posteriors_site={}, posteriors_node=posteriors,
                          rst_file=cache/'rst', tree_file=cache/'tree_paml.nwk',
                          taxid_to_species=alignment_data.taxid_to_species, node_id_map=mapping), 'hit'
    driver._read_cached = fixed_cache
    def forbidden(*args, **kwargs):
        raise AssertionError('PAML execution is forbidden in this audit')
    driver.run_asr_pipeline = forbidden
    ctx = driver.load_gene_context(gene, str(inp), str(inp/'pruned_tree_file.nwk'),
                                  str(inp/'inputs/taxid/species_taxid_map.tsv'), 'lg',
                                  str(inp/'audit_lock_only'), 0.1,
                                  alignment_path=inp/(gene+'.fa'), alignment_format='fasta')
    a.context.write_bytes(pickle.dumps(ctx))
else:
    ctx = pickle.loads(a.context.read_bytes())

with (inp/(a.gene+'.discovery.tsv')).open() as f:
    raw = list(csv.DictReader(f, delimiter='\t'))
pairs = read_trait_pairs(inp/'traitfiles_ok_dir')
weights = observed_pss(inp/'traitfiles_ok_dir/contrast_hypotheses_pairs.tsv')
rows, diagnostics = analyze_observed(ctx, gene, observed_entries(gene,raw), pairs, weights, 0.1)
out = [{'position': r.position, 'side': r.side, 'score': r.asr_path_score,
        'type': r.convergence_type, 'hypotheses': r.n_hypotheses} for r in rows]
a.output.write_text(json.dumps({'rows': out, 'diagnostics': diagnostics}, default=str, indent=2))
print(a.output, len(out))
