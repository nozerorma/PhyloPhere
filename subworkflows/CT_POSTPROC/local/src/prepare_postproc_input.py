#!/usr/bin/env python3
# prepare_postproc_input.py — Normalize the ct_disambiguation master table for CT post-processing.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Input preparation of CT post-processing: standardizes the key columns of the disambiguation
master table and, when the alignments and the contrast design are given, adds the extant-species
residue tally (residue_descriptors.py).

Called by:  CAAS_PREPARE_POSTPROC_INPUT process (ctpp_clustfilter.nf)
Inputs:     --input  disambiguation master CSV/TSV (the separator is detected from the header)
            --alignment, --alignment-format, --hyp-pairs  optional, for the tally
Outputs:    --output  normalized TSV (Gene, Position, ...)
            --removed-output  TSV of precluster removals (header only: no row is removed)
"""

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

from residue_descriptors import add_species_tally, _read_hyp_pairs
from postproc_io import load_observed, coverage


def _normalize_schema(df: pd.DataFrame) -> pd.DataFrame:
    # Only the structural keys are renamed (gene/msa_pos → Gene/Position, the names the
    # downstream steps read). The other columns (convergence_type, pvalue, caas)
    # keep their lowercase names.
    rename_map = {}
    if "gene" in df.columns:
        rename_map["gene"] = "Gene"
    if "msa_pos" in df.columns:
        rename_map["msa_pos"] = "Position"
    if rename_map:
        df = df.rename(columns=rename_map)

    return df


def _digest(path):
    import hashlib
    value = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1048576), b''):
            value.update(chunk)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare ct_disambiguation master CSV for CT post-processing."
    )
    parser.add_argument("--input", required=True, help="Disambiguation master CSV/TSV")
    parser.add_argument("--completion", required=True, help="observed_scoring.json")
    parser.add_argument(
        "--output",
        default="postproc_disambiguation_input.tsv",
        help="Normalized output TSV",
    )
    parser.add_argument(
        "--removed-output",
        default="removed_patterns_precluster.tsv",
        help="Precluster removal output TSV",
    )
    # Optional extant-species residue tally (top/bottom_species_residues, n_top/bottom_species,
    # n_conserved_pairs). It needs the alignment directory and the contrast design; without them
    # the columns stay empty.
    parser.add_argument("--alignment", default=None,
                        help="Alignment directory (recursive, gene = basename without alignment suffix; preserves internal '.')")
    parser.add_argument("--alignment-format", default="fasta",
                        help="Bio.AlignIO format for --alignment (default: fasta)")
    parser.add_argument("--hyp-pairs", default=None,
                        help="contrast_hypotheses_pairs.tsv (hypothesis_id, species1 = top, species2 = bottom)")
    parser.add_argument("--taxid-mapping", default=None)
    args = parser.parse_args()

    df, record = load_observed(args.input, args.completion)
    if args.hyp_pairs:
        pairs = _read_hyp_pairs(args.hyp_pairs)
        supplied = {h: [list(pair) for pair in items] for h, items in pairs.items()}
        for batch in record['batches']:
            if 'design' in batch:
                expected = {f'H{int(h)}': items for h, items in batch['design'].items()}
                if supplied != expected:
                    raise ValueError('Tally design differs from the original observed design')
    cleaned = _normalize_schema(df.copy())
    cleaned = add_species_tally(cleaned, args.alignment, args.hyp_pairs,
                                ali_format=args.alignment_format, taxid_mapping=args.taxid_mapping)
    # The quality universe contains every candidate, including side=none. ASR
    # categories/scores do not select it; the completion ledger checks coverage.
    Path('quality_universe.json').write_text(json.dumps({
        'version': 1, 'universe': 'all_ct_candidates',
        'alignments': {g:v for b in record['batches'] for g,v in b.get('alignments',{}).items()},
        'strict': bool(record['batches'] and record['batches'][0].get('strict',False)),
        'scientific': {k: record['batches'][0].get(k) for k in ('tree_sha256','model','posterior_threshold','taxid_sha256','weights_sha256')} if record['batches'] else None,
        'genes': coverage(cleaned, 'Gene', 'Position'),
        'excluded_genes': record.get('excluded_genes', []),
        'score_rule': 'caas_pos_species_pss_v4',
        'master_sha256': _digest(args.input),
        'completion_sha256': _digest(args.completion),
    }, indent=2))

    removed = pd.DataFrame(columns=[*df.columns, "removal_reason"])

    cleaned.to_csv(args.output, sep="\t", index=False)
    removed.to_csv(args.removed_output, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    print("Precluster MRCA posterior pruning retired; all input rows kept for postproc filtering")
    print(f"Rows kept for postproc filtering: {len(cleaned)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
